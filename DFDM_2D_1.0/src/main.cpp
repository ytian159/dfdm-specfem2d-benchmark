#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include "CLI.hpp"
#include "simulation.hpp"
#include "mpi_helper.hpp"
#include "domain.hpp"
#include "ProcessData.hpp"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/udp_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

using namespace std;


int main(int argc, char **argv){
    DFDM::mpi_init();
    auto my_rank = DFDM::get_rank_id();
    /// logger init

    CLI::App app("DFDM 2D unstructured solver");
    bool print_debug = false; // default
    app.add_option("--print-debug", print_debug, "Set to print out debug information in log files"); 

    std::string log_file_dir = "logs/"; // default
    app.add_option("--log-files", log_file_dir, "log files directory"); // Check if the file exists
        std::string config_file = "../config/config.toml"; // default
    app.add_option("-c,--config-file", config_file, "Input Simulation Congurations file")->required(); // Check if the file exists
    // app.add_option("-d,--double", value, "Some Value");

    std::string output_dir; "./output/";
    app.add_option("-o,--output-directory", output_dir, "Directory to output all the simulation data")->required(); // Check if the file exists

    std::string domain_in_dir; "./domain_in/";
    app.add_option("-m,--mesh-input-directory", domain_in_dir, "Mesh input directory generated from mesher")->required(); // Check if the file exists

    std::string domain_out_dir; "./domain_out/";
    app.add_option("-d,--domain-output", domain_out_dir, "Directory where the solver will keep a copy of domain decomposition")->required(); // Check if the file exists

    bool tests_run =false;
    app.add_option("--enable-tests-run", tests_run, "When set to true, this will cause the solver to exit before starting the time step loop"); // Check if the file exists

    CLI11_PARSE(app, argc, argv);
 

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);
    console_sink->set_pattern("[DFDM::LOGGER::] [%^%l%$] %v");

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file_dir+"/multisink_" + std::to_string(my_rank) + ".txt", true);
    file_sink->set_level(print_debug ? spdlog::level::trace : spdlog::level::off);

    auto logger_mult = std::make_shared<spdlog::logger>("multi_sink");

    logger_mult->set_level(spdlog::level::trace);
    logger_mult->sinks().push_back(console_sink);
    logger_mult->sinks().push_back(file_sink);
    logger_mult->flush_on(spdlog::level::trace);
    /// logger init end

    if(my_rank == 0){
        if(print_debug){
            logger_mult->info("Debug mode is on, printing debug information to log files");
            logger_mult->info("Using log file directory: {}", log_file_dir);
        }else{
            logger_mult->info("Debug mode is off, not printing debug information to log files");
        }
    }
    DFDM::ProcessData my_data(my_rank, logger_mult);

    config_file = std::filesystem::absolute(config_file);
    output_dir = std::filesystem::absolute(output_dir);
    domain_in_dir = std::filesystem::absolute(domain_in_dir);
    domain_out_dir = std::filesystem::absolute(domain_out_dir);

    if(my_rank == 0){
        logger_mult->info("Using config file: {}", config_file);
        logger_mult->info("Using output directory: {}", output_dir);
        logger_mult->info("Using domain in directory: {}", domain_in_dir);
        logger_mult->info("Using domain out directory: {}", domain_out_dir);
    }


    double init_time = MPI_Wtime(); // start of initialization timer

    my_data.process_id = my_rank;
    my_data.read_data(domain_in_dir+"process_"+std::to_string(my_rank)+"_data.txt");
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(my_rank == 0){
        logger_mult->info("Mesh reading complete across all CPUs, total CPUs: {}", DFDM::get_total_ranks());
    }
    
    std::filesystem::create_directory(output_dir);
    std::filesystem::create_directory(domain_out_dir);
    #ifdef ENABLE_TEST
        std::filesystem::create_directory("./tests/generated_mats/");
    #endif
    if(my_rank == 0){
        logger_mult->debug("DFDM::Main: Output directories created");
    }
    
    auto total_ranks = DFDM::get_total_ranks();
    // reading simulation file, initializing simulation object
    DFDM::simulation my_sim(config_file,logger_mult);
    // std::cout << "mu: " << my_sim.model.mu[0] << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        logger_mult->info("DFDM::Main: Simulation initialized");
    }
    //updating list of local element (elements local to this process) ids.
    DFDM::Domain sim_domain(my_data, domain_out_dir, logger_mult);
    // MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        logger_mult->info("DFDM::Main: Simulation domain initialized");
    }
    //create elements local to this process
    std::vector<DFDM::Element2D> local_elements;
    if(my_rank == 0){
        logger_mult->info("DFDM::Main: Initializing local elements...");
    }   
    sim_domain.initialize_elements(my_sim, local_elements, my_data);
    // MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        logger_mult->info("DFDM::Main: Element initialization complete");
    }   
    DFDM::Element2D::update_dt_nt(local_elements, my_sim, logger_mult);
    if(my_rank == 0) logger_mult->info("DFDM::Main: Time step and number of time steps updated");
    // source init
    DFDM::source DFDM_source;
    // generate source only in the rank that has the source element
    if(DFDM::Element2D::is_local_element(my_sim.source_elem_id, local_elements)){
        decltype(auto) source_elem = DFDM::Element2D::get_local_element(my_sim.source_elem_id, local_elements);
        logger_mult->info("DFDM::Main: Source element found in rank:{}, source_elem:{}", my_rank, source_elem.global_id);
        logger_mult->debug("DFDM::Main: Source element Nx1: {}, Nz1: {}", source_elem.Nx1, source_elem.Nz1);
        DFDM_source.source_init(my_sim.source_elem_id, source_elem.Nx1, source_elem.Nz1, source_elem.order_b1, my_sim.time_steps, logger_mult);
        //generate source x in b1 basis, also gen source for z dim but with different knots, need to check this with Chao i.e. why do we initialize with element 0?
        DFDM_source.gen_source();
        //use the oeprators from source element here, need to
        DFDM_source.src_transform(source_elem.ops_x, source_elem.ops_z);
        DFDM_source.src_time_gen(my_sim.time_steps, my_sim.frequency, my_sim.delta_t);
        DFDM_source.invMsg_gen();
        DFDM_source.get_source_location(source_elem.element_grid);
        logger_mult->info("DFDM::Main: Source generated in rank:{}", my_rank);
    }

    // receivers initialization
    // generate receivers only in the ranks that have the receiver element
    std::vector<DFDM::receiver> DFDM_receivers;
    std::vector<uint32_t> local_receiver_ids;
    //initialize the receiver objects only in those ranks where the receiver elements are present
    for(uint32_t i = 0; i < my_sim.receiver_elem_ids.size(); i++){
        if(DFDM::Element2D::is_local_element(my_sim.receiver_elem_ids[i], local_elements)){
            local_receiver_ids.push_back(my_sim.receiver_elem_ids[i]);
            logger_mult->debug("DFDM::Main: Receiver element found in rank:{}, receiver_elem:{}", my_rank, my_sim.receiver_elem_ids[i]);
        }
    }
    
    logger_mult->debug("DFDM::Main: Number of receivers in rank {}: {}", my_rank, local_receiver_ids.size());

    // for each receiver
    for(uint32_t i = 0; i < local_receiver_ids.size(); i++){
        if(DFDM::Element2D::is_local_element(local_receiver_ids[i], local_elements)){
            decltype(auto) receiver_elem = DFDM::Element2D::get_local_element(local_receiver_ids[i], local_elements);
            logger_mult->debug("DFDM::Main: Receiver element found in rank:{}, receiver_elem:{}", my_rank, receiver_elem.global_id);
            DFDM::receiver DFDM_receiver;
            DFDM_receiver.receiver_init(local_receiver_ids[i], receiver_elem.Nx1, receiver_elem.Nz1, receiver_elem.order_b1, my_sim.time_steps, logger_mult);
            DFDM_receiver.gen_receiver(my_sim.delta_t);
            DFDM_receiver.get_receiver_location(receiver_elem.element_grid);
            logger_mult->info("DFDM::Main: Receiver location in rank:{}, receiver_elem:{}, xglobal:{}, zglobal:{}", my_rank, receiver_elem.global_id, DFDM_receiver.xglobal, DFDM_receiver.zglobal);

            DFDM_receivers.push_back(DFDM_receiver);
            logger_mult->debug("DFDM::Main: Receiver generated in rank:{}", my_rank);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // wait for all elements to be initialized across all ranks
    if(my_rank == 0){
        logger_mult->info("DFDM::Main: Simulation initialization complete");
    }
    init_time = MPI_Wtime() - init_time; // end of initialization timer

    double total_compute_time = 0;
    DFDM::utils::double_int local_max_di, local_min_di;
    //print out the grids:
    for(uint32_t z = 0; z < local_elements.size(); z++){
        local_elements[z].element_grid.x2d11.print_file(output_dir+"/grid_x_"+std::to_string(local_elements[z].global_id));
        local_elements[z].element_grid.z2d11.print_file(output_dir+"/grid_z_"+std::to_string(local_elements[z].global_id));
    }

    for(uint32_t z = 0; z < local_elements.size(); z++){
        local_elements[z].element_grid.mu11.print_file(output_dir+"/mu11_"+std::to_string(local_elements[z].global_id));
        local_elements[z].element_grid.mu22.print_file(output_dir+"/mu22_"+std::to_string(local_elements[z].global_id));
        local_elements[z].element_grid.rho12.print_file(output_dir+"/rho12_"+std::to_string(local_elements[z].global_id));
        local_elements[z].element_grid.rho21.print_file(output_dir+"/rho21_"+std::to_string(local_elements[z].global_id));
    }
    if(my_rank == 0)logger_mult->info("DFDM::Main: Simulation started");
    bool debug_run = false;
    uint32_t ndebug_steps = 1000;
    if(tests_run == true){
        return 1;
    }
    for(uint64_t steps = 0; steps < my_sim.time_steps; steps++){
        if(debug_run && steps > ndebug_steps){
            break;
        }
        double temp_time = MPI_Wtime();
        sim_domain.simulate_timestep(steps, my_sim, DFDM_source, DFDM_receivers, sim_domain, local_elements);

        total_compute_time += MPI_Wtime() - temp_time;

        // MPI_Barrier(MPI_COMM_WORLD);
        double local_max = -999999999.0;
        for(uint32_t z = 0; z < local_elements.size(); z++){
            auto max_temp = DFDM::utils::max2d(local_elements[z].state.Umid[steps]);
            if(max_temp > local_max){
                local_max = max_temp;
            }
        }
        


        if(steps%100 == 0){
            if(my_rank == 0){
                std::cout<<" PROCESS:" << my_rank << " Computing Max, Time Step:" << steps <<  std::endl;
            }
            std::vector<double> incoming_max;
            if(my_rank == 0){
                incoming_max.resize(total_ranks);
            }
            MPI_Gather(&local_max, 1, MPI_DOUBLE, incoming_max.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if(my_rank == 0){
                double g_max = *max_element(std::begin(incoming_max), std::end(incoming_max));
                std::cout<<" PROCESS:" << my_rank << " Time Step:" << steps << " G_MAX:" << g_max << std::endl;
            }  
        }

        if(steps%50 == 0){
            for(uint32_t z = 0; z < local_elements.size(); z++){
                local_elements[z].state.Umid[steps].print_file(output_dir+"/elem_"+std::to_string(local_elements[z].global_id)+"_"+std::to_string(steps)+".out");
                // local_elements[z].element_grid.x2d11.print_file(output_dir+"/grid_x_"+std::to_string(local_elements[z].global_id));
                // local_elements[z].element_grid.z2d11.print_file(output_dir+"/grid_z_"+std::to_string(local_elements[z].global_id));
            }
        }
    }
    local_max_di.val = total_compute_time;
    local_max_di.rank = my_rank;

    DFDM::utils::double_int max_c_time, min_c_time;
    double max_i_time = 0, min_i_time = 0, average_runtime = 0;
    MPI_Reduce(&local_max_di, &max_c_time, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_max_di, &min_c_time, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);    

    MPI_Reduce(&local_max_di.val, &average_runtime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);   

    average_runtime = average_runtime/total_ranks;

    MPI_Reduce(&init_time, &max_i_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&init_time, &min_i_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
        std::cout << "Max Time Spent in timestep loop:" << max_c_time.val << "\tFrom rank:"<< max_c_time.rank << "\tnMin Time spent in timestep loop:" << min_c_time.val << "\tFrom rank:"<< min_c_time.rank << std::endl;
        std::cout << "Average timestep loop runtime:" << average_runtime << std::endl;
        std::cout << "load balance factor for time step loop:" << max_c_time.val/average_runtime << std::endl;
        std::cout << "Max Init time:" << max_i_time << "\tMin Init time:" << min_i_time << std::endl;
    }

    // print out the receiver data
    logger_mult->debug("DFDM::Main: Printing receiver data");
    for (const auto& receiver : DFDM_receivers) {
        // print out the receiver data to file
        std::string filename = output_dir + "/receiver_elem_" + std::to_string(receiver.relement_id) + "_recorded_values.out";
        receiver.ur.print_file(filename);

        filename = output_dir + "/receiver_elem_" + std::to_string(receiver.relement_id) + "_stress.out";
        receiver.sigr.print_file(filename); // pressure is the mean of the two diagonal stress components
    }

    DFDM::mpi_finalize();
    return 0;
}
