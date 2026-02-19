#include "simulation.hpp"
#include "model.hpp"

DFDM::simulation::simulation()
                            : frequency(1), duration(10)
                            , order_b1(5), gauss_order(5), delta_t(0.001), source_elem_id(0)
                            , receiver_elem_ids({}), time_steps(40000)
                            {}

DFDM::simulation::simulation(std::string toml_file, std::shared_ptr<spdlog::logger> logger){ // reads from toml file and process data (per CPU domain)
    auto data = toml::parse(toml_file);

    std::string Model  = toml::find<std::string>(data, "Model");
    std::string model_file = "";

    if(Model != "uniform" && Model != "ak135-f"){
        throw std::runtime_error("Model not supported, please use 'uniform' or 'ak135-f'");
    }else if(Model != "uniform"){
        model_file = toml::find<std::string>(data, "model_file");
    }else{
        model_file = ""; // no model file needed for uniform model
    }
    // DFDM::model model;
    if (Model.compare("uniform")==0){
        double vp  = toml::find<double>(data, "Vp");
        double vs = toml::find<double>(data, "Vs");
        double rho = toml::find<double>(data, "rho");
        model.Model = "uniform";
        model.vp = {vp};
        model.rho = {rho};
        model.mu = {rho*vp*vp};
        model.num_layers = 1;
        model.logger = logger;
        
        std::vector<double> rad = toml::find<std::vector<double>>(data, "rads");

        model.r = {*std::max_element(rad.begin(), rad.end())};

        int i=0;
        // std::cout << "Layer " << i << ", r=" << model.r[i] << " rho=" << model.rho[i] << " vp=" << model.vp[i] << " mu=" << model.mu[i] << std::endl;

    }
    else if (Model.compare("ak135-f")==0){
        // DFDM::model model;
        model.model_init("ak135-f", model_file, logger);
    }

    
    // grid_length_x = toml::find<double>(data, "grid_size_x");
    // grid_length_z = toml::find<double>(data, "grid_size_z");
    frequency = toml::find<double>(data, "frequency");
    duration = toml::find<double>(data, "duration");

    // order_b2 = toml::find<uint32_t>(data, "order_b2");
    order_b1 = toml::find<uint32_t>(data, "order_b1");

    gauss_order = toml::find<uint32_t>(data, "gauss_order");

    delta_t = toml::find<double>(data, "delta_t");
    time_steps = toml::find<uint64_t>(data, "time_steps");

    source_elem_id = toml::find<uint32_t>(data, "source_element");
    // receiver_elem_id = toml::find<uint32_t>(data, "receiver_element");
    receiver_elem_ids = toml::find<std::vector<uint32_t>>(data, "receiver_elements");

    auto my_rank = DFDM::get_rank_id();
    if(my_rank == 0){
        std::cout << "DFDM Simulation Configs:" << std::endl;
        // std::cout << "\tvelocity: "<<Vp<< std::endl;
        std::cout << "\tfrequency: "<< frequency<< std::endl;
        // std::cout << "\torder_b2: "<< order_b2<< std::endl;
        std::cout << "\torder_b1: "<< order_b1<< std::endl;
        std::cout << "\tgauss_order: "<<gauss_order<< std::endl;
        
        std::cout << "\treceiver_elements: ";
        for (const auto& id : receiver_elem_ids) {
            std::cout << id << " ";
        }
        std::cout << std::endl;
        // std::cout << "\trho: "<< rho<< std::endl;
        std::cout << "\tdelta_t: "<< delta_t<< std::endl;
        std::cout << "\ttime_steps: "<< time_steps<< std::endl;
        std::cout << "\tModel: "<< Model<< std::endl;
        std::cout << "\tModel file: "<< model_file<< std::endl;

    }
   

}




