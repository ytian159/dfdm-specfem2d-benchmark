#include "ProcessData.hpp"
#include "ElementData.hpp"

DFDM::ProcessData::ProcessData(uint32_t idx, std::shared_ptr<spdlog::logger> logger_){
    logger = logger_;
    logger->debug("Logger initiated for rank:{}", idx);
    process_id = idx;
}

void DFDM::ProcessData::print_data(std::string file_name) {
    std::ofstream outfile(file_name);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file " << file_name << std::endl;
        return;
    }
    outfile << "Process_id:" << std::endl;
    outfile << process_id << std::endl;
    outfile << "Local_elements_count:" << std::endl;
    outfile << element_list.size() << std::endl;


    for (auto& element : element_list) {
        outfile << "Element_global_id:" << std::endl;
        outfile << element.global_id << std::endl;
        outfile << "Region_id:" << std::endl;
        outfile << element.region_id << std::endl;
        outfile << "Owner_CPU:" << std::endl;
        outfile << element.owner_cpu << std::endl;
        
        // Print neighbors
        outfile << "Element_neighbors (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.neighbors[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        //print face connections
        outfile << "Element_face_connections_nbrs (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.face_connections_nbrs[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_face_connections_self (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.face_connections_self[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;


        outfile << "Element_rotation_nbr (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.rotation_nbr[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_rotation_self (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.rotation_self[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_neighbor_Nx1 (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.neighbor_Nx1[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_neighbor_Nz1 (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.neighbor_Nz1[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_neighbor_processes (z_left, z_right, x_bottom, x_top):" << std::endl;
        for (size_t i = 0; i < 4; ++i) {
            outfile << element.neighbor_process[i];
            if (i < 3) {
                outfile << ",";
            }
        }
        outfile << std::endl;

        outfile << "Element_dimensions (NAX, NAZ):" << std::endl;
        outfile << element.NAX << "," << element.NAZ << std::endl;

        outfile << "Element_dimensions (nax, naz):" << std::endl;
        outfile << element.nax << "," << element.naz << std::endl;

        // Print xa and za values
        outfile << "xa_values, za_values:" << std::endl;
        for (uint32_t i = 0; i < element.nax; ++i) {
            for (uint32_t j = 0; j < element.naz; ++j) {
                outfile << element.xa(i, j) << "," << element.za(i, j) << std::endl;
            }
        }
    }
    
    outfile.close();
}

void DFDM::ProcessData::read_data(std::string file_name) {
    logger->debug("Rank {} attempting to open file {}", process_id, file_name);
    std::ifstream infile(file_name);
    if (!infile.is_open()) {
        logger->error("Unable to open file {}", file_name);
        throw std::runtime_error("Unable to open file " + file_name);
    }else{
        logger->debug("Rank {} has opened File {} for reading the mesh", process_id, file_name);
    }

    std::string line;
    std::getline(infile, line); // Skip "process_id:" line
    infile >> process_id;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line
    logger->debug("Process {} has read process id {}", process_id, process_id);
    // Read Local_elements_count
    std::getline(infile, line); // Skip "Local_elements_count:" line
    infile >> elements_count;
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line

    logger->debug("Process {} has {} elements", process_id, elements_count);

    element_list.clear(); // Clear existing elements

    while (std::getline(infile, line)) {
        if (line == "Element_global_id:") {
            DFDM::ElementData element;
            
            infile >> element.global_id;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::getline(infile, line); // Skip "Region_id:" line
            infile >> element.region_id;
            logger->debug("Element {} has region id {}", element.global_id, element.region_id);
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line

            std::getline(infile, line); // Skip "Domain_id:" line
            infile >> element.domain_id;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has domain id {}", element.global_id, element.domain_id);

            std::getline(infile, line); // Skip "Element_owner_cpu:" line
            infile >> element.owner_cpu;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has owner cpu {}", element.global_id, element.owner_cpu);

            std::getline(infile, line); // Skip "inflat:" line
            infile >> element.inflat;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has inflat {}", element.global_id, element.inflat);

            char comma;

            std::getline(infile, line); // Skip "theta:" line
            infile >> element.theta1 >> comma >> element.theta2;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has theta1 {} theta2 {}", element.global_id, element.theta1, element.theta2);

            std::getline(infile, line); // Skip "r:" line
            infile >> element.r1 >> comma >> element.r2;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has r1 {} r2 {}", element.global_id, element.r1, element.r2);

            std::getline(infile, line); // Skip "rx:" line
            infile >> element.rx1 >> comma >> element.rx2;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has rx1 {} rx2 {}", element.global_id, element.rx1, element.rx2);

            std::getline(infile, line); // Skip "rz:" line
            infile >> element.rz1 >> comma >> element.rz2;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has rz1 {} rz2 {}", element.global_id, element.rz1, element.rz2);

            std::getline(infile, line); // Skip "rad:" line
            infile >> element.rad0 >> comma >> element.rad1;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has rad0 {} rad1 {}", element.global_id, element.rad0, element.rad1);

            std::getline(infile, line); // Skip "coeff:" line
            infile >> element.coeff;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has coeff {}", element.global_id, element.coeff);

            std::getline(infile, line); // Skip "npart:" line
            infile >> element.npart;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has npart {}", element.global_id, element.npart);

            // Read neighbors
            std::getline(infile, line); // Skip "Element_neighbors" line
            infile >> element.neighbors[0] >> comma
                   >> element.neighbors[1] >> comma
                   >> element.neighbors[2] >> comma
                   >> element.neighbors[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger->debug("Element {} has neighbors {} {} {} {}", element.global_id, element.neighbors[0], element.neighbors[1], element.neighbors[2], element.neighbors[3]);
            //read face connections
            std::getline(infile, line); // Skip "Element_face_connections_nbrs" line
            infile >> element.face_connections_nbrs[0] >> comma
                   >> element.face_connections_nbrs[1] >> comma
                   >> element.face_connections_nbrs[2] >> comma
                   >> element.face_connections_nbrs[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            //Read face connections self
            std::getline(infile, line); // Skip "Element_face_connections_self" line
            infile >> element.face_connections_self[0] >> comma
                   >> element.face_connections_self[1] >> comma
                   >> element.face_connections_self[2] >> comma
                   >> element.face_connections_self[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


            // Read rotation
            std::getline(infile, line); // Skip "Element_rotation_nbr" line
            infile >> element.rotation_nbr[0] >> comma
                   >> element.rotation_nbr[1] >> comma
                   >> element.rotation_nbr[2] >> comma
                   >> element.rotation_nbr[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


            //Read rotation self
            std::getline(infile, line); // Skip "Element_rotation_self" line
            infile >> element.rotation_self[0] >> comma
                   >> element.rotation_self[1] >> comma
                   >> element.rotation_self[2] >> comma
                   >> element.rotation_self[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            
            //Read neighbor_Nx1
            std::getline(infile, line); // Skip "neighbor_Nx1:" line
            infile >> element.neighbor_Nx1[0] >> comma
                   >> element.neighbor_Nx1[1] >> comma
                   >> element.neighbor_Nx1[2] >> comma
                   >> element.neighbor_Nx1[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            //Read neighbor_Nz1
            std::getline(infile, line); // Skip "neighbor_Nz1:" line
            infile >> element.neighbor_Nz1[0] >> comma
                   >> element.neighbor_Nz1[1] >> comma
                   >> element.neighbor_Nz1[2] >> comma
                   >> element.neighbor_Nz1[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // Read neighbor processes
            std::getline(infile, line); // Skip "Element_neighbor_processes" line
            infile >> element.neighbor_process[0] >> comma
                   >> element.neighbor_process[1] >> comma
                   >> element.neighbor_process[2] >> comma
                   >> element.neighbor_process[3];
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            std::getline(infile, line); // Skip "NAX, NAZ:" line
            infile >> element.NAX >> comma >> element.NAZ;
            logger->debug("Element {} has refined dimensions {}x{}", element.global_id, element.NAX, element.NAZ);
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line
            // Read dimensions
            std::getline(infile, line); // Skip "Element_dimensions" line
            infile >> element.nax >> comma >> element.naz;
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            // Initialize xa and za matrices
            element.xa.resize(element.nax, element.naz);
            element.za.resize(element.nax, element.naz);
            logger->debug("Element {} has unrefined dimensions {}x{}", element.global_id, element.nax, element.naz);
            // Read xa and za values
            std::getline(infile, line); // Skip "xa_values, za_values:" line
            for (uint32_t i = 0; i < element.nax; ++i) {
                for (uint32_t j = 0; j < element.naz; ++j) {
                    double xa_val, za_val;
                    infile >> xa_val >> comma >> za_val;
                    element.xa.value_set(i, j, xa_val);
                    element.za.value_set(i, j, za_val);
                    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }
            }
            logger->debug("Element {} pushed", element.global_id);
            element_list.push_back(element);
        }
    }
    logger->debug("Process {} has read all elements", process_id);
    infile.close();
}


