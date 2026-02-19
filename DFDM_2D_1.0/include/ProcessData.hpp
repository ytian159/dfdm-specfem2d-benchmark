#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include "ElementData.hpp"

#include "spdlog/spdlog.h"
namespace DFDM{
    class ProcessData{
        public:
            uint32_t process_id;
            uint32_t elements_count;
            std::shared_ptr<spdlog::logger> logger;
            // these are neighboring CPUs for this domain/Process
            int x_neighbor_left; 
            int x_neighbor_right;
            int z_neighbor_top;
            int z_neighbor_bottom;

            std::vector<DFDM::ElementData> element_list;

            ProcessData(uint32_t idx, std::shared_ptr<spdlog::logger> logger);
            void print_data(std::string file_name);
            void read_data(std::string file_name);

    };
    
}