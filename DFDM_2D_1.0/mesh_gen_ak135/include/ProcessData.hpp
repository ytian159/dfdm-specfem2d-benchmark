#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include "element.hpp"

namespace Mesh{
    class ProcessData{
        public:
            uint32_t process_id;
            uint32_t elements_count;
            std::vector<Mesh::Element> element_list;

            ProcessData(uint32_t idx);
            void print_data(std::string file_name);
            void read_data(std::string file_name);

    };
    
}