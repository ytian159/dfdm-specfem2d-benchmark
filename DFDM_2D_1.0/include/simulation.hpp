#pragma once

#include <string>
#include <toml.hpp>
#include <iostream>
#include <vector>
#include "mpi_helper.hpp"
#include "model.hpp"
#include "spdlog/spdlog.h"

namespace DFDM{
    class simulation{
        public:
            // uint64_t total_elements;
            // double velocity; // m/s
            DFDM::model model;
            // double grid_length_x; // m
            // double grid_length_z; // m
            double frequency; // Hz
            double duration;
            // uint64_t num_elements_x; // elements in x dir
            // uint64_t num_elements_z; // elements in z dir
            // uint32_t order_b2; // order for B2 basis function (p)
            uint32_t order_b1; // order for B1 basis function (k = p+1)
            // uint64_t elem_size_x; // Nx, this will be per element now for unstructured mesh
            // uint64_t elem_size_z; // Nz, this will be per element now for unstructured mesh
            uint32_t gauss_order; // order of gaussian interpolation

            // double alpha_left;
            // double alpha_right;
            // double alpha_bottom;
            // double alpha_top;

            double delta_t;
            uint32_t source_elem_id;
            std::vector<uint32_t> receiver_elem_ids;
            uint64_t time_steps;


            simulation();
            simulation(std::string toml_file, std::shared_ptr<spdlog::logger> logger);
        };
}