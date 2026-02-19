#pragma once

#include <string>
#include <toml.hpp>
#include <iostream>
#include <vector>
#include "spdlog/spdlog.h"

namespace DFDM{
    class model{
        public:
            // uint64_t total_elements;
            // double velocity; // m/s
            std::shared_ptr<spdlog::logger> logger;
            std::string Model;
            int num_layers;
            std::vector<double> vp;
            std::vector<double> rho;
            std::vector<double> r;
            std::vector<double> mu;

            double R_EARTH = 6371000.0; // radius of the Earth in meters

            model();
            void model_init(std::string Model_, std::string model_file, std::shared_ptr<spdlog::logger> logger_);
            double get_value(std::string parameter,double x_global, double z_global);
        };
}