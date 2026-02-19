#pragma once

#include "matrix.hpp"
#include <vector>
#include "operators.hpp"
#include "spdlog/spdlog.h"
#include "grid.hpp"

namespace DFDM{
    class source{
        public:
            std::shared_ptr<spdlog::logger> logger;
            uint32_t selement_id; //
            // uint64_t sloc; // (elem_size+1)/2  (middle of element)
            // double src_grid; //source value on grid x
            double xglobal;
            double zglobal;
            double xlocal;
            double zlocal;
            uint64_t Nx1;
            uint64_t Nx2;
            uint64_t Nz1;
            uint64_t Nz2;
            uint32_t px1;
            uint32_t px2;
            uint32_t pz1;
            uint32_t pz2;
            uint32_t kx1;
            uint32_t kx2;
            uint32_t kz1;
            uint32_t kz2;

            DFDM::matrix<double> gbx1; // sgb1
            DFDM::matrix<double> gbx2; // sgb2
            DFDM::matrix<double> gbz1; // sgb3
            DFDM::matrix<double> gbz2; // sgb4

            DFDM::matrix<double> gbx_inv1;
            DFDM::matrix<double> gbx_inv2;
            DFDM::matrix<double> gbz_inv1;
            DFDM::matrix<double> gbz_inv2; 

            std::vector<double> tx1;
            std::vector<double> tx2;
            std::vector<double> tz1;
            std::vector<double> tz2;
            // DFDM::matrix<double> inv_sourcegb1; // source converted to b1 basis

            DFDM::matrix<double> re;
            DFDM::matrix<double> Ft;

            DFDM::matrix<double> invMsg12;
            DFDM::matrix<double> invMsg21;
            source();
            void source_init(uint32_t global_id, uint32_t Nx1_, uint32_t Nz1_, uint32_t order_b1, uint64_t time_steps, std::shared_ptr<spdlog::logger> logger_);
            // double get_src_grid(DFDM::matrix<double>& grid);
            void gen_source();
            void src_transform(const DFDM::Operators& ops_x, const DFDM::Operators& ops_z);
            void src_time_gen(uint64_t time_steps, double frequency, double delta_t);
            void invMsg_gen();
            void get_source_location(const DFDM::Grid& element_grid);
    };
}