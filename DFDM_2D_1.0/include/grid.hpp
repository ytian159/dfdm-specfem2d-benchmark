#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <string>    
#include "matrix.hpp"
#include "ElementData.hpp"
#include "operators.hpp"
#include "simulation.hpp"
#include "spdlog/spdlog.h"
#include "model.hpp"

namespace DFDM{   
    class Grid{
        public:
            std::shared_ptr<spdlog::logger> logger;
            uint32_t element_id;
            double inflat;
            int region_id;
            double theta1, theta2;
            double r1,r2;
            double rx1,rx2;
            double rz1,rz2;
            double coeff;
            double rad0, rad1;
            int npart;

            uint64_t Nxt; // size of unrefined grid
            uint64_t Nzt;

            uint64_t Nx1;
            uint64_t Nz1;

            uint64_t Nx2;
            uint64_t Nz2;

            uint32_t px1, pz1;

            DFDM::matrix<double> x2d11t; // the unrefined grid vals (read from the mesh files)
            DFDM::matrix<double> z2d11t;                
            DFDM::matrix<double> x2d11; // refined grid vals
            DFDM::matrix<double> z2d11;

            DFDM::matrix<double> dxpdx11;
            DFDM::matrix<double> dzpdx11;
            DFDM::matrix<double> dxpdz11;
            DFDM::matrix<double> dzpdz11;
            DFDM::matrix<double> dxpdx12;
            DFDM::matrix<double> dzpdx12;
            DFDM::matrix<double> dxpdz12;
            DFDM::matrix<double> dzpdz12;
            DFDM::matrix<double> dxpdx21;
            DFDM::matrix<double> dzpdx21;
            DFDM::matrix<double> dxpdz21;
            DFDM::matrix<double> dzpdz21;
            DFDM::matrix<double> dxpdx22;
            DFDM::matrix<double> dzpdx22;
            DFDM::matrix<double> dxpdz22;
            DFDM::matrix<double> dzpdz22;
            DFDM::matrix<double> Jac12;
            DFDM::matrix<double> Jac21;
            DFDM::matrix<double> Jac11;
            DFDM::matrix<double> Jac22;

            DFDM::matrix<double> rho12;
            DFDM::matrix<double> rho21;
            DFDM::matrix<double> mu11;
            DFDM::matrix<double> mu22;

            void grid_init(const DFDM::ElementData& el_data, const DFDM::simulation sim, std::shared_ptr<spdlog::logger> logger_);
            void refine(int p, int pt, DFDM::model model, uint32_t ord_gi, DFDM::Operators& xops, DFDM::Operators& zops);
            void local_to_global(double x_local, double z_local, double& x_global, double& z_global);
            // double elem_start;
            // double elem_end;
            // // grid, initiated in a func, do we really need this?
            // DFDM::matrix<double> element_grid; // same as in matlab
            // void init_ends(uint64_t total_elements, double domain_len, uint64_t offset_ind);
            // void init_grid(uint64_t elem_size);
        private:
            std::vector<double> hatpointsx1;
            std::vector<double> hatpointsx2;
            std::vector<double> hatpointsz1;
            std::vector<double> hatpointsz2;
            void hat_points(std::vector<double>& hatpoints1, std::vector<double>& hatpoints2, DFDM::Operators& ops, uint32_t N, uint32_t p, uint32_t ord_gi);
    };
}