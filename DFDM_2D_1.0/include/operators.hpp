#pragma once

#include <vector>
#include <cstdint>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

#include "matrix.hpp"
#include "utils.hpp"
#include "spdlog/logger.h"
#include <omp.h>
#include <set>

namespace DFDM{
    struct basis_param{
                uint32_t nb;
                uint32_t pb;
                std::vector<double> ps;
                std::vector<double> t;
        };

    enum OP_TYPE{X,Z};
    enum NBR{LEFT, RIGHT, LOWER, UPPER};
    class Operators{
        public:
                std::shared_ptr<spdlog::logger> logger;
                uint64_t element_id;
                OP_TYPE op_dim;
                uint64_t p1, p2, N1, N2, k1, k2;
                uint32_t ord_gi;
                DFDM::matrix<double> KK12; // N,N-1
                DFDM::matrix<double>  KK21; // N-1, N
                DFDM::matrix<double>  MM11; // N,N
                DFDM::matrix<double>  MM22; // N-1, N-1
                DFDM::matrix<double>  B1; //N,N
                DFDM::matrix<double>  B1_t;
                DFDM::matrix<double>  B2; //N,N
                DFDM::matrix<double>  B2_t;

                DFDM::matrix<double>  L11;
                DFDM::matrix<double>  L11_t;
                DFDM::matrix<double>  L22;
                DFDM::matrix<double>  L22_t;

                DFDM::matrix<double>  invL11;
                DFDM::matrix<double>  invL11_t;
                DFDM::matrix<double>  invL22;
                DFDM::matrix<double>  invL22_t;

                //Tranfor matrices, to be used when neighbors have mismatched points.
                //TAs are for bottom when x, and left when z
                //TBs are for top when x, and right when z
                DFDM::matrix<double> TA21;
                DFDM::matrix<double> TA12;
                DFDM::matrix<double> TA11;
                DFDM::matrix<double> TA22;
                DFDM::matrix<double> TB21;
                DFDM::matrix<double> TB12;
                DFDM::matrix<double> TB11;
                DFDM::matrix<double> TB22;

                //knot vectors
                std::vector<double> t1;
                std::vector<double> t2;

                //these can be made private later
                std::vector<std::vector<double>> xint; // gaussian int nodes for complete element
                std::vector<std::vector<double>> wxint; // gaussian int weights for complete element
                
                Operators(OP_TYPE op_dim_);
                void init_operators(uint64_t N_, uint64_t p_, uint32_t gauss_order, uint64_t global_id, std::shared_ptr<spdlog::logger> logger_);
                void gen_gauss_params();
                void gen_bspline();
                void mass_matrix(uint32_t option, DFDM::matrix<double>& mm);
                void stiffness_matrix(uint32_t kind, DFDM::matrix<double>& kk);
                void compute_operators();
                void gen_knots();
                std::vector<DFDM::matrix<double>> gen_connected_ops(int32_t neighbor_id, int32_t face_id, uint32_t nbr_Nx1, uint32_t nbr_Nz1, uint32_t nbr_px1, uint32_t nbr_pz1);
        private:
                void gaussian_int(uint32_t N, double a, double b, std::vector<double> &xints, std::vector<double> &wints); // compute weights
                double get_maxdiff(std::vector<double>& vec, std::vector<double>& yo);
                double dbspline_c(std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo, uint32_t d_order);
                double bspline_c(std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo);
                DFDM::basis_param setup_basis(uint32_t N, uint32_t p);
                DFDM::matrix<double> inner_product2(const DFDM::basis_param& basis1, const DFDM::basis_param& basis2);
        };

        // template<class T>
        // std::vector<std::vector<DFDM::matrix<double>>> Mat4D_prod_scalar(std::vector<std::vector<DFDM::matrix<double>>>& inmat_4d, T scalar);
        // std::vector<std::vector<DFDM::matrix<double>>> Mat4D_add(std::vector<std::vector<DFDM::matrix<double>>>& inmat_4d_1, std::vector<std::vector<DFDM::matrix<double>>>& inmat_4d_2);
        // std::vector<std::vector<DFDM::matrix<double>>> Mat4D_sub(std::vector<std::vector<DFDM::matrix<double>>>& inmat_4d_1, std::vector<std::vector<DFDM::matrix<double>>>& inmat_4d_2);
        // std::vector<std::vector<DFDM::matrix<double>>> Mat_Mat4D_prod(DFDM::matrix<double>& mat, std::vector<std::vector<DFDM::matrix<double>>>& mat_4d);
        // std::vector<std::vector<DFDM::matrix<double>>> transpose_4D(std::vector<std::vector<DFDM::matrix<double>>>& mat_4D);
        // void simulate_timesteps2D(simulation& curr_sim, source& sim_source_x, source& sim_source_z, mesh_grid& grid_x, mesh_grid& grid_z, 
        //                           mat_operators& ops_x, mat_operators& ops_z, boundary& boundary_x, boundary_t& boundary_z, DFDM::matrix<double>& invMsg);
}