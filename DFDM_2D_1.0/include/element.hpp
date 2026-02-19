// #pragma once
// #ifndef ELEMENT2D_HPP
// #define ELEMENT2D_HPP

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <string>

#include "matrix.hpp"
#include "grid.hpp"
#include "knots.hpp"
#include "operators.hpp"
#include "State.hpp"
#include "source.hpp"
#include "receiver.hpp"
#include "simulation.hpp"
#include "mpi_helper.hpp"
#include "utils.hpp"
#include "spdlog/spdlog.h"

namespace DFDM{
    // class simulation;
    class Element2D{
        public:
            std::shared_ptr<spdlog::logger> logger;
            uint32_t element_id;// for legacy code this is till used
            uint32_t global_id;
            uint32_t domain_id; // this is basically the domain id or local rank id
            uint32_t region_id; // this is basically the domain id or local rank id
            uint32_t owner_cpu;

            uint32_t nax; // size of unrefined grid
            uint32_t naz;

            uint32_t Nx1; // size of refined grid
            uint32_t Nz1;

            uint32_t Nx2;
            uint32_t Nz2;

            uint32_t order_b1; // number of basis functions p1?
            // uint32_t order_b2; // p2?

            uint64_t gauss_order;

            DFDM::Grid element_grid;

            DFDM::Operators ops_x;
            DFDM::Operators ops_z;

            DFDM::State state;
            uint32_t synch_needed;
            double alpha_mo, alpha_po, alpha_om, alpha_op;

            //connectivity lists
            int32_t neighbors[4]; // top, bottom, left ,right
            int32_t neighbor_process[4];// change this to left, right, bottom , top
            int32_t face_connections_nbrs[4]; // the connected faces in corresponding element (the face to be received from nbr
            int32_t face_connections_self[4]; // which of my faces are connected to which neighbor (face to be sent)
            bool rotation_nbr[4]; // 0 no rotation, 1 rotation, rotation of neighbor's face relative to my face
            bool rotation_self[4]; // 0 no rotation, 1 rotation, rotation of my face relative to neighbor's face
            uint32_t neighbor_Nx1[4];
            uint32_t neighbor_Nz1[4];
            uint32_t neighbor_px1[4];
            uint32_t neighbor_pz1[4];
            bool locality[4] = {0, 0, 0, 0};

            // rotation matrices
            DFDM::matrix<double> rot_mat_mo;
            DFDM::matrix<double> rot_mat_po;
            DFDM::matrix<double> rot_mat_om;
            DFDM::matrix<double> rot_mat_op;

            std::vector<MPI_Status> s_status;
            std::vector<MPI_Request> s_request;
            std::vector<MPI_Status> r_status;
            std::vector<MPI_Request> r_request;
       
            Element2D(DFDM::simulation& sim_in, DFDM::ElementData& el_data, std::shared_ptr<spdlog::logger> logger_);
            void initialize(const DFDM::simulation& sim, const DFDM::ElementData& el_data);
            void grid_refine(DFDM::simulation& sim);
            void gen_dfdm_matrices();
            bool is_receiver(const DFDM::simulation& sim);
            bool is_source(const DFDM::simulation& sim);

            void update_wavefield(uint64_t time_step, const DFDM::simulation& sim);
            void compute_Uboundary_local(uint64_t time_step);
            void send_Uboundary(std::vector<DFDM::Element2D>& loc_elements, uint64_t time_step);
            void receive_Uboundary(std::vector<DFDM::Element2D>& loc_elements, uint64_t time_step);
            void apply_Uboundary(uint64_t time_step);
            void compute_KU2dA(uint64_t time_step);
            void compute_Sboundary_local();
            void Sboundary_rotate();
            void compute_boundary_Svalue_inn_elem(DFDM::matrix<double>& xinvL_t, DFDM::matrix<double>& zinvL_t, 
                                                  DFDM::matrix<double>& S_, DFDM::matrix<double>& Smo, DFDM::matrix<double>& Spo, 
                                                  DFDM::matrix<double>& Som, DFDM::matrix<double>& Sop);
            void send_Sboundary(std::vector<DFDM::Element2D>& loc_elements);
            void receive_Sboundary(std::vector<DFDM::Element2D>& loc_elements);
            void apply_Sboundary();
            void compute_KS2dA(uint64_t time_step);
            void source_injection(const DFDM::source& source, const DFDM::simulation& sim, uint64_t it);
            void receiver_record(std::vector<DFDM::receiver>& receivers, const DFDM::simulation& sim, uint64_t it);
            void wavefield_record(uint64_t it);
            void rotate_sij(const DFDM::matrix<double>& sxx, const DFDM::matrix<double>& szz, 
                            DFDM::matrix<double>& sxx_r, DFDM::matrix<double>& szz_r, const DFDM::matrix<double>& rot_mat);
            DFDM::matrix<double> flip1D_boundary(DFDM::matrix<double>& face, uint32_t i1, uint32_t i2, uint32_t j1, uint32_t j2, bool orient);
            DFDM::matrix<double> get_face(DFDM::matrix<double>& Umo, DFDM::matrix<double>& Upo, DFDM::matrix<double>& Uom, DFDM::matrix<double>& Uop, uint32_t jFace, bool orient);

            void assign_locality(std::vector<uint32_t>& local_elements);
            DFDM::Element2D& getElement(uint32_t target_el, std::vector<DFDM::Element2D>& loc_elements);
            void gen_rotation_matrices();
            DFDM::receiver& get_local_receiver(uint64_t global_id, std::vector<DFDM::receiver>& receivers);

            static bool is_local_element(uint64_t global_id, std::vector<DFDM::Element2D>& elements);
            static DFDM::Element2D& get_local_element(uint64_t global_id, std::vector<DFDM::Element2D>& elements);
            static void update_dt_nt(const std::vector<DFDM::Element2D>& elements, DFDM::simulation sim, std::shared_ptr<spdlog::logger> logger_);
        private:

    };

}

// #endif // ELEMENT2D_HPP
