#pragma once
#include <mpi.h>
#include <vector>
#include "simulation.hpp"
#include "mpi_helper.hpp"
#include "ProcessData.hpp"
#include "element.hpp"
#include "spdlog/spdlog.h"
// Domain is composed of elements for each rank for e.g. each rank will have a domain of elements
// that it is responsible for.
namespace DFDM{
    // class simulation;
    // class Element2D;
    class Domain{
        public:
            std::shared_ptr<spdlog::logger> logger;
            std::vector<uint32_t> local_element_ids;
            std::vector<uint32_t> outer_elements, inner_elements;

            int total_ranks;
            int my_rank;
            // these are neighboring CPUs for this domain
            int x_neighbor_left; 
            int x_neighbor_right;
            int z_neighbor_top;
            int z_neighbor_bottom;

            Domain(DFDM::simulation& sim, std::string domain_file);
            Domain(DFDM::ProcessData& my_data, std::string domain_file, std::shared_ptr<spdlog::logger> logger_);
            void initialize_elements(DFDM::simulation& my_sim, std::vector<DFDM::Element2D>& local_elements, DFDM::ProcessData& my_data);
            void refine_grid(DFDM::simulation& my_sim, std::vector<DFDM::Element2D>& local_elements);
            bool is_local(const std::vector<uint32_t>& local_ids, int32_t neighbor_id);
            int rank_map_test(uint64_t elements, uint32_t ranks);
            void domain_dims(uint64_t tile_area, int& tile_dim_x, int& tile_dim_z);
            void domain_decomposition(uint64_t elements_x_dim, uint64_t elements_z_dim, uint32_t tile_dim_x, uint32_t tile_dim_z, uint32_t rank_id, std::vector<uint32_t>& domain_elements);
            void simulate_timestep(uint64_t step, const DFDM::simulation& sim, const DFDM::source& source, std::vector<DFDM::receiver>& receivers, const DFDM::Domain& curr_domain, std::vector<DFDM::Element2D>& local_elements);
            void domain_decomposition_rec(uint64_t elements_x_dim, uint64_t elements_z_dim, uint32_t tile_dim_x, uint32_t tile_dim_z, uint32_t rank_id, std::vector<uint32_t>& domain_elements);
    };
}
