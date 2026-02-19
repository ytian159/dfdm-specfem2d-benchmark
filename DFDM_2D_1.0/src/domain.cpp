#include "domain.hpp"
#include "mpi_helper.hpp"

DFDM::Domain::Domain(DFDM::ProcessData& my_data, std::string domain_file, std::shared_ptr<spdlog::logger> logger_){
  logger = logger_;
  total_ranks = DFDM::get_total_ranks();
  my_rank = DFDM::get_rank_id();
  std::ofstream domain_out(domain_file+std::to_string(my_rank)+".out");
  logger->debug("initiating element copy to domain object");
  // copying element ids
  for (const auto& element : my_data.element_list) {
    local_element_ids.push_back(element.global_id);
  }
  logger->debug("element ids copied to domain object");

  // Write the local element ids to the domain output file
  domain_out << "RANK:" << my_rank << "\tNumber of local elements:" << local_element_ids.size() << std::endl;
  for (const auto& id : local_element_ids) {
    domain_out << "RANK:" << my_rank << "\telem:" << id << std::endl;
  }
  logger->debug("Element ids from process {} copied to domain output file:{}", my_rank, domain_file);
  domain_out.flush();
}

// DFDM::Domain::Domain(DFDM::simulation& sim, std::string domain_file){
//     total_ranks = DFDM::get_total_ranks();
//     my_rank = DFDM::get_rank_id();
//     std::ofstream domain_out(domain_file+std::to_string(my_rank)+".out");

//     auto rank_test = rank_map_test(sim.total_elements, total_ranks);
//     if(rank_test == -1){
//         std::cerr << "domain decomposition not possible with the given number of MPI ranks" << std::endl;
//         exit(-1);
//     }

//     tile_area = (sim.num_elements_x*sim.num_elements_z)/total_ranks;
//     domain_dims(tile_area, tile_dim_x, tile_dim_z);
//     auto larger_dim = (tile_dim_x > tile_dim_z) ? tile_dim_x : tile_dim_z;
//     if(std::abs(tile_dim_x - tile_dim_z) == larger_dim - 1 && (tile_dim_x - tile_dim_z != 0)){
//       if(my_rank == 0)
//         std::cout << "Using dedicated rectangular decomposition, this is not optimal, please try different rank combination. Square is optimal!" << std::endl;
//       domain_decomposition_rec(sim.num_elements_x, sim.num_elements_z, tile_dim_x, tile_dim_z, my_rank, local_element_ids);
//     }else{
//       domain_decomposition(sim.num_elements_x, sim.num_elements_z, tile_dim_x, tile_dim_z, my_rank, local_element_ids);
//     }


//     int tiles_along_z_dim = sim.num_elements_z/tile_dim_z;
//     int tiles_along_x_dim = sim.num_elements_x/tile_dim_x;

// // computing neighboring ranks in both directions
//     int32_t my_col = my_rank/tiles_along_x_dim;
//     int32_t my_row = my_rank%tiles_along_x_dim;
//     x_neighbor_left = my_rank - 1;
//     x_neighbor_right = my_rank + 1;
//     domain_out<<"RANK:"<<my_rank<<"\ttile_dim_x:"<<tile_dim_x<< " tile_dim_z:"<<tile_dim_z << std::endl;
//     for(uint32_t t = 0; t < local_element_ids.size(); t++){
//         domain_out << "RANK:" << my_rank << "\telem:" << local_element_ids[t]<< std::endl;
//     }
 
//     if((x_neighbor_right/tiles_along_x_dim) > my_col) x_neighbor_right = -1;
//     if((x_neighbor_left/tiles_along_x_dim) < my_col) x_neighbor_left = -1;


//     z_neighbor_top = my_rank + tiles_along_x_dim;
//     if(z_neighbor_top > (total_ranks - 1) || ((z_neighbor_top%tiles_along_x_dim) > my_row)) z_neighbor_top = -1;
//     z_neighbor_bottom  = my_rank - tiles_along_x_dim;
//     if(((z_neighbor_bottom%tiles_along_x_dim) < my_row)) z_neighbor_bottom = -1;
    
//     domain_out << "RANK:" << my_rank << "\trank_x_right:" << x_neighbor_right<< "\trank_x_left:" << x_neighbor_left<<"\trank_z_bottom:" << z_neighbor_bottom<<"\trank_z_top:" << z_neighbor_top<<std::endl;     domain_out.flush();
//     domain_out.close();
// }

bool DFDM::Domain::is_local(const std::vector<uint32_t>& local_ids, int32_t neighbor_id){
  if(neighbor_id >= 0){
    for(uint32_t t = 0; t < local_ids.size(); t++){
        if(local_ids[t] == neighbor_id){
            return true;
        }
      }
  }else{
    return true;
  }
  return false;
}

void DFDM::Domain::refine_grid(DFDM::simulation& my_sim, std::vector<DFDM::Element2D>& local_elements){
  for(auto& element : local_elements){
      element.grid_refine(my_sim);
  }
}

void DFDM::Domain::initialize_elements(DFDM::simulation& my_sim, std::vector<DFDM::Element2D>& local_elements, DFDM::ProcessData& my_data){
    logger->debug("Domain::initialize_elements:: Initializing elements list in rank {}", my_rank);
    for(uint32_t ele_count = 0; ele_count < local_element_ids.size(); ele_count++){
        local_elements.emplace_back(DFDM::Element2D(my_sim, my_data.element_list[ele_count], logger));
    }
    logger->debug("Domain::initialize_elements::elements list initialization in rank {} completed", my_rank);

    logger->debug("Domain::initialize_elements::Getting ready for elements setup loop...");
    for(uint32_t ele_count = 0; ele_count < local_elements.size(); ele_count++){
        logger->debug("Domain::initialize_elements::Setting up element {} in rank {}", my_data.element_list[ele_count].global_id, my_rank);
        local_elements[ele_count].initialize(my_sim, my_data.element_list[ele_count]); // this will read data from files and populate element
        logger->debug("Domain::initialize_elements::generating rotation matrices in element {} in rank {}",  my_data.element_list[ele_count].global_id, my_rank);
        local_elements[ele_count].gen_rotation_matrices();
        logger->debug("Domain::initialize_elements::generating dfdm matrices in element {} in rank {}",  my_data.element_list[ele_count].global_id, my_rank);
        local_elements[ele_count].gen_dfdm_matrices(); // compute dfdm matrices
        logger->debug("Domain::initialize_elements::Refining grid in element {} in rank {}",  my_data.element_list[ele_count].global_id, my_rank);
        local_elements[ele_count].grid_refine(my_sim);
        logger->debug("Domain::initialize_elements::Assigning locality to element {} in rank {}",  my_data.element_list[ele_count].global_id, my_rank);
        local_elements[ele_count].assign_locality(local_element_ids);
    }
    logger->debug("Domain::initialize_elements::Elements setup loop completed in rank {}", my_rank);

}

int DFDM::Domain::rank_map_test(uint64_t elements, uint32_t ranks){
    int n = (double)elements/ranks; 
    auto sqr = std::sqrt(n);
    if((n == sqr*sqr) || (elements%ranks == 0)){ // testing the below printed condition
      return n;
    }else{
        std::cerr << "ERROR: Ideally, number of MPI tasks (n) should chosen such that total number of elements (E) is completely divisble by n. Or, E/n should be a complete square." << std::endl;
        return -1;
    }
}

void DFDM::Domain::domain_dims(uint64_t tile_area, int& tile_dim_x, int& tile_dim_z){
  std::vector<int> dim_x, dim_z;
  for(int z = 1; z <= tile_area; z++){
    for(int x = 1; x <= tile_area; x++){
      if(x*z == tile_area){
        dim_x.push_back(x);
        dim_z.push_back(z);
      }
    }
  }

  int min = 10000;
  for(int t = 0; t < dim_x.size(); t++){
    if(min >= std::abs(dim_x[t] - dim_z[t])){
      min = std::abs(dim_x[t] - dim_z[t]);
      tile_dim_x = dim_x[t];
      tile_dim_z = dim_z[t];
    }
  }

}

// for rectangular domains
void DFDM::Domain::domain_decomposition_rec(uint64_t elements_x_dim, uint64_t elements_z_dim, uint32_t tile_dim_x, uint32_t tile_dim_z, uint32_t rank_id, std::vector<uint32_t>& domain_elements){

  for(uint32_t x = 0; x < tile_dim_x; x++){
    for(uint32_t z = 0; z < tile_dim_z; z++){
      uint32_t domains_z_dim = elements_z_dim / tile_dim_z; // number of domains in z dir of size tile
      uint32_t domains_x_dim = elements_x_dim / tile_dim_x;

      uint32_t d_row = rank_id % domains_x_dim; // which col and row the current domain is going to be, note that this row and col is based on number of domains along each dim.
      uint32_t d_col = rank_id / domains_x_dim;

      uint32_t domain_id = rank_id;

      int offset = 1;
      if(tile_dim_x > tile_dim_z){
        offset = (tile_dim_x > tile_dim_z) ? tile_dim_x : tile_dim_z;
      }
      uint32_t first_elem = (d_col * offset) + (d_row * offset);

      if(offset == 1) 
        offset = (tile_dim_x > tile_dim_z) ? tile_dim_x : tile_dim_z;
      else
        offset = 1;

      uint32_t curr_elem = first_elem + z*offset + x*offset;
      domain_elements.push_back(curr_elem);
    }
  }
}

// this works for square domains.
void DFDM::Domain::domain_decomposition(uint64_t elements_x_dim, uint64_t elements_z_dim, uint32_t tile_dim_x, uint32_t tile_dim_z, uint32_t rank_id, std::vector<uint32_t>& domain_elements){

  for(uint32_t x = 0; x < tile_dim_z; x++){
    for(uint32_t z = 0; z < tile_dim_x; z++){
      uint32_t domains_z_dim = elements_z_dim / tile_dim_x; // number of domains in z dir of size tile
      uint32_t domains_x_dim = elements_x_dim / tile_dim_z;

      uint32_t d_row = rank_id % domains_z_dim; // which col and row the current domain is going to be, note that this row and col is based on number of domains along each dim.
      uint32_t d_col = rank_id / domains_z_dim;

      uint32_t domain_id = d_col * domains_z_dim + d_row;
      uint32_t first_elem = (d_col * tile_dim_z) * elements_z_dim + (d_row * tile_dim_x);

      uint32_t curr_elem = first_elem + z + x*elements_z_dim;
      domain_elements.push_back(curr_elem);
    }
  }
}

void DFDM::Domain::simulate_timestep(uint64_t step, const DFDM::simulation& sim, const DFDM::source& source, std::vector<DFDM::receiver>& receivers, const DFDM::Domain& curr_domain, std::vector<DFDM::Element2D>& local_elements){
    logger->debug("Domain::simulate_timestep::Starting simulation timestep:{} in rank:{}", step, my_rank);
    // update wavefield
    logger->debug("Domain::simulate_timestep::Updating wavefield starting in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.update_wavefield(step, sim);
    }
    logger->debug("Domain::simulate_timestep::Wavefield updated for all elements in rank:{}", my_rank);


    // compute and copy local boundaries for U matrix
    logger->debug("Domain::simulate_timestep::Computing U boundary for all elements in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.compute_Uboundary_local(step);
    }
    logger->debug("Domain::simulate_timestep::U boundary computed for all elements in rank:{}", my_rank);

    // sends U boundary to neighbors, if neighbor is local, then copy to neighbor's U boundary
    // if neighor is remote, then send to remote through MPI
    logger->debug("Domain::simulate_timestep::Sending U boundary to neighbors from rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.send_Uboundary(local_elements, step);
    }
    // mpi_barrier();
    // wait for all sends to complete
    logger->debug("Domain::simulate_timestep::U boundary share completed to neighbors from rank:{}", my_rank);

    // receive U boundary from neighbors, if neighbor is local, then boundary has already been copied
    // if neighor is remote, then receive from remote through MPI
    logger->debug("Domain::simulate_timestep::Receiving U boundary from neighbors in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.receive_Uboundary(local_elements, step);
    }
    logger->debug("Domain::simulate_timestep::U boundary received from neighbors in rank:{}", my_rank);
    // mpi_barrier();

    logger->debug("Domain::simulate_timestep::applying U boundaries in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.apply_Uboundary(step);
    }
    logger->debug("Domain::simulate_timestep::U boundaries applied in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep:: computing KU2dA in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.compute_KU2dA(step);
    }
    logger->debug("Domain::simulate_timestep::KU2dA computed in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Computing and copying local boundaries for S matrix in rank:{}", my_rank);

    for(auto &loc : local_elements){
      loc.compute_Sboundary_local();
    }
    logger->debug("Domain::simulate_timestep:: completed Computing and copying local boundaries for S matrix in rank:{}", my_rank);

    for(auto &loc : local_elements){
      loc.Sboundary_rotate();
    }
    logger->debug("Domain::simulate_timestep:: Sboundary rotation completed in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Sending S boundary to neighbors from rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.send_Sboundary(local_elements);
    }
    // mpi_barrier();
    // wait for all sends to complete
    logger->debug("Domain::simulate_timestep::S boundary share completed to neighbors from rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Receiving S boundary from neighbors in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.receive_Sboundary(local_elements);
    }
    // mpi_barrier();
    // wait for all receives to complete
    logger->debug("Domain::simulate_timestep::S boundary received from neighbors in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::applying S boundaries in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.apply_Sboundary();
    }
    logger->debug("Domain::simulate_timestep::S boundaries applied in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep:: computing KS2dA in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.compute_KS2dA(step);
    }
    logger->debug("Domain::simulate_timestep:: KS2dA completed in rankm:{}", my_rank);

    logger->debug("Domain::simulate_timestep:: source injection starting in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.source_injection(source, sim, step);
    }
    logger->debug("Domain::simulate_timestep:: source injection completed in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Starting receiver record loop in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.receiver_record(receivers, sim, step);
    }
    logger->debug("Domain::simulate_timestep::Receiver record loop completed in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Starting wavefield record loop in rank:{}", my_rank);
    for(auto &loc : local_elements){
      loc.wavefield_record(step);
    }
    logger->debug("Domain::simulate_timestep::Wavefield record loop completed in rank:{}", my_rank);

    logger->debug("Domain::simulate_timestep::Time step:{} completed in rank:{}",step, my_rank);
    
}