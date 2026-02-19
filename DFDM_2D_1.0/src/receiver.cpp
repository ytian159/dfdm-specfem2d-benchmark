#include "receiver.hpp"

DFDM::receiver::receiver(){
    
}
void DFDM::receiver::receiver_init(uint32_t global_id, uint32_t Nx1_, uint32_t Nz1_, uint32_t order_b1, uint64_t time_steps, std::shared_ptr<spdlog::logger> logger_){
    logger = logger_;
    relement_id = global_id; // 2nd element

    // local coordinates range from 0 to 1.0
    xlocal = 0.5;
    zlocal = 0.5;

    Nx1 = Nx1_;  
    Nx2 = Nx1 - 1;
    Nz1 = Nz1_;  
    Nz2 = Nz1 - 1;

    px1 = order_b1;   
    px2 = px1 - 1;
    pz1 = order_b1;   
    pz2 = pz1 - 1;

    kx1 = px1 + 1;  
    kz1 = pz1 + 1; 
    kx2 = px2 + 1;
    kz2 = pz2 + 1;

    tx1 = DFDM::utils::get_knot_vector(Nx1, kx1);
    tx2 = DFDM::utils::get_knot_vector(Nx2, kx2);
    tz1 = DFDM::utils::get_knot_vector(Nz1, kz1);
    tz2 = DFDM::utils::get_knot_vector(Nz2, kz2);

    gbx1.resize(1, Nx1);
    gbx2.resize(1, Nx2);
    gbz1.resize(1, Nz1);
    gbz2.resize(1, Nz2);

    ur.resize(time_steps, 3); // for benchmarking
    sigr.resize(time_steps, 5); // for benchmarking
    logger->debug("Receiver::receiver_init initialized");
}


void DFDM::receiver::gen_receiver(double delta_t){
    for (uint64_t i = 0; i < ur.rows; ++i) {
        ur.value_set(i, 0, (i+1) * delta_t); 
        sigr.value_set(i, 0, (i+1) * delta_t); 
    }
    logger->debug("Receiver::gen_receiver ur value set");
    for (uint64_t i = 0; i < Nx1; i++) {
        gbx1.value_set(0, i, DFDM::utils::bspline_c(tx1, Nx1, i, kx1, xlocal));
    }
    for (uint64_t i = 0; i < Nx2; i++) {
        gbx2.value_set(0, i, DFDM::utils::bspline_c(tx2, Nx2, i, kx2, xlocal));
    }
    for (uint64_t j = 0; j < Nz1; j++) {
        gbz1.value_set(0, j, DFDM::utils::bspline_c(tz1, Nz1, j, kz1, zlocal));
    }
    for (uint64_t j = 0; j < Nz2; j++) {
        gbz2.value_set(0, j, DFDM::utils::bspline_c(tz2, Nz2, j, kz2, zlocal));
    }

    logger->debug("Receiver::gen_receiver generated");
}


void DFDM::receiver::get_receiver_location(const DFDM::Grid& element_grid){

    xglobal= DFDM::utils::bilinear_interpolate(element_grid.x2d11, Nx1, Nz1, xlocal, zlocal);
    zglobal= DFDM::utils::bilinear_interpolate(element_grid.z2d11, Nx1, Nz1, xlocal, zlocal);

    logger->debug("receiver::get_receiver_location: xlocal: {} , zlocal: {}", xlocal, zlocal);
    logger->debug("receiver::get_receiver_location: Nx1: {} , Nz1: {}", Nx1, Nz1);

    logger->debug("receiver::get_receiver_location: xglobal: {} , zglobal: {}", xglobal, zglobal);
}