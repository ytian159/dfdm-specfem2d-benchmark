#include "source.hpp"
#include "operators.hpp"
#include "grid.hpp"

DFDM::source::source(){
}
void DFDM::source::source_init(uint32_t global_id, uint32_t Nx1_, uint32_t Nz1_, uint32_t order_b1, uint64_t time_steps, std::shared_ptr<spdlog::logger> logger_ ){

    //determining the location of source element.
    // s_element is the global id of the element containing the source.
    // if it is local to this process then initialize the source
    logger = logger_;
    selement_id = global_id;
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

    gbx1.resize(Nx1, 1);
    gbx2.resize(Nx2, 1);
    gbz1.resize(Nz1, 1);
    gbz2.resize(Nz2, 1);

    re.resize(time_steps,1);
    Ft.resize(time_steps,1);
    logger->debug("DFDM::Source::source_init: Source initialized, Nx1: {}, Nz1: {}", Nx1, Nz1);
}

void DFDM::source::gen_source(){
    logger->debug("DFDM::Source::gen_source: Generating source");
    for (uint64_t i = 0; i < Nx1; i++) {
        gbx1.value_set(i, 0, DFDM::utils::bspline_c(tx1, Nx1, i, kx1, xlocal));
    }
    logger->debug("DFDM::Source::gen_source: gbx1 generated");
    for (uint64_t i = 0; i < Nx2; i++) {
        gbx2.value_set(i, 0, DFDM::utils::bspline_c(tx2, Nx2, i, kx2, xlocal));
    }
    logger->debug("DFDM::Source::gen_source: gbx2 generated");
    for (uint64_t j = 0; j < Nz1; j++) {
        gbz1.value_set(j, 0, DFDM::utils::bspline_c(tz1, Nz1, j, kz1, zlocal));
    }
    logger->debug("DFDM::Source::gen_source: gbz1 generated");
    for (uint64_t j = 0; j < Nz2; j++) {
        gbz2.value_set(j, 0, DFDM::utils::bspline_c(tz2, Nz2, j, kz2, zlocal));
    }
    logger->debug("DFDM::Source::gen_source: gbz2 generated");
}

void DFDM::source::src_transform(const DFDM::Operators& ops_x, const DFDM::Operators& ops_z){
    // Source in Grid 1
    logger->debug("DFDM::Source::src_transform: ");
    logger->debug("DFDM::Source::src_transform: gbx1 dimensions: {}x{}", gbx1.rows, gbx1.cols);
    logger->debug("DFDM::Source::src_transform: ops_x.invL11 dimensions: {}x{}", ops_x.invL11.rows, ops_x.invL11.cols);
    logger->debug("DFDM::Source::src_transform: gbz2 dimensions: {}x{}", gbz2.rows, gbz2.cols);
    logger->debug("DFDM::Source::src_transform: ops_z.invL11 dimensions: {}x{}", ops_z.invL11.rows, ops_z.invL11.cols);

    gbx_inv1 = ops_x.invL11.matprod(gbx1);
    gbz_inv2 = ops_z.invL22.matprod(gbz2);
    logger->debug("DFDM::Source::src_transform: gbx_inv1, gbz_inv2 generated");
    gbx_inv2 = ops_x.invL22.matprod(gbx2);
    gbz_inv1 = ops_z.invL11.matprod(gbz1);
    logger->debug("DFDM::Source::src_transform: gbx_inv2, gbz_inv1 generated");
    // inv_sourcegb1 = ops_m.matprod(sourcegb1);
}

void DFDM::source::src_time_gen(uint64_t time_steps, double frequency, double delta_t){
    logger->debug("DFDM::Source::src_time_gen: ");
    for(uint32_t i = 0; i < time_steps; i++){
        re.value_set(i,0, M_PI*frequency*((i+1)*delta_t - 1.5/frequency));
        double resquared = re(i,0)*re(i,0);
        Ft.value_set(i,0,  1e5 * (1.0 - 2.0*resquared)*std::exp(-resquared)); 
    }
    // Ft.print_file("Ft.out");   
    logger->debug("DFDM::Source::src_time_gen: re, Ft generated");
}

void DFDM::source::invMsg_gen(){
    logger->debug("DFDM::Source::invMsg_gen: ");
    invMsg12.resize(Nx1, Nz2);
    invMsg21.resize(Nx2, Nz1);
    for (uint32_t i = 0; i < Nx1; i++) {
        for (uint32_t j = 0; j < Nz2; j++) {
            invMsg12.value_set(i, j, gbx_inv1(i, 0) * gbz_inv2(j, 0));
        }
    }

    logger->debug("DFDM::Source::invMsg_gen: invMsg12 generated");
    for (uint32_t i = 0; i < Nx2; i++) {
        for (uint32_t j = 0; j < Nz1; j++) {
            invMsg21.value_set(i, j, gbx_inv2(i, 0) * gbz_inv1(j, 0));
        }
    }
    // invMsg12.print_file("invMsg12");
    // invMsg21.print_file("invMsg21");
    logger->debug("DFDM::Source::invMsg_gen: invMsg21 generated");
}

void DFDM::source::get_source_location(const DFDM::Grid& element_grid){

    xglobal= DFDM::utils::bilinear_interpolate(element_grid.x2d11, Nx1, Nz1, xlocal, zlocal);
    zglobal= DFDM::utils::bilinear_interpolate(element_grid.z2d11, Nx1, Nz1, xlocal, zlocal);

    logger->info("DFDM::Source::get_source_location: xlocal: {} , zlocal: {}", xlocal, zlocal);
    logger->debug("DFDM::Source::get_source_location: Nx1: {} , Nz1: {}", Nx1, Nz1);

    logger->info("DFDM::Source::get_source_location: xglobal: {} , zglobal: {}", xglobal, zglobal);
}