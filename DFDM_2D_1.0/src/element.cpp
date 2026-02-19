#include "element.hpp"

DFDM::Element2D::Element2D(DFDM::simulation& sim_in, DFDM::ElementData& el_data, std::shared_ptr<spdlog::logger> logger_) :ops_x(DFDM::OP_TYPE::X), ops_z(DFDM::OP_TYPE::Z){
    logger = logger_;
    element_id = el_data.global_id;
    global_id = el_data.global_id;
    region_id = el_data.region_id;
    domain_id = el_data.domain_id;
   
    nax = el_data.nax;
    naz = el_data.naz;
    Nx1 = el_data.NAX;
    Nz1 = el_data.NAZ;
    Nx2 = el_data.NAX - 1;
    Nz2 = el_data.NAZ - 1;


    owner_cpu = el_data.owner_cpu;

    order_b1 = sim_in.order_b1;
    // order_b2 = sim_in.order_b2;

    gauss_order = sim_in.gauss_order;
    

    for(uint32_t i = 0; i < 4; i++){
        neighbors[i] = el_data.neighbors[i];
        neighbor_process[i] = el_data.neighbor_process[i];
        face_connections_nbrs[i] = el_data.face_connections_nbrs[i];
        face_connections_self[i] = el_data.face_connections_self[i];
        rotation_nbr[i] = el_data.rotation_nbr[i];
        rotation_self[i] = el_data.rotation_self[i];
        neighbor_Nx1[i] = el_data.neighbor_Nx1[i];
        neighbor_Nz1[i] = el_data.neighbor_Nz1[i];
        neighbor_px1[i] = order_b1;//el_data.neighbor_px1[i]; this will later read from mesh
        neighbor_pz1[i] = order_b1;//el_data.neighbor_pz1[i];
    }

    // setting alpha values
    for (uint32_t i = 0; i < 4; i++) {
        if (neighbors[i] != -1) {
            switch (i) {
                case 0:
                    alpha_mo = 0.5;
                    break;
                case 1:
                    alpha_po = 0.5;
                    break;
                case 2:
                    alpha_om = 0.5;
                    break;
                case 3:
                    alpha_op = 0.5;
                    break;
            }
        } else {
            switch (i) {
                case 0:
                    alpha_mo = 0;
                    break;
                case 1:
                    alpha_po = 0;
                    break;
                case 2:
                    alpha_om = 0;
                    break;
                case 3:
                    alpha_op = 0;
                    break;
            }
        }
    }


    // each send requires 2 requests and 2 statuses
    // each receive requires 2 requests and 2 statuses
    // each elements sends and receives at most 4 matrix boundaries.
    // each boundary requires 2 sends and 2 receives
    // so 2*4*2 = 16
    // 0-3 for left, 4-7 for right, 8-11  for bottom, 12-15 for top
    for(uint32_t i = 0; i < 16; i++){
        MPI_Status temp_status;
        MPI_Request temp_request;
        s_status.push_back(temp_status);
        s_request.push_back(temp_request);
    }

    for(uint32_t i = 0; i < 16; i++){
        MPI_Status temp_status;
        MPI_Request temp_request;
        r_status.push_back(temp_status);
        r_request.push_back(temp_request);
    }
}


void DFDM::Element2D::assign_locality(std::vector<uint32_t>& local_elements){
    // if the neighbor is in the local elements, set the locality to true
    for (uint32_t i = 0; i < 4; i++) {
        if (std::find(local_elements.begin(), local_elements.end(), neighbors[i]) != local_elements.end()) {
            locality[i] = 1;
        }
    }
}

DFDM::Element2D& DFDM::Element2D::getElement(uint32_t target_el, std::vector<DFDM::Element2D>& loc_elements){
    for(uint32_t t = 0; t < loc_elements.size(); t++){
        if(loc_elements[t].global_id == target_el){
            logger->debug("DFDM::Element2D::getElement:: Requested local element found: {}", loc_elements[t].global_id);
            return loc_elements[t];
        }
    }
    std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
    ": Requested local element not found" << std::endl;
    exit(-1);
}

void DFDM::Element2D::initialize(const DFDM::simulation& sim, const DFDM::ElementData& el_in){
    logger->debug("initialize:: Setting up element grid in element: {}", global_id);
    // initialize element grid
    element_grid.grid_init(el_in, sim, logger);
    // order_b1 and order_b2 right now are both same, px1, pz1 in matlab
    logger->debug("initialize::Setting up x ops in element: {}", global_id);
    ops_x.init_operators(Nx1, order_b1, sim.gauss_order, global_id, logger); // initialize matrix operators for x dim
    logger->debug("initialize:: Setting up z ops in element: {}", global_id);
    ops_z.init_operators(Nz1, order_b1, sim.gauss_order, global_id, logger); // initialize matrix operators for z dim
    logger->debug("initialize:: Setting up State in element: {}", global_id);
    state.init(Nx1, Nz1, sim.time_steps, logger); // initialize State matrices
    logger->debug("initialize:: Setting element: {} completed", global_id);

    //boundaries are now inside state.
}

void DFDM::Element2D::gen_dfdm_matrices(){
    //first gen knots
    logger->debug("gen_dfdm_matrices:: gen_knots from xops in element: {}", global_id);
    ops_x.gen_knots();
    logger->debug("gen_dfdm_matrices:: gen_knots from zops in element: {}", global_id);
    ops_z.gen_knots();
    logger->debug("gen_dfdm_matrices:: gen_knots completed for element: {}", global_id);
    // compute gaussian quadrature params
    logger->debug("gen_dfdm_matrices:: gen_gauss_params from xops in element: {}", global_id);
    ops_x.gen_gauss_params();
    logger->debug("gen_dfdm_matrices:: gen_gauss_params from zops in element: {}", global_id);
    ops_z.gen_gauss_params();
    logger->debug("gen_dfdm_matrices:: gen_gauss_params completed for element: {}", global_id);
    // compute bsplines
    logger->debug("gen_dfdm_matrices:: gen_bspline from xops in element: {}", global_id);
    ops_x.gen_bspline();
    
    logger->debug("gen_dfdm_matrices:: gen_bspline from zops in element: {}", global_id);
    ops_z.gen_bspline();
    logger->debug("gen_dfdm_matrices:: gen_bspline completed for element: {}", global_id);

    // compute matrix operators using gaussian quadrature.
    logger->debug("gen_dfdm_matrices:: compute_operators from xops in element: {}", global_id);
    ops_x.compute_operators();
    logger->debug("gen_dfdm_matrices:: compute_operators from zops in element: {}", global_id);
    ops_z.compute_operators();
    logger->debug("gen_dfdm_matrices:: compute_operators completed for element: {}", global_id);

    // compute connection matrices for neighbors (op, x dim)
    logger->debug("gen_dfdm_matrices:: gen_connected_ops for element: {}", global_id);
    {   //compute connection matrices for left neighbor (op, z dim)
        int nbr_id = 0;
        std::vector<DFDM::matrix<double>> conn_matrices;
        if(neighbors[nbr_id] != -1){
            conn_matrices = ops_z.gen_connected_ops(neighbors[nbr_id], face_connections_nbrs[nbr_id], neighbor_Nx1[nbr_id], neighbor_Nz1[nbr_id], neighbor_px1[nbr_id], neighbor_pz1[nbr_id]);
            ops_z.TA21 = conn_matrices[0];
            ops_z.TA12 = conn_matrices[1];
            ops_z.TA11 = conn_matrices[2];
            ops_z.TA22 = conn_matrices[3];   
            logger->debug("gen_dfdm_matrices::gen_connected_ops completed for nbr_id:{} for element:{}. Dimensions of ops_z.TA21: {}x{}, ops_z.TA12: {}x{}, ops_z.TA11: {}x{}, ops_z.TA22: {}x{}", 
            nbr_id, global_id, ops_z.TA21.rows, ops_z.TA21.cols, ops_z.TA12.rows, ops_z.TA12.cols, ops_z.TA11.rows, ops_z.TA11.cols, ops_z.TA22.rows, ops_z.TA22.cols);
        }else{
            logger->debug("gen_dfdm_matrices::gen_connected_ops skipped for nbr_id:{} for element:{} because nbr doesn't exist", nbr_id, global_id);
        }
   
        //compute connection matrices for right neighbor (op, z dim)
        nbr_id = 1;
        if(neighbors[nbr_id] != -1){        
            conn_matrices = ops_z.gen_connected_ops(neighbors[nbr_id], face_connections_nbrs[nbr_id], neighbor_Nx1[nbr_id], neighbor_Nz1[nbr_id], neighbor_px1[nbr_id], neighbor_pz1[nbr_id]);
            ops_z.TB21 = conn_matrices[0];
            ops_z.TB12 = conn_matrices[1];
            ops_z.TB11 = conn_matrices[2];
            ops_z.TB22 = conn_matrices[3];
            logger->debug("gen_dfdm_matrices::gen_connected_ops completed for nbr_id:{} for element:{}. Dimensions of ops_z.TB21: {}x{}, ops_z.TB12: {}x{}, ops_z.TB11: {}x{}, ops_z.TB22: {}x{}", 
            nbr_id, global_id, ops_z.TB21.rows, ops_z.TB21.cols, ops_z.TB12.rows, ops_z.TB12.cols, ops_z.TB11.rows, ops_z.TB11.cols, ops_z.TB22.rows, ops_z.TB22.cols);
        }else{
            logger->debug("gen_dfdm_matrices::gen_connected_ops skipped for nbr_id:{} for element:{} because nbr doesn't exist", nbr_id, global_id);
        }
   
        //compute connection matrices for lower neighbor (op, x dim)
        nbr_id = 2;
        if(neighbors[nbr_id] != -1){        
            conn_matrices = ops_x.gen_connected_ops(neighbors[nbr_id], face_connections_nbrs[nbr_id], neighbor_Nx1[nbr_id], neighbor_Nz1[nbr_id], neighbor_px1[nbr_id], neighbor_pz1[nbr_id]);
            ops_x.TA21 = conn_matrices[0];
            ops_x.TA12 = conn_matrices[1];
            ops_x.TA11 = conn_matrices[2];
            ops_x.TA22 = conn_matrices[3];      
            logger->debug("gen_dfdm_matrices::gen_connected_ops completed for nbr_id:{} for element:{}. Dimensions of ops_x.TA21: {}x{}, ops_x.TA12: {}x{}, ops_x.TA11: {}x{}, ops_x.TA22: {}x{}", 
            nbr_id, global_id, ops_x.TA21.rows, ops_x.TA21.cols, ops_x.TA12.rows, ops_x.TA12.cols, ops_x.TA11.rows, ops_x.TA11.cols, ops_x.TA22.rows, ops_x.TA22.cols);
         }else{
            logger->debug("gen_dfdm_matrices::gen_connected_ops skipped for nbr_id:{} for element:{} because nbr doesn't exist", nbr_id, global_id);
        }   
//compute connection matrices for upper neighbor (op, x dim)
        nbr_id = 3;
        if(neighbors[nbr_id] != -1){        
            conn_matrices = ops_x.gen_connected_ops(neighbors[nbr_id], face_connections_nbrs[nbr_id], neighbor_Nx1[nbr_id], neighbor_Nz1[nbr_id], neighbor_px1[nbr_id], neighbor_pz1[nbr_id]);
            ops_x.TB21 = conn_matrices[0];
            ops_x.TB12 = conn_matrices[1];
            ops_x.TB11 = conn_matrices[2];
            ops_x.TB22 = conn_matrices[3];
            logger->debug("gen_dfdm_matrices::gen_connected_ops completed for nbr_id:{} for element:{}. Dimensions of  ops_x.TB21: {}x{}, ops_x.TB12: {}x{}, ops_x.TB11: {}x{}, ops_x.TB22: {}x{}", 
            nbr_id, global_id,  ops_x.TB21.rows, ops_x.TB21.cols, ops_x.TB12.rows, ops_x.TB12.cols, ops_x.TB11.rows, ops_x.TB11.cols, ops_x.TB22.rows, ops_x.TB22.cols);
        }else{
            logger->debug("gen_dfdm_matrices::gen_connected_ops skipped for nbr_id:{} for element:{} because nbr doesn't exist", nbr_id, global_id);
        }   

    }
    // if(global_id == 0){
    //     ops_z.TA21.print_file("./output/ops_z.TA21");
    //     ops_z.TA12.print_file("./output/ops_z.TA12");
    //     ops_z.TA11.print_file("./output/ops_z.TA11");
    //     ops_z.TA22.print_file("./output/ops_z.TA22");

    //     ops_z.TB21.print_file("./output/ops_z.TB21");
    //     ops_z.TB12.print_file("./output/ops_z.TB12");
    //     ops_z.TB11.print_file("./output/ops_z.TB11");
    //     ops_z.TB22.print_file("./output/ops_z.TB22");

    //     ops_x.TA21.print_file("./output/ops_x.TA21");
    //     ops_x.TA12.print_file("./output/ops_x.TA12");
    //     ops_x.TA11.print_file("./output/ops_x.TA11");
    //     ops_x.TA22.print_file("./output/ops_x.TA22");

    //     ops_x.TB21.print_file("./output/ops_x.TB21");
    //     ops_x.TB12.print_file("./output/ops_x.TB12");
    //     ops_x.TB11.print_file("./output/ops_x.TB11");
    //     ops_x.TB22.print_file("./output/ops_x.TB22");
    // }
    logger->debug("gen_dfdm_matrices:: DFDM matrices for element: {} completed", global_id);
}


void DFDM::Element2D::update_wavefield(uint64_t time_step, const DFDM::simulation& sim){
    // if(global_id == 1){
    //     state.U21.print_file("U21_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    //     state.U12.print_file("U12_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    //     state.dU2dtt12.print_file("dU2dtt12_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    // }
    auto deltasq = (sim.delta_t * sim.delta_t);
    state.U12 = (state.U12_1*2 - state.U12_0) + state.dU2dtt12*deltasq;
    state.U12_0 = state.U12_1;
    state.U12_1 = state.U12;

    state.U21 = (state.U21_1*2 - state.U21_0) + state.dU2dtt21*deltasq;
    state.U21_0 = state.U21_1;
    state.U21_1 = state.U21;

    // if(global_id == 0 || global_id == 1 || global_id == 4){
        // state.U12.print_file("./output/state.U12_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.dU2dtt12.print_file("./output/state.dU2dtt12_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U21.print_file("./output/state.U21_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.dU2dtt21.print_file("./output/state.dU2dtt21_"+std::to_string(global_id)+"_"+std::to_string(time_step));
    // }

    // if(global_id == 1){
    //     state.U21.print_file("after_U21_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    //     state.U12.print_file("after_U12_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    //     state.dU2dtt12.print_file("after_dU2dtt12_"+std::to_string(time_step)+"_"+std::to_string(global_id));
    // }
}


void DFDM::Element2D::compute_Uboundary_local(uint64_t time_step){
    // compute U boundary for x direction
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: computing x direction");
    auto tmpU12 = ops_x.invL11_t.matprod(state.U12); // convert to B1 basis
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: coverted to B1 basis, tmpU12 size: {}x{}, U12mo_inn size: {}x{}, element:{}", tmpU12.rows, tmpU12.cols, state.U12mo_inn.rows, state.U12mo_inn.cols, global_id);
    DFDM::matrix<double>::row_copy(tmpU12, 0, state.U12mo_inn, 0); // copy left boundary
    DFDM::matrix<double>::row_copy(tmpU12, tmpU12.rows-1, state.U12po_inn, 0); // copy right boundary
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: copied left and right boundaries");

    tmpU12 = state.U12.transpose_inplace();
    tmpU12 = ops_z.invL22_t.matprod(tmpU12); // convert to B2 basis
    tmpU12 = tmpU12.transpose_inplace(); 
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: coverted to B2 basis");
    DFDM::matrix<double>::col_copy(tmpU12, 0, state.U12om_inn, 0); // copy bottom boundary
    DFDM::matrix<double>::col_copy(tmpU12, tmpU12.cols - 1, state.U12op_inn, 0); // copy top boundary
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: copied bottom and top boundaries");

    logger->debug("DFDM::Element2D::compute_Uboundary_local:: computing z direction");
    // compute U boundary for z direction
    auto tmpU21 = ops_x.invL22_t.matprod(state.U21); // convert to B1 basis
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: coverted to B1 basis");
    DFDM::matrix<double>::row_copy(tmpU21, 0, state.U21mo_inn, 0); // copy left boundary
    DFDM::matrix<double>::row_copy(tmpU21, tmpU21.rows-1, state.U21po_inn, 0); // copy right boundary
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: copied left and right boundaries");
    tmpU21 = state.U21.transpose_inplace();
    tmpU21 = ops_z.invL11_t.matprod(tmpU21); // convert to B2 basis
    tmpU21 = tmpU21.transpose_inplace(); 
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: coverted to B2 basis");
    DFDM::matrix<double>::col_copy(tmpU21, 0, state.U21om_inn, 0); // copy bottom boundary
    DFDM::matrix<double>::col_copy(tmpU21, tmpU21.cols - 1, state.U21op_inn, 0); // copy top boundary
    logger->debug("DFDM::Element2D::compute_Uboundary_local:: copied bottom and top boundaries");
    // if(global_id == 0 || global_id == 1 || global_id == 4){
        // state.U12mo_inn.print_file("./output/state.U12mo_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U12po_inn.print_file("./output/state.U12po_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U12om_inn.print_file("./output/state.U12om_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U12op_inn.print_file("./output/state.U12op_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));

        // state.U21mo_inn.print_file("./output/state.U21mo_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U21po_inn.print_file("./output/state.U21po_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U21om_inn.print_file("./output/state.U21om_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));
        // state.U21op_inn.print_file("./output/state.U21op_inn_"+std::to_string(global_id)+"_"+std::to_string(time_step));  
    // }
}



void DFDM::Element2D::send_Uboundary(std::vector<DFDM::Element2D>& loc_elements, uint64_t time_step){
// make sure to consider the difference between local and remote elements.
// only send if its remote
// loop over al the neighbors, prepare boundary for sending, if it is local 
// then copy it locally in the neighbor element directly

    std::vector<DFDM::matrix<double>> to_send_U12_vec(4);
    std::vector<DFDM::matrix<double>> to_send_U21_vec(4);
    logger->debug("DFDM::Element2D::send_Uboundary:: sending U boundary");
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            // auto to_send_U12 = get_face(state.U12mo_inn, state.U12po_inn, state.U12om_inn, state.U12op_inn, face_connections_self[nbr], rotation_self[nbr]);
            auto to_send_U12 = get_face(state.U12mo_inn, state.U12po_inn, state.U12om_inn, state.U12op_inn, nbr, rotation_self[nbr]);
            auto to_send_U21 = get_face(state.U21mo_inn, state.U21po_inn, state.U21om_inn, state.U21op_inn, nbr, rotation_self[nbr]);
            logger->debug("DFDM::Element2D::send_Uboundary:: extracted faces to send for U12 and U21 for neighbor:{}", neighbors[nbr]);
            if(nbr == 0 || nbr == 1){
                to_send_U12 = to_send_U12.transpose_inplace();
                to_send_U12 = ops_z.invL22_t.matprod(to_send_U12);
                to_send_U21 = to_send_U21.transpose_inplace();
                to_send_U21 = ops_z.invL11_t.matprod(to_send_U21);
            }else if(nbr == 2 || nbr == 3){
                to_send_U12 = ops_x.invL11_t.matprod(to_send_U12);        
                to_send_U21 = ops_x.invL22_t.matprod(to_send_U21);
            }else{
                std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                ": Invalid face connection when preparing U boundary for sending " << std::endl;
                exit(-1);
            }
            if(locality[nbr] == 1){ // if its local to this rank
                // get the element
                logger->debug("DFDM::Element2D::send_Uboundary:: neighbor:{} for current element:{} is local", neighbors[nbr], global_id);
                decltype(auto) target_el = get_local_element(neighbors[nbr], loc_elements);
                switch(face_connections_nbrs[nbr]){// which of nbr's face am I facing?
                    case 0:
                        target_el.state.U12mo_out = to_send_U12;
                        target_el.state.U21mo_out = to_send_U21;
                        logger->debug("DFDM::Element2D::send_Uboundary:: copied U12 ({}x{}) and U21 ({}x{}) to local neighbor:{}'s face:{}, self_id:{}, my_face:{}, target_id:{}", target_el.state.U12mo_out.rows, target_el.state.U12mo_out.cols, to_send_U21.rows, to_send_U21.cols, neighbors[nbr], face_connections_nbrs[nbr], global_id, nbr, target_el.global_id);
                        break;
                    case 1:
                        target_el.state.U12po_out = to_send_U12;
                        target_el.state.U21po_out = to_send_U21;
                        logger->debug("DFDM::Element2D::send_Uboundary:: copied U12 ({}x{}) and U21 ({}x{}) to local neighbor:{}'s face:{}, self_id:{}, my_face:{}", to_send_U12.rows, to_send_U12.cols, to_send_U21.rows, to_send_U21.cols, neighbors[nbr], face_connections_nbrs[nbr], global_id, nbr);
                        break;
                    case 2:
                        logger->debug("DFDM::Element2D::send_Uboundary:: copying U12 ({}x{}) and U21 ({}x{}) to local neighbor:{}'s face:{}, self_id:{}, my_face:{}", to_send_U12.rows, to_send_U12.cols, to_send_U21.rows, to_send_U21.cols, neighbors[nbr], face_connections_nbrs[nbr], global_id, nbr);
                        target_el.state.U12om_out = to_send_U12;
                        target_el.state.U21om_out = to_send_U21;
                        logger->debug("DFDM::Element2D::send_Uboundary:: copied U12 ({}x{}) and U21 ({}x{}) to local neighbor:{}'s face:{}, self_id:{}, my_face:{}", to_send_U12.rows, to_send_U12.cols, to_send_U21.rows, to_send_U21.cols, neighbors[nbr], face_connections_nbrs[nbr], global_id, nbr);
                        break;
                    case 3:
                        target_el.state.U12op_out = to_send_U12;
                        target_el.state.U21op_out = to_send_U21;
                        // if(global_id == 1){
                        //     if(neighbors[nbr] == 0)target_el.state.U12op_out.print_file("./output/copied_target_el.state.U12op_out_"+std::to_string(time_step));
                        // }
                        logger->debug("DFDM::Element2D::send_Uboundary:: copied U12 ({}x{}) and U21 ({}x{}) to local neighbor:{}'s face:{}, self_id:{}, my_face:{}", to_send_U12.rows, to_send_U12.cols, to_send_U21.rows, to_send_U21.cols, neighbors[nbr], face_connections_nbrs[nbr], global_id, nbr);
                        break;
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when preparing U boundary for sending " << std::endl;
                        exit(-1);
                }
            }else{
                // send to remote
                // when sending generate a unique tag depending on which face is being sent
                // use the target element's face that is facing you to generate the tag
                logger->debug("DFDM::Element2D::send_Uboundary:: neighbor:{} for current element:{} is remote on proces:{}", neighbors[nbr], global_id, neighbor_process[nbr]);
                switch(face_connections_nbrs[nbr]){ // need the switch statement to generate unique tags for boundary types
                    case 0:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::MO, neighbor_process[nbr], neighbors[nbr]);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::MO, neighbor_process[nbr], neighbors[nbr]);

                        to_send_U12_vec[nbr] = to_send_U12;
                        to_send_U21_vec[nbr] = to_send_U21;

                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU12, to_send_U12_vec[nbr], s_request[nbr]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U12 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, from elem:{} to nbr's face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U12.rows, to_send_U12.cols, tagU12, global_id, face_connections_nbrs[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU21, to_send_U21_vec[nbr], s_request[nbr+4]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U21 to process:{} for targ_el:{} with dims:{}x{} and tag:{} from elem:{} to nbr's face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U21.rows, to_send_U21.cols, tagU21, global_id, face_connections_nbrs[nbr]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        // mpi_wait(s_request[0], s_request[2]);
                        // mpi_wait(s_request[2], s_request[3]);
                        break;
                    }
                    case 1:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::PO, neighbor_process[nbr], neighbors[nbr]);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::PO, neighbor_process[nbr], neighbors[nbr]);
                        
                        to_send_U12_vec[nbr] = to_send_U12;
                        to_send_U21_vec[nbr] = to_send_U21;

                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU12, to_send_U12_vec[nbr], s_request[nbr]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U12 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U12.rows, to_send_U12.cols, tagU12, face_connections_nbrs[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU21, to_send_U21_vec[nbr], s_request[nbr+4]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U21 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U21.rows, to_send_U21.cols, tagU21, face_connections_nbrs[nbr]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        // mpi_wait(s_request[4], s_request[5]);
                        // mpi_wait(s_request[6], s_request[7]);
                        break;
                    }
                    case 2:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::OM, neighbor_process[nbr], neighbors[nbr]);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::OM, neighbor_process[nbr], neighbors[nbr]);
                        
                        to_send_U12_vec[nbr] = to_send_U12;
                        to_send_U21_vec[nbr] = to_send_U21;

                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU12, to_send_U12_vec[nbr], s_request[nbr]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U12 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U12.rows, to_send_U12.cols, tagU12, face_connections_nbrs[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU21, to_send_U21_vec[nbr], s_request[nbr+4]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U21 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U21.rows, to_send_U21.cols, tagU21, face_connections_nbrs[nbr]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        // DFDM::send_mat_nb(neighbor_process[nbr], tagU12, to_send_U12, s_request[nbr], s_request[nbr+1]);
                        // DFDM::send_mat_nb(neighbor_process[nbr], tagU21, to_send_U21, s_request[nbr+2], s_request[nbr+3]);
                        // mpi_wait(s_request[8], s_request[9]);
                        // mpi_wait(s_request[10], s_request[11]);
                        break;
                    }
                    case 3:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::OP, neighbor_process[nbr], neighbors[nbr]);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::OP, neighbor_process[nbr], neighbors[nbr]);
                        
                        to_send_U12_vec[nbr] = to_send_U12;
                        to_send_U21_vec[nbr] = to_send_U21;

                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU12, to_send_U12_vec[nbr], s_request[nbr]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U12 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U12.rows, to_send_U12.cols, tagU12, face_connections_nbrs[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagU21, to_send_U21_vec[nbr], s_request[nbr+4]);
                        logger->debug("DFDM::Element2D::send_Uboundary:: sent U21 to process:{} for targ_el:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], neighbors[nbr], to_send_U21.rows, to_send_U21.cols, tagU21, face_connections_nbrs[nbr]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        // DFDM::send_mat_nb(neighbor_process[nbr], tagU12, to_send_U12, s_request[nbr], s_request[nbr+1]);
                        // DFDM::send_mat_nb(neighbor_process[nbr], tagU21, to_send_U21, s_request[nbr+2], s_request[nbr+3]);
                        // mpi_wait(s_request[12], s_request[13]);
                        // mpi_wait(s_request[14], s_request[15]);
                        break;
                    }
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when preparing U boundary for sending " << std::endl;
                        exit(-1);
                }
            }

        }
        logger->debug("DFDM::Element2D::send_Uboundary:: completed sending U boundary for nbr:{}", nbr);
    }   
}

void DFDM::Element2D::receive_Uboundary(std::vector<DFDM::Element2D>& loc_elements, uint64_t time_step){

    std::vector<DFDM::matrix<double>> to_recv_U12;
    std::vector<DFDM::matrix<double>> to_recv_U21;
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        to_recv_U12.push_back(DFDM::matrix<double>());
        to_recv_U21.push_back(DFDM::matrix<double>());
    }

    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            if(locality[nbr] == 0){ // if its remote to this rank (1 means local, 0 means remote)
                //then I need to receive the boundary from remote
                // when receiving use your own process id, element id and your face that is facing
                // the neighbor to generate the tag (self)
                // depending on the face connection, assign the received boundary to the correct variable
                switch(face_connections_self[nbr]){
                    case 0:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::MO, owner_cpu, global_id); //
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::MO, owner_cpu, global_id);
                        
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_U12[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_U21[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_U12[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_U21[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                        }

                        logger->debug("DFDM::Element2D::receive_Uboundary:: receiving U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU12, to_recv_U12[nbr], r_request[nbr]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U21 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], global_id, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols, tagU21, face_connections_self[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU21, to_recv_U21[nbr], r_request[nbr+4]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U21 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{}", neighbor_process[nbr], global_id, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols, tagU21, face_connections_self[nbr]);
                        break;
                    }
                    case 1:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::PO, owner_cpu, global_id);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::PO, owner_cpu, global_id);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_U12[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_U21[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_U12[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_U21[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                        }
                        logger->debug("DFDM::Element2D::receive_Uboundary:: receiving U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU12, to_recv_U12[nbr], r_request[nbr]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU21, to_recv_U21[nbr], r_request[nbr+4]);
                        break;
                    }
                    case 2:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::OM, owner_cpu, global_id);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::OM, owner_cpu, global_id);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_U12[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_U21[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_U12[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_U21[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                        }
                        logger->debug("DFDM::Element2D::receive_Uboundary:: receiving U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU12, to_recv_U12[nbr], r_request[nbr]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU21, to_recv_U21[nbr], r_request[nbr+4]);
                        
                        break;
                    }
                    case 3:{
                        auto tagU12 = DFDM::generateTag(DFDM::BdryType::U12, DFDM::BdryDirection::OP, owner_cpu, global_id);
                        auto tagU21 = DFDM::generateTag(DFDM::BdryType::U21, DFDM::BdryDirection::OP, owner_cpu, global_id);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_U12[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_U21[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_U12[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_U21[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                        }
                        logger->debug("DFDM::Element2D::receive_Uboundary:: receiving U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU12, to_recv_U12[nbr], r_request[nbr]);
                        
                        logger->debug("DFDM::Element2D::receive_Uboundary:: received U12 from process:{} for self_element_id:{} with dims:{}x{} and tag:{}, face:{} from elem:{}", neighbor_process[nbr], global_id, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, tagU12, face_connections_self[nbr], neighbors[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagU21, to_recv_U21[nbr], r_request[nbr+4]);                       
                        
                        break;
                    }
                }
            }else{
                //if neighbor is local then we already have the boundary in *out variables so do nothing
                logger->debug("DFDM::Element2D::receive_Uboundary:: neighbor:{} for current element:{} is local on proces:{}, so no remote receive requested", neighbors[nbr], global_id, neighbor_process[nbr]);
            }  
        }
    }
    logger->debug("DFDM::Element2D::receive_Uboundary::completed the non blocking receive loop");
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            if(locality[nbr] == 0){
                switch(nbr){
                    case 0:{
                        logger->debug("DFDM::Element2D::receive_Uboundary:: waiting for U12 and U21 for self_element_id:{}, nbr:{}", global_id, nbr);
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: completed waiting for U12 and U21 for self_element_id:{}, nbr:{}, to_recv_U12.rows:{}, to_recv_U12.cols:{} and to_recv_U21.rows:{}, to_recv_U21.cols:{}", global_id, nbr, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols);

                        state.U12mo_out = to_recv_U12[nbr];
                        state.U21mo_out = to_recv_U21[nbr];
                        break;
                    }
                    case 1:{
                        logger->debug("DFDM::Element2D::receive_Uboundary:: waiting for U12 and U21 for self_element_id:{}, nbr:{}", global_id, nbr);
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: completed waiting for U12 and U21 for self_element_id:{}, nbr:{} to_recv_U12.rows:{}, to_recv_U12.cols:{} and to_recv_U21.rows:{}, to_recv_U21.cols:{}", global_id, nbr, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols);
                        state.U12po_out = to_recv_U12[nbr];
                        state.U21po_out = to_recv_U21[nbr];
                        break;
                    }
                    case 2:{
                        logger->debug("DFDM::Element2D::receive_Uboundary:: waiting for U12 and U21 for self_element_id:{}, nbr:{}", global_id, nbr);
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: completed waiting for U12 and U21 for self_element_id:{}, nbr:{} to_recv_U12.rows:{}, to_recv_U12.cols:{} and to_recv_U21.rows:{}, to_recv_U21.cols:{}", global_id, nbr, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols);
                        state.U12om_out = to_recv_U12[nbr];
                        state.U21om_out = to_recv_U21[nbr];
                        break;
                    }
                    case 3:{
                        logger->debug("DFDM::Element2D::receive_Uboundary:: waiting for U12 and U21 for self_element_id:{}, nbr:{}", global_id, nbr);
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: completed waiting for U12 and U21 for self_element_id:{}, nbr:{} to_recv_U12.rows:{}, to_recv_U12.cols:{} and to_recv_U21.rows:{}, to_recv_U21.cols:{}", global_id, nbr, to_recv_U12[nbr].rows, to_recv_U12[nbr].cols, to_recv_U21[nbr].rows, to_recv_U21[nbr].cols);
                        state.U12op_out = to_recv_U12[nbr];
                        state.U21op_out = to_recv_U21[nbr];
                        break;
                    }
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when synchronizing U boundaries in receiver " << std::endl;
                        exit(-1);
                }
            }
        }
    }
    logger->debug("DFDM::Element2D::receive_Uboundary:: completed the blocking receive loop");

    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            // applying transformation to the boundary and converting basis
            switch(nbr){ // here nbr can also be replaced with face_connections_self[nbr] as they are the same
                case 0:
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        // using ops_z.TA* because we are getting from left neighbor (z dim)
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12mo_out for self_element_id:{}. Dimensions of ops_z.TA22: {}x{}, state.U12mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.TA22.rows, ops_z.TA22.cols, state.U12mo_out.rows, state.U12mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U12mo_out = ops_z.TA22.matprod(state.U12mo_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12mo_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.U12mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.U12mo_out.rows, state.U12mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U12mo_out = ops_z.invL22.matprod(state.U12mo_out);

                        // state.U21mo_out = to_recv_U21;
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21mo_out for self_element_id:{}. Dimensions of ops_z.TA11: {}x{}, state.U21mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.TA11.rows, ops_z.TA11.cols, state.U21mo_out.rows, state.U21mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U21mo_out = ops_z.TA11.matprod(state.U21mo_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21mo_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.U21mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.U21mo_out.rows, state.U21mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U21mo_out = ops_z.invL11.matprod(state.U21mo_out);
                        // to match with the in boundaries
                        state.U12mo_out.reshape(1,state.U12mo_out.data_vec.size());
                        state.U21mo_out.reshape(1,state.U21mo_out.data_vec.size());
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        // using ops_z.TA* because we are getting from left neighbor (z dim)
                        // state.U12mo_out = to_recv_U12;
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12mo_out for self_element_id:{}. Dimensions of ops_z.TA21: {}x{}, state.U12mo_out: {}x{}, face_conn:{}, my_face:{}, nbr_id:{}", global_id, ops_z.TA21.rows, ops_z.TA21.cols, state.U12mo_out.rows, state.U12mo_out.cols, face_connections_nbrs[nbr], nbr, neighbors[nbr]);
                        state.U12mo_out = ops_z.TA21.matprod(state.U12mo_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12mo_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.U12mo_out: {}x{}, face_conn:{}, my_face:{}, nbr_id:{}", global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.U12mo_out.rows, state.U12mo_out.cols, face_connections_nbrs[nbr], nbr, neighbors[nbr]);
                        state.U12mo_out = ops_z.invL22.matprod(state.U12mo_out);

                        // state.U21mo_out = to_recv_U21;
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21mo_out for self_element_id:{}. Dimensions of ops_z.TA12: {}x{}, state.U21mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.TA12.rows, ops_z.TA12.cols, state.U21mo_out.rows, state.U21mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U21mo_out = ops_z.TA12.matprod(state.U21mo_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21mo_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.U21mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.U21mo_out.rows, state.U21mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.U21mo_out = ops_z.invL11.matprod(state.U21mo_out);

                        // to match with the in boundaries
                        state.U12mo_out.reshape(1,state.U12mo_out.data_vec.size());
                        state.U21mo_out.reshape(1,state.U21mo_out.data_vec.size());
                    }
                    break;
                case 1:
                        // using ops_z.TA* because we are getting from right neighbor (z dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12po_out for self_element_id:{}. Dimensions of ops_z.TB22: {}x{}, state.U12po_out: {}x{}, face_conn:{}", global_id, ops_z.TB22.rows, ops_z.TB22.cols, state.U12po_out.rows, state.U12po_out.cols, face_connections_nbrs[nbr]);
                        state.U12po_out = ops_z.TB22.matprod(state.U12po_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12po_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.U12po_out: {}x{}, face_conn:{}", global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.U12po_out.rows, state.U12po_out.cols, face_connections_nbrs[nbr]);
                        state.U12po_out = ops_z.invL22.matprod(state.U12po_out);

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21po_out for self_element_id:{}. Dimensions of ops_z.TB11: {}x{}, state.U21po_out: {}x{}, face_conn:{}", global_id, ops_z.TB11.rows, ops_z.TB11.cols, state.U21po_out.rows, state.U21po_out.cols, face_connections_nbrs[nbr]);
                        state.U21po_out = ops_z.TB11.matprod(state.U21po_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21po_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.U21po_out: {}x{}, face_conn:{}", global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.U21po_out.rows, state.U21po_out.cols, face_connections_nbrs[nbr]);
                        state.U21po_out = ops_z.invL11.matprod(state.U21po_out);    
                        
                        state.U12po_out.reshape(1,state.U12po_out.data_vec.size());
                        state.U21po_out.reshape(1,state.U21po_out.data_vec.size());
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12po_out for self_element_id:{}. Dimensions of ops_z.TB21: {}x{}, state.U12po_out: {}x{}, face_conn:{}", global_id, ops_z.TB21.rows, ops_z.TB21.cols, state.U12po_out.rows, state.U12po_out.cols, face_connections_nbrs[nbr]);
                        state.U12po_out = ops_z.TB21.matprod(state.U12po_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12po_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.U12po_out: {}x{}, face_conn:{}", global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.U12po_out.rows, state.U12po_out.cols, face_connections_nbrs[nbr]);
                        state.U12po_out = ops_z.invL22.matprod(state.U12po_out);

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21po_out for self_element_id:{}. Dimensions of ops_z.TB12: {}x{}, state.U21po_out: {}x{}, face_conn:{}", global_id, ops_z.TB12.rows, ops_z.TB12.cols, state.U21po_out.rows, state.U21po_out.cols, face_connections_nbrs[nbr]);
                        state.U21po_out = ops_z.TB12.matprod(state.U21po_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21po_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.U21po_out: {}x{}, face_conn:{}", global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.U21po_out.rows, state.U21po_out.cols, face_connections_nbrs[nbr]);
                        state.U21po_out = ops_z.invL11.matprod(state.U21po_out);

                        state.U12po_out.reshape(1,state.U12po_out.data_vec.size());
                        state.U21po_out.reshape(1,state.U21po_out.data_vec.size());
                    }
                    break;
                case 2:
                    // using ops_x.TA* because we are getting from bottom neighbor (x dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12om_out for self_element_id:{}. Dimensions of ops_x.TA12: {}x{}, state.U12om_out: {}x{}, face_conn:{}", global_id, ops_x.TA12.rows, ops_x.TA12.cols, state.U12om_out.rows, state.U12om_out.cols, face_connections_nbrs[nbr]);
                        state.U12om_out = ops_x.TA12.matprod(state.U12om_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12om_out for self_element_id:{}. Dimensions of ops_x.invL11: {}x{}, state.U12om_out: {}x{}, face_conn:{}", global_id, ops_x.invL11.rows, ops_x.invL11.cols, state.U12om_out.rows, state.U12om_out.cols, face_connections_nbrs[nbr]);
                        state.U12om_out = ops_x.invL11.matprod(state.U12om_out);

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21om_out for self_element_id:{}. Dimensions of ops_x.TA21: {}x{}, state.U21om_out: {}x{}, face_conn:{}", global_id, ops_x.TA21.rows, ops_x.TA21.cols, state.U21om_out.rows, state.U21om_out.cols, face_connections_nbrs[nbr]);
                        state.U21om_out = ops_x.TA21.matprod(state.U21om_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21om_out for self_element_id:{}. Dimensions of ops_x.invL22: {}x{}, state.U21om_out: {}x{}, face_conn:{}", global_id, ops_x.invL22.rows, ops_x.invL22.cols, state.U21om_out.rows, state.U21om_out.cols, face_connections_nbrs[nbr]);
                        state.U21om_out = ops_x.invL22.matprod(state.U21om_out);        

                        state.U12om_out.reshape(state.U12om_out.data_vec.size(), 1);
                        state.U21om_out.reshape(state.U21om_out.data_vec.size(), 1);                      
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12om_out for self_element_id:{}. Dimensions of ops_x.TA11: {}x{}, state.U12om_out: {}x{}, face_conn:{}", global_id, ops_x.TA11.rows, ops_x.TA11.cols, state.U12om_out.rows, state.U12om_out.cols, face_connections_nbrs[nbr]);
                        state.U12om_out = ops_x.TA11.matprod(state.U12om_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12om_out for self_element_id:{}. Dimensions of ops_x.invL11: {}x{}, state.U12om_out: {}x{}, face_conn:{}", global_id, ops_x.invL11.rows, ops_x.invL11.cols, state.U12om_out.rows, state.U12om_out.cols, face_connections_nbrs[nbr]);
                        state.U12om_out = ops_x.invL11.matprod(state.U12om_out);

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21om_out for self_element_id:{}. Dimensions of ops_x.TA22: {}x{}, state.U21om_out: {}x{}, face_conn:{}", global_id, ops_x.TA22.rows, ops_x.TA22.cols, state.U21om_out.rows, state.U21om_out.cols, face_connections_nbrs[nbr]);
                        state.U21om_out = ops_x.TA22.matprod(state.U21om_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21om_out for self_element_id:{}. Dimensions of ops_x.invL22: {}x{}, state.U21om_out: {}x{}, face_conn:{}", global_id, ops_x.invL22.rows, ops_x.invL22.cols, state.U21om_out.rows, state.U21om_out.cols, face_connections_nbrs[nbr]);
                        state.U21om_out = ops_x.invL22.matprod(state.U21om_out);

                        state.U12om_out.reshape(state.U12om_out.data_vec.size(), 1);
                        state.U21om_out.reshape(state.U21om_out.data_vec.size(), 1);
                    }
                    break;
                case 3:
                    // using ops_x.TA* because we are getting from top neighbor (x dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12op_out for self_element_id:{}. Dimensions of ops_x.TB12: {}x{}, state.U12op_out: {}x{}, face_conn:{}", global_id, ops_x.TB12.rows, ops_x.TB12.cols, state.U12op_out.rows, state.U12op_out.cols, face_connections_nbrs[nbr]);
                        state.U12op_out = ops_x.TB12.matprod(state.U12op_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12op_out for self_element_id:{}. Dimensions of ops_x.invL11: {}x{}, state.U12op_out: {}x{}, face_conn:{}", global_id, ops_x.invL11.rows, ops_x.invL11.cols, state.U12op_out.rows, state.U12op_out.cols, face_connections_nbrs[nbr]);
                        state.U12op_out = ops_x.invL11.matprod(state.U12op_out);

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21op_out for self_element_id:{}. Dimensions of ops_x.TB21: {}x{}, state.U21op_out: {}x{}, face_conn:{}", global_id, ops_x.TB21.rows, ops_x.TB21.cols, state.U21op_out.rows, state.U21op_out.cols, face_connections_nbrs[nbr]);
                        state.U21op_out = ops_x.TB21.matprod(state.U21op_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21op_out for self_element_id:{}. Dimensions of ops_x.invL22: {}x{}, state.U21op_out: {}x{}, face_conn:{}", global_id, ops_x.invL22.rows, ops_x.invL22.cols, state.U21op_out.rows, state.U21op_out.cols, face_connections_nbrs[nbr]);
                        state.U21op_out = ops_x.invL22.matprod(state.U21op_out);

                        state.U12op_out.reshape(state.U12op_out.data_vec.size(), 1);
                        state.U21op_out.reshape(state.U21op_out.data_vec.size(), 1);
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){

                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12op_out for self_element_id:{}. Dimensions of ops_x.TB11: {}x{}, state.U12op_out: {}x{}, face_conn:{}", global_id, ops_x.TB11.rows, ops_x.TB11.cols, state.U12op_out.rows, state.U12op_out.cols, face_connections_nbrs[nbr]);
                        state.U12op_out = ops_x.TB11.matprod(state.U12op_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U12op_out for self_element_id:{}. Dimensions of ops_x.invL11: {}x{}, state.U12op_out: {}x{}, face_conn:{}", global_id, ops_x.invL11.rows, ops_x.invL11.cols, state.U12op_out.rows, state.U12op_out.cols, face_connections_nbrs[nbr]);
                        state.U12op_out = ops_x.invL11.matprod(state.U12op_out);
                        
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21op_out for self_element_id:{}. Dimensions of ops_x.TB22: {}x{}, state.U21op_out: {}x{}, face_conn:{}", global_id, ops_x.TB22.rows, ops_x.TB22.cols, state.U21op_out.rows, state.U21op_out.cols, face_connections_nbrs[nbr]);
                        state.U21op_out = ops_x.TB22.matprod(state.U21op_out);
                        logger->debug("DFDM::Element2D::receive_Uboundary:: applying transformation to U21op_out for self_element_id:{}. Dimensions of ops_x.invL22: {}x{}, state.U21op_out: {}x{}, face_conn:{}", global_id, ops_x.invL22.rows, ops_x.invL22.cols, state.U21op_out.rows, state.U21op_out.cols, face_connections_nbrs[nbr]);
                        state.U21op_out = ops_x.invL22.matprod(state.U21op_out);

                        state.U12op_out.reshape(state.U12op_out.data_vec.size(), 1);
                        state.U21op_out.reshape(state.U21op_out.data_vec.size(), 1);

                    }
                    break;
            }
        }
    }
}

void DFDM::Element2D::apply_Uboundary(uint64_t time_step){
    // all elements have received boundary now apply to the state
    logger->debug("DFDM::Element2D::apply_Uboundary:: applying boundary to U12mo for self_element_id:{}. Dimensions of state.U12mo_inn: {}x{}, state.U12mo_out: {}x{}, state.U21mo_inn: {}x{}, state.U21mo_out: {}x{}", global_id, state.U12mo_inn.rows, state.U12mo_inn.cols, state.U12mo_out.rows, state.U12mo_out.cols, state.U21mo_inn.rows, state.U21mo_inn.cols, state.U21mo_out.rows, state.U21mo_out.cols);
    state.U12mo = state.U12mo_inn * (1- alpha_mo) + state.U12mo_out * alpha_mo;
    state.U21mo = state.U21mo_inn * (1- alpha_mo) + state.U21mo_out * alpha_mo;
    logger->debug("DFDM::Element2D::apply_Uboundary:: applying boundary to U12po for self_element_id:{}. Dimensions of state.U12po_inn: {}x{}, state.U12po_out: {}x{}, state.U21po_inn: {}x{}, state.U21po_out: {}x{}", global_id, state.U12po_inn.rows, state.U12po_inn.cols, state.U12po_out.rows, state.U12po_out.cols, state.U21po_inn.rows, state.U21po_inn.cols, state.U21po_out.rows, state.U21po_out.cols);
    state.U12po = state.U12po_inn * (1- alpha_po) + state.U12po_out * alpha_po;
    state.U21po = state.U21po_inn * (1- alpha_po) + state.U21po_out * alpha_po;
    logger->debug("DFDM::Element2D::apply_Uboundary:: applying boundary to U12om for self_element_id:{}. Dimensions of state.U12om_inn: {}x{}, state.U12om_out: {}x{}, state.U21om_inn: {}x{}, state.U21om_out: {}x{}", global_id, state.U12om_inn.rows, state.U12om_inn.cols, state.U12om_out.rows, state.U12om_out.cols, state.U21om_inn.rows, state.U21om_inn.cols, state.U21om_out.rows, state.U21om_out.cols);
    state.U12om = state.U12om_inn * (1- alpha_om) + state.U12om_out * alpha_om;
    state.U21om = state.U21om_inn * (1- alpha_om) + state.U21om_out * alpha_om;
    logger->debug("DFDM::Element2D::apply_Uboundary:: applying boundary to U12op for self_element_id:{}. Dimensions of state.U12op_inn: {}x{}, state.U12op_out: {}x{}, state.U21op_inn: {}x{}, state.U21op_out: {}x{}", global_id, state.U12op_inn.rows, state.U12op_inn.cols, state.U12op_out.rows, state.U12op_out.cols, state.U21op_inn.rows, state.U21op_inn.cols, state.U21op_out.rows, state.U21op_out.cols);
    state.U12op = state.U12op_inn * (1- alpha_op) + state.U12op_out * alpha_op;
    state.U21op = state.U21op_inn * (1- alpha_op) + state.U21op_out * alpha_op;
}

void DFDM::Element2D::compute_KU2dA(uint64_t time_step){


    auto dUdxp11 = DFDM::utils::dFdxp(ops_x.KK12, state.U21, ops_x.invL22_t, ops_x.invL11, state.U21mo, state.U21po);
    auto dUdxp22 = DFDM::utils::dFdxp(ops_x.KK21, state.U12, ops_x.invL11_t, ops_x.invL22, state.U12mo, state.U12po);

    auto dUdzp11 = DFDM::utils::dFdzp(ops_z.KK12, state.U12, ops_z.invL22_t, ops_z.invL11, state.U12om, state.U12op);
    auto dUdzp22 = DFDM::utils::dFdzp(ops_z.KK21, state.U21, ops_z.invL11_t, ops_z.invL22, state.U21om, state.U21op);


    auto dUdx11 = DFDM::matrix<double>::elementwise_mult(dUdxp11, element_grid.dxpdx11) + 
                  DFDM::matrix<double>::elementwise_mult(dUdzp11, element_grid.dzpdx11);
    auto dUdx22 = DFDM::matrix<double>::elementwise_mult(dUdxp22, element_grid.dxpdx22) + 
                  DFDM::matrix<double>::elementwise_mult(dUdzp22, element_grid.dzpdx22);

    auto dUdz11 = DFDM::matrix<double>::elementwise_mult(dUdxp11, element_grid.dxpdz11) + 
                  DFDM::matrix<double>::elementwise_mult(dUdzp11, element_grid.dzpdz11);
    auto dUdz22 = DFDM::matrix<double>::elementwise_mult(dUdxp22, element_grid.dxpdz22) + 
                  DFDM::matrix<double>::elementwise_mult(dUdzp22, element_grid.dzpdz22);

    // txx
    auto txx11 = DFDM::matrix<double>::elementwise_mult(dUdx11, element_grid.mu11);
    auto txx22 = DFDM::matrix<double>::elementwise_mult(dUdx22, element_grid.mu22);
    // tzz
    auto tzz11 = DFDM::matrix<double>::elementwise_mult(dUdz11, element_grid.mu11);
    auto tzz22 = DFDM::matrix<double>::elementwise_mult(dUdz22, element_grid.mu22);

    // // txx
    // // std::cout << "mu in practice: " << element_grid.mu11(0,0) << std::endl;
    // auto txx11 = dUdx11 * element_grid.mu11(0,0);
    // auto txx22 = dUdx22 * element_grid.mu22(0,0);
    // // tzz
    // auto tzz11 = dUdz11 * element_grid.mu11(0,0);
    // auto tzz22 = dUdz22 * element_grid.mu22(0,0);

    auto Sxx11_add = DFDM::matrix<double>::elementwise_mult(txx11, element_grid.dxpdx11) +
                     DFDM::matrix<double>::elementwise_mult(tzz11, element_grid.dxpdz11);
    state.Sxx11 = DFDM::matrix<double>::elementwise_mult(Sxx11_add, element_grid.Jac11);

    auto Szz11_add = DFDM::matrix<double>::elementwise_mult(txx11, element_grid.dzpdx11) +
                     DFDM::matrix<double>::elementwise_mult(tzz11, element_grid.dzpdz11);
    state.Szz11 = DFDM::matrix<double>::elementwise_mult(Szz11_add, element_grid.Jac11);

    auto Sxx22_add = DFDM::matrix<double>::elementwise_mult(txx22, element_grid.dxpdx22) +
                     DFDM::matrix<double>::elementwise_mult(tzz22, element_grid.dxpdz22);
    state.Sxx22 = DFDM::matrix<double>::elementwise_mult(Sxx22_add, element_grid.Jac22);

    auto Szz22_add = DFDM::matrix<double>::elementwise_mult(txx22, element_grid.dzpdx22) +
                     DFDM::matrix<double>::elementwise_mult(tzz22, element_grid.dzpdz22);
    state.Szz22 = DFDM::matrix<double>::elementwise_mult(Szz22_add, element_grid.Jac22);
}

void DFDM::Element2D::compute_Sboundary_local(){
    compute_boundary_Svalue_inn_elem(ops_x.invL11_t, ops_z.invL11_t, state.Sxx11, state.Sxx11mo_inn, state.Sxx11po_inn, state.Sxx11om_inn, state.Sxx11op_inn);
    compute_boundary_Svalue_inn_elem(ops_x.invL22_t, ops_z.invL22_t, state.Sxx22, state.Sxx22mo_inn, state.Sxx22po_inn, state.Sxx22om_inn, state.Sxx22op_inn);
    compute_boundary_Svalue_inn_elem(ops_x.invL11_t, ops_z.invL11_t, state.Szz11, state.Szz11mo_inn, state.Szz11po_inn, state.Szz11om_inn, state.Szz11op_inn);
    compute_boundary_Svalue_inn_elem(ops_x.invL22_t, ops_z.invL22_t, state.Szz22, state.Szz22mo_inn, state.Szz22po_inn, state.Szz22om_inn, state.Szz22op_inn);
}


void DFDM::Element2D::Sboundary_rotate(){
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){
            switch(nbr){
                case 0:
                    rotate_sij(state.Sxx11mo_inn, state.Szz11mo_inn, state.Sxx11mo_inn_r, state.Szz11mo_inn_r, rot_mat_mo);
                    rotate_sij(state.Sxx22mo_inn, state.Szz22mo_inn, state.Sxx22mo_inn_r, state.Szz22mo_inn_r, rot_mat_mo);
                    break;
                case 1:
                    rotate_sij(state.Sxx11po_inn, state.Szz11po_inn, state.Sxx11po_inn_r, state.Szz11po_inn_r, rot_mat_po);
                    rotate_sij(state.Sxx22po_inn, state.Szz22po_inn, state.Sxx22po_inn_r, state.Szz22po_inn_r, rot_mat_po);
                    break;
                case 2:
                    rotate_sij(state.Sxx11om_inn, state.Szz11om_inn, state.Sxx11om_inn_r, state.Szz11om_inn_r, rot_mat_om);
                    rotate_sij(state.Sxx22om_inn, state.Szz22om_inn, state.Sxx22om_inn_r, state.Szz22om_inn_r, rot_mat_om);
                    break;
                case 3:
                    rotate_sij(state.Sxx11op_inn, state.Szz11op_inn, state.Sxx11op_inn_r, state.Szz11op_inn_r, rot_mat_op);
                    rotate_sij(state.Sxx22op_inn, state.Szz22op_inn, state.Sxx22op_inn_r, state.Szz22op_inn_r, rot_mat_op);
                    break;
            }
        }
        
    }
}

void DFDM::Element2D::compute_boundary_Svalue_inn_elem(DFDM::matrix<double>& xinvL_t, DFDM::matrix<double>& zinvL_t, DFDM::matrix<double>& S_, DFDM::matrix<double>& Smo, DFDM::matrix<double>& Spo, DFDM::matrix<double>& Som, DFDM::matrix<double>& Sop){
    auto tmp = xinvL_t.matprod(S_);
    auto row = tmp.row_get(0);
    Smo.row_set(0, row);
    // Smo = tmp.row_get(0);
    row = tmp.row_get(tmp.rows-1);
    Spo.row_set(0, row);
    // Spo = tmp.row_get(tmp.rows-1);
    tmp = S_.transpose_inplace();
    tmp = zinvL_t.matprod(tmp);
    tmp = tmp.transpose_inplace();
    auto col = tmp.col_get(0);
    Som.col_set(0, col);
    // Som = tmp.col_get(0);
    col = tmp.col_get(tmp.cols-1);
    Sop.col_set(0, col);
    // Sop = tmp.col_get(tmp.cols-1);
}

void DFDM::Element2D::send_Sboundary(std::vector<DFDM::Element2D>& loc_elements){
// make sure to consider the difference between local and remote elements.
// only send if its remote
// loop over al the neighbors, prepare boundary for sending, if it is local 
// then copy it locally in the neighbor element directly

    std::vector<DFDM::matrix<double>> to_send_S11(4), to_send_S22(4);
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){// this determines which face of nbr's am I facing
                to_send_S11[nbr] = get_face(state.Sxx11mo_inn_r, state.Sxx11po_inn_r, state.Sxx11om_inn_r, state.Sxx11op_inn_r, face_connections_self[nbr], rotation_self[nbr]);
                to_send_S22[nbr] = get_face(state.Sxx22mo_inn_r, state.Sxx22po_inn_r, state.Sxx22om_inn_r, state.Sxx22op_inn_r, face_connections_self[nbr], rotation_self[nbr]);

                if(face_connections_self[nbr] == 0 || face_connections_self[nbr] == 1){
                    to_send_S11[nbr] = to_send_S11[nbr].transpose_inplace();
                    to_send_S11[nbr] = ops_z.invL11_t.matprod(to_send_S11[nbr]);
                    to_send_S22[nbr] = to_send_S22[nbr].transpose_inplace();
                    to_send_S22[nbr] = ops_z.invL22_t.matprod(to_send_S22[nbr]);
                }else if(face_connections_self[nbr] == 2 || face_connections_self[nbr] == 3){
                    auto temp = to_send_S11[nbr];
                    to_send_S11[nbr] = ops_x.invL22_t.matprod(to_send_S22[nbr]);
                    to_send_S22[nbr] = ops_x.invL11_t.matprod(temp);
                }else{
                    std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                    ": Invalid face connection when preparing S boundary for sending " << std::endl;
                    exit(-1);
                }
            }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                  to_send_S11[nbr] = get_face(state.Szz11mo_inn_r, state.Szz11po_inn_r, state.Szz11om_inn_r, state.Szz11op_inn_r, face_connections_self[nbr], rotation_self[nbr]);
                  to_send_S22[nbr] = get_face(state.Szz22mo_inn_r, state.Szz22po_inn_r, state.Szz22om_inn_r, state.Szz22op_inn_r, face_connections_self[nbr], rotation_self[nbr]);
                  logger->debug("DFDM::Element2D::send_Sboundary: element:{} extracted to_send_S11({}x{}) and to_send_S22({}x{}) for nbr:{}",
                    global_id, to_send_S11[nbr].rows, to_send_S11[nbr].cols, to_send_S22[nbr].rows, to_send_S22[nbr].cols, neighbors[nbr]);
                  if(face_connections_self[nbr] == 0 || face_connections_self[nbr] == 1){
                    auto temp11 = to_send_S11[nbr];
                    auto temp22 = to_send_S22[nbr];
                    to_send_S11[nbr] = to_send_S22[nbr].transpose_inplace();
                    to_send_S11[nbr] = ops_z.invL22_t.matprod(to_send_S11[nbr]);
                    
                    to_send_S22[nbr] = temp11.transpose_inplace();
                    to_send_S22[nbr] = ops_z.invL11_t.matprod(to_send_S22[nbr]);
                   }else if(face_connections_self[nbr] == 2 || face_connections_self[nbr] == 3){
                    to_send_S11[nbr] = ops_x.invL11_t.matprod(to_send_S11[nbr]);
                    to_send_S22[nbr] = ops_x.invL22_t.matprod(to_send_S22[nbr]);
                   }else{
                    std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                    ": Invalid face connection when preparing S boundary for sending " << std::endl;
                    exit(-1);
                  }
            }else{
                std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                ": Invalid face connection when preparing S boundary for sending " << std::endl;
                exit(-1);
            }              
            if(locality[nbr] == 1){ // if its local to this rank
                // get the element
                DFDM::Element2D& target_el = get_local_element(neighbors[nbr], loc_elements);
                logger->debug("DFDM::Element2D::send_Sboundary: nbr:{} for element:{} is local, nbr's face:{}", neighbors[nbr], global_id, face_connections_nbrs[nbr]); 
                        
                switch(face_connections_nbrs[nbr]){// which of nbr's face am I facing?
                    case 0:
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S11 ({}x{}) to state.Sxx11mo_out (from:{}, to:{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Sxx11mo_out = to_send_S11[nbr];
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S22 ({}x{}) to state.Sxx22mo_out (from:{}, to:{})",to_send_S22[nbr].rows, to_send_S22[nbr].cols, global_id, target_el.global_id);
                        target_el.state.Sxx22mo_out = to_send_S22[nbr];
                        break;
                    case 1:
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S11 ({}x{}) to state.Sxx11po_out (from:{}, to:{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Sxx11po_out = to_send_S11[nbr];
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S11 ({}x{}) to state.Sxx22po_out (from:{}, to:{})",to_send_S22[nbr].rows, to_send_S22[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Sxx22po_out = to_send_S22[nbr];
                        break;
                    case 2:
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S11 ({}x{}) to state.Szz11om_out (from:{}, to:{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Szz11om_out = to_send_S11[nbr];
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S22 ({}x{}) to state.Szz22om_out (from:{}, to:{})",to_send_S22[nbr].rows, to_send_S22[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Szz22om_out = to_send_S22[nbr];
                        break;
                    case 3:
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S11 ({}x{}) to state.Szz11op_out (from:{}, to:{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Szz11op_out = to_send_S11[nbr];
                        logger->debug("DFDM::Element2D::send_Sboundary: copying to_Send_S22 ({}x{}) to state.Szz22op_out (from:{}, to:{})",to_send_S22[nbr].rows, to_send_S22[nbr].cols, global_id, target_el.global_id); 
                        target_el.state.Szz22op_out = to_send_S22[nbr];
                        break;
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when preparing U boundary for sending " << std::endl;
                        exit(-1);
                }
            }else{
                // send to remote
                // when sending generate a unique tag depending on which face is being sent
                // use the target element's face that is facing you to generate the tag
                switch(face_connections_nbrs[nbr]){
                    case 0:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SXX11, DFDM::BdryDirection::MO, neighbor_process[nbr], neighbors[nbr]);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SXX22, DFDM::BdryDirection::MO, neighbor_process[nbr], neighbors[nbr]);
                        logger->debug("DFDM::Element2D::send_Sboundary: element:{} sending SXX11 and SXX22 to nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        logger->debug("DFDM::Element2D::send_Sboundary: sending to_send_S11 ({}x{}) and to_send_S22 ({}x{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, to_send_S22[nbr].rows, to_send_S22[nbr].cols);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS11, to_send_S11[nbr], s_request[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS22, to_send_S22[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        break;
                    }
                    case 1:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SXX11, DFDM::BdryDirection::PO, neighbor_process[nbr], neighbors[nbr]);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SXX22, DFDM::BdryDirection::PO, neighbor_process[nbr], neighbors[nbr]);
                        logger->debug("DFDM::Element2D::send_Sboundary: element:{} sending SXX11 and SXX22 to nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        logger->debug("DFDM::Element2D::send_Sboundary: sending to_send_S11 ({}x{}) and to_send_S22 ({}x{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, to_send_S22[nbr].rows, to_send_S22[nbr].cols);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS11, to_send_S11[nbr], s_request[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS22, to_send_S22[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        break;
                    }
                    case 2:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SZZ11, DFDM::BdryDirection::OM, neighbor_process[nbr], neighbors[nbr]);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SZZ22, DFDM::BdryDirection::OM, neighbor_process[nbr], neighbors[nbr]);
                        logger->debug("DFDM::Element2D::send_Sboundary: element:{} sending SZZ11 and SZZ22 to nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        logger->debug("DFDM::Element2D::send_Sboundary: sending to_send_S11 ({}x{}) and to_send_S22 ({}x{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, to_send_S22[nbr].rows, to_send_S22[nbr].cols);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS11, to_send_S11[nbr], s_request[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS22, to_send_S22[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        break;
                    }
                    case 3:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SZZ11, DFDM::BdryDirection::OP, neighbor_process[nbr], neighbors[nbr]);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SZZ22, DFDM::BdryDirection::OP, neighbor_process[nbr], neighbors[nbr]);
                        logger->debug("DFDM::Element2D::send_Sboundary: element:{} sending SZZ11 and SZZ22 to nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        logger->debug("DFDM::Element2D::send_Sboundary: sending to_send_S11 ({}x{}) and to_send_S22 ({}x{})",to_send_S11[nbr].rows, to_send_S11[nbr].cols, to_send_S22[nbr].rows, to_send_S22[nbr].cols);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS11, to_send_S11[nbr], s_request[nbr]);
                        DFDM::send_matnb_dim(neighbor_process[nbr], tagS22, to_send_S22[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        break;
                    }
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when preparing S boundary for sending " << std::endl;
                        exit(-1);
                }
            }
        }
    }
}   

void DFDM::Element2D::receive_Sboundary(std::vector<DFDM::Element2D>& loc_elements){

    std::vector<DFDM::matrix<double>> to_recv_S11(4);
    std::vector<DFDM::matrix<double>> to_recv_S22(4);
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            if(locality[nbr] == 0){ // if its remote to this rank (1 means local, 0 means remote)
                //then I need to receive the boundary from remote
                // when receiving use your own process id, element id and your face that is facing
                // the neighbor to generate the tag (self)
                // depending on the face connection, assign the received boundary to the correct variable
                switch(face_connections_self[nbr]){
                    case 0:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SXX11, DFDM::BdryDirection::MO, owner_cpu, global_id);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SXX22, DFDM::BdryDirection::MO, owner_cpu, global_id);
                        logger->debug("DFDM::Element2D::receive_Sboundary: element:{} receiving S11 and S22 from nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_S11[nbr].resize(neighbor_Nz1[nbr], 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nz1[nbr] - 1, 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_S11[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nx1[nbr], 1); // 
                        }
                        
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS11, to_recv_S11[nbr], r_request[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS22, to_recv_S22[nbr], r_request[nbr+4]);
                        // DFDM::mpi_wait(r_request[1], r_request[3]);
                        // state.Sxx11mo_out = to_recv_S11[nbr];
                        // state.Sxx22mo_out = to_recv_S22[nbr];
                        break;
                    }
                    case 1:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SXX11, DFDM::BdryDirection::PO, owner_cpu, global_id);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SXX22, DFDM::BdryDirection::PO, owner_cpu, global_id);
                        logger->debug("DFDM::Element2D::receive_Sboundary: element:{} receiving S11 and S22 from nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_S11[nbr].resize(neighbor_Nz1[nbr], 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nz1[nbr] - 1, 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_S11[nbr].resize(neighbor_Nx1[nbr] - 1, 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nx1[nbr], 1); // 
                        }

                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS11, to_recv_S11[nbr], r_request[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS22, to_recv_S22[nbr], r_request[nbr+4]);
                        // DFDM::mpi_wait(r_request[5], r_request[7]);
                        // state.Sxx11po_out = to_recv_S11[nbr];
                        // state.Sxx22po_out = to_recv_S22[nbr];
                        break;
                    }
                    case 2:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SZZ11, DFDM::BdryDirection::OM, owner_cpu, global_id);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SZZ22, DFDM::BdryDirection::OM, owner_cpu, global_id);
                        logger->debug("DFDM::Element2D::receive_Sboundary: element:{} receiving S11 and S22 from nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_S11[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_S11[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nx1[nbr] - 1, 1); // 
                        }

                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS11, to_recv_S11[nbr], r_request[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS22, to_recv_S22[nbr], r_request[nbr+4]);
                        // DFDM::mpi_wait(r_request[9], r_request[11]);
                        // state.Szz11om_out = to_recv_S11[nbr];
                        // state.Szz22om_out = to_recv_S22[nbr];
                        break;
                    }
                    case 3:{
                        auto tagS11 = DFDM::generateTag(DFDM::BdryType::SZZ11, DFDM::BdryDirection::OP, owner_cpu, global_id);
                        auto tagS22 = DFDM::generateTag(DFDM::BdryType::SZZ22, DFDM::BdryDirection::OP, owner_cpu, global_id);
                        logger->debug("DFDM::Element2D::receive_Sboundary: element:{} receiving S11 and S22 from nbr:{} with tagS11:{}, tagS22:{}", global_id, neighbors[nbr], tagS11, tagS22);
                        if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                            to_recv_S11[nbr].resize(neighbor_Nz1[nbr] - 1, 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nz1[nbr], 1); // 
                        }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                            to_recv_S11[nbr].resize(neighbor_Nx1[nbr], 1); // 
                            to_recv_S22[nbr].resize( neighbor_Nx1[nbr] - 1, 1); // 
                        }
                        
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS11, to_recv_S11[nbr], r_request[nbr]);
                        DFDM::recv_matnb_dim(neighbor_process[nbr], tagS22, to_recv_S22[nbr], r_request[nbr+4]);
                        // DFDM::mpi_wait(r_request[13], r_request[15]);
                        // state.Szz11op_out = to_recv_S11[nbr];
                        // state.Szz22op_out = to_recv_S22[nbr];
                        break;
                    }
                }
            }else{
                //if neighbor is local then we already have the boundary in *out variables so do nothing
            }
        }  
    }
    logger->debug("DFDM::Element2D:receive_Sboundary: completed the non-blocking receive loop for element:{}", global_id);
    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            if(locality[nbr] == 0){
                switch(nbr){
                    case 0:
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        state.Sxx11mo_out = to_recv_S11[nbr];
                        state.Sxx22mo_out = to_recv_S22[nbr];
                        break;
                    case 1:
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        state.Sxx11po_out = to_recv_S11[nbr];
                        state.Sxx22po_out = to_recv_S22[nbr];
                        break;
                    case 2:
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        state.Szz11om_out = to_recv_S11[nbr];
                        state.Szz22om_out = to_recv_S22[nbr];
                        break;
                    case 3:
                        // DFDM::mpi_wait(s_request[nbr], s_request[nbr+4]);
                        DFDM::mpi_wait(r_request[nbr], r_request[nbr+4]);
                        state.Szz11op_out = to_recv_S11[nbr];
                        state.Szz22op_out = to_recv_S22[nbr];
                        break;
                    default:
                        std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
                        ": Invalid neighbor index when receiving S boundaries " << std::endl;
                        exit(-1);
                }
            }
        }
    }

    for(uint32_t nbr = 0; nbr < 4; nbr++){
        if(neighbors[nbr] != -1){ // if neighbor exists
            switch(nbr){ // here nbr can also be replaced with face_connections_self[nbr] as they are the same
                case 0:
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        // using ops_z.TA* because we are getting from left neighbor (z dim)
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx11mo_out for self_element_id:{}. Dimensions of ops_z.TA11: {}x{}, state.Sxx11mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.TA11.rows, ops_z.TA11.cols, state.Sxx11mo_out.rows, state.Sxx11mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx11mo_out = ops_z.TA11.matprod(state.Sxx11mo_out);
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx11mo_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.Sxx11mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.Sxx11mo_out.rows, state.Sxx11mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx11mo_out = ops_z.invL11.matprod(state.Sxx11mo_out);
                        
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx22mo_out for self_element_id:{}. Dimensions of ops_z.TA22: {}x{}, state.Sxx22mo_out: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.TA22.rows, ops_z.TA22.cols, state.Sxx22mo_out.rows, state.Sxx22mo_out.cols, face_connections_nbrs[nbr], nbr); 
                        state.Sxx22mo_out = ops_z.TA22.matprod(state.Sxx22mo_out);
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx22mo_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.invL22: {}x{}, face_conn:{}, my_face:{}", global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.Sxx22mo_out.rows, state.Sxx22mo_out.cols, face_connections_nbrs[nbr], nbr); 
                        state.Sxx22mo_out = ops_z.invL22.matprod(state.Sxx22mo_out);

                        state.Sxx11mo_out.reshape(1, state.Sxx11mo_out.data_vec.size());
                        state.Sxx22mo_out.reshape(1, state.Sxx22mo_out.data_vec.size());
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        // using ops_z.TA* because we are getting from left neighbor (z dim)
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx11mo_out for self_element_id:{}. Dimensions of ops_z.TA12: {}x{}, state.Sxx11mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.TA12.rows, ops_z.TA12.cols, state.Sxx11mo_out.rows, state.Sxx11mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx11mo_out = ops_z.TA12.matprod(state.Sxx11mo_out);
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx11mo_out for self_element_id:{}. Dimensions of ops_z.invL11: {}x{}, state.Sxx11mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.invL11.rows, ops_z.invL11.cols, state.Sxx11mo_out.rows, state.Sxx11mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx11mo_out = ops_z.invL11.matprod(state.Sxx11mo_out);

                        // state.U21mo_out = to_recv_U21;
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx11mo_out for self_element_id:{}. Dimensions of ops_z.TA21: {}x{}, state.Sxx22mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.TA21.rows, ops_z.TA21.cols, state.Sxx22mo_out.rows, state.Sxx22mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx22mo_out = ops_z.TA21.matprod(state.Sxx22mo_out);
                        logger->debug("DFDM::Element2D::receive_Sboundary:: applying transformation to Sxx22mo_out for self_element_id:{}. Dimensions of ops_z.invL22: {}x{}, state.Sxx11mo_out: {}x{}, face_conn:{}, my_face:{}", 
                        global_id, ops_z.invL22.rows, ops_z.invL22.cols, state.Sxx22mo_out.rows, state.Sxx22mo_out.cols, face_connections_nbrs[nbr], nbr);
                        state.Sxx22mo_out = ops_z.invL22.matprod(state.Sxx22mo_out);


                        state.Sxx11mo_out.reshape(1, state.Sxx11mo_out.data_vec.size());
                        state.Sxx22mo_out.reshape(1, state.Sxx22mo_out.data_vec.size());
                    }
                    break;
                case 1:
                        // using ops_z.TA* because we are getting from right neighbor (z dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        state.Sxx11po_out = ops_z.TB11.matprod(state.Sxx11po_out);
                        state.Sxx11po_out = ops_z.invL11.matprod(state.Sxx11po_out);

                        state.Sxx22po_out = ops_z.TB22.matprod(state.Sxx22po_out);
                        state.Sxx22po_out = ops_z.invL22.matprod(state.Sxx22po_out);   

                        state.Sxx11po_out.reshape(1, state.Sxx11po_out.data_vec.size());
                        state.Sxx22po_out.reshape(1, state.Sxx22po_out.data_vec.size()); 
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        state.Sxx11po_out = ops_z.TB12.matprod(state.Sxx11po_out);
                        state.Sxx11po_out = ops_z.invL11.matprod(state.Sxx11po_out);

                        state.Sxx22po_out = ops_z.TB21.matprod(state.Sxx22po_out);
                        state.Sxx22po_out = ops_z.invL22.matprod(state.Sxx22po_out);

                        state.Sxx11po_out.reshape(1, state.Sxx11po_out.data_vec.size());
                        state.Sxx22po_out.reshape(1, state.Sxx22po_out.data_vec.size());
                    }
                    break;
                case 2:
                    // using ops_x.TA* because we are getting from bottom neighbor (x dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        state.Szz11om_out = ops_x.TA12.matprod(state.Szz11om_out);
                        state.Szz11om_out = ops_x.invL11.matprod(state.Szz11om_out);
                        
                        state.Szz22om_out = ops_x.TA21.matprod(state.Szz22om_out);
                        state.Szz22om_out = ops_x.invL22.matprod(state.Szz22om_out);

                        state.Szz11om_out.reshape(state.Szz11om_out.data_vec.size(),1);
                        state.Szz22om_out.reshape(state.Szz22om_out.data_vec.size(),1);
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        state.Szz11om_out = ops_x.TA11.matprod(state.Szz11om_out);
                        state.Szz11om_out = ops_x.invL11.matprod(state.Szz11om_out);

                        state.Szz22om_out = ops_x.TA22.matprod(state.Szz22om_out);
                        state.Szz22om_out = ops_x.invL22.matprod(state.Szz22om_out); 

                        state.Szz11om_out.reshape(state.Szz11om_out.data_vec.size(),1);
                        state.Szz22om_out.reshape(state.Szz22om_out.data_vec.size(),1);  
                    }
                    break;
                case 3:
                    // using ops_x.TA* because we are getting from top neighbor (x dim)
                    if(face_connections_nbrs[nbr] == 0 || face_connections_nbrs[nbr] == 1){
                        state.Szz11op_out = ops_x.TB12.matprod(state.Szz11op_out);
                        state.Szz11op_out = ops_x.invL11.matprod(state.Szz11op_out);

                        state.Szz22op_out = ops_x.TB21.matprod(state.Szz22op_out);
                        state.Szz22op_out = ops_x.invL22.matprod(state.Szz22op_out);

                        state.Szz11op_out.reshape(state.Szz11op_out.data_vec.size(),1);
                        state.Szz22op_out.reshape(state.Szz22op_out.data_vec.size(),1);
                    }else if(face_connections_nbrs[nbr] == 2 || face_connections_nbrs[nbr] == 3){
                        state.Szz11op_out = ops_x.TB11.matprod(state.Szz11op_out);
                        state.Szz11op_out = ops_x.invL11.matprod(state.Szz11op_out);
                        
                        state.Szz22op_out = ops_x.TB22.matprod(state.Szz22op_out);
                        state.Szz22op_out = ops_x.invL22.matprod(state.Szz22op_out);

                        state.Szz11op_out.reshape(state.Szz11op_out.data_vec.size(),1);
                        state.Szz22op_out.reshape(state.Szz22op_out.data_vec.size(),1);
                    }
                    break;
            }
        }
    }
        // applying transformation to the boundary and converting basis
}

void DFDM::Element2D::apply_Sboundary(){
    // all elements have received boundary now apply to the state
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Sxx11mo for self_element_id:{}. Dimensions of state.Sxx11mo_inn: {}x{}, state.Sxx11mo_out: {}x{}, alpha_mo: {}", 
    global_id, state.Sxx11mo_inn.rows, state.Sxx11mo_inn.cols, state.Sxx11mo_out.rows, state.Sxx11mo_out.cols, alpha_mo);
    state.Sxx11mo = state.Sxx11mo_inn * alpha_mo + state.Sxx11mo_out * (1-alpha_mo);
    
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Sxx22mo for self_element_id:{}. Dimensions of state.Sxx22mo_inn: {}x{}, state.Sxx22mo_out: {}x{}, alpha_mo: {}", 
    global_id, state.Sxx22mo_inn.rows, state.Sxx22mo_inn.cols, state.Sxx22mo_out.rows, state.Sxx22mo_out.cols, alpha_mo);
    state.Sxx22mo = state.Sxx22mo_inn * alpha_mo + state.Sxx22mo_out * (1-alpha_mo);
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Sxx11po for self_element_id:{}. Dimensions of state.Sxx11po_inn: {}x{}, state.Sxx11po_out: {}x{}, alpha_po: {}", global_id, state.Sxx11po_inn.rows, state.Sxx11po_inn.cols, state.Sxx11po_out.rows, state.Sxx11po_out.cols, alpha_po);

    state.Sxx11po = state.Sxx11po_inn * alpha_po + state.Sxx11po_out * (1-alpha_po);
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Sxx22po for self_element_id:{}. Dimensions of state.Sxx22po_inn: {}x{}, state.Sxx22po_out: {}x{}, alpha_po: {}", global_id, state.Sxx22po_inn.rows, state.Sxx22po_inn.cols, state.Sxx22po_out.rows, state.Sxx22po_out.cols, alpha_po);

    state.Sxx22po = state.Sxx22po_inn * alpha_po + state.Sxx22po_out * (1-alpha_po);
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Szz11om for self_element_id:{}. Dimensions of state.Szz11om_inn: {}x{}, state.Szz11om_out: {}x{}, alpha_om: {}", global_id, state.Szz11om_inn.rows, state.Szz11om_inn.cols, state.Szz11om_out.rows, state.Szz11om_out.cols, alpha_om);

    state.Szz11om = state.Szz11om_inn * alpha_om + state.Szz11om_out * (1-alpha_om);
    
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Szz22om for self_element_id:{}. Dimensions of state.Szz22om_inn: {}x{}, state.Szz22om_out: {}x{}, alpha_om: {}", global_id, state.Szz22om_inn.rows, state.Szz22om_inn.cols, state.Szz22om_out.rows, state.Szz22om_out.cols, alpha_om);
    state.Szz22om = state.Szz22om_inn * alpha_om + state.Szz22om_out * (1-alpha_om);
    
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Szz11op for self_element_id:{}. Dimensions of state.Szz11op_inn: {}x{}, state.Szz11op_out: {}x{}, alpha_op: {}", global_id, state.Szz11op_inn.rows, state.Szz11op_inn.cols, state.Szz11op_out.rows, state.Szz11op_out.cols, alpha_op);
    state.Szz11op = state.Szz11op_inn * alpha_op + state.Szz11op_out * (1-alpha_op);
    
    logger->debug("DFDM::Element2D::apply_Sboundary:: applying boundary to Szz22op for self_element_id:{}. Dimensions of state.Szz22op_inn: {}x{}, state.Szz22op_out: {}x{}, alpha_op: {}", global_id, state.Szz22op_inn.rows, state.Szz22op_inn.cols, state.Szz22op_out.rows, state.Szz22op_out.cols, alpha_op);
    state.Szz22op = state.Szz22op_inn * alpha_op + state.Szz22op_out * (1-alpha_op);
}

void DFDM::Element2D::compute_KS2dA(uint64_t time_step){

    auto dSxxdxp21 = DFDM::utils::dFdxp(ops_x.KK21, state.Sxx11, ops_x.invL11_t, ops_x.invL22, state.Sxx11mo, state.Sxx11po);
    auto dSxxdxp12 = DFDM::utils::dFdxp(ops_x.KK12, state.Sxx22, ops_x.invL22_t, ops_x.invL11, state.Sxx22mo, state.Sxx22po);

    auto dSzzdzp21 = DFDM::utils::dFdzp(ops_z.KK12, state.Szz22, ops_z.invL22_t, ops_z.invL11, state.Szz22om, state.Szz22op);
    auto dSzzdzp12 = DFDM::utils::dFdzp(ops_z.KK21, state.Szz11, ops_z.invL11_t, ops_z.invL22, state.Szz11om, state.Szz11op);


    state.dU2dtt21 = dSxxdxp21 + dSzzdzp21;
    state.dU2dtt12 = dSxxdxp12 + dSzzdzp12;
    

}

void DFDM::Element2D::source_injection(const DFDM::source& source, const DFDM::simulation& sim, uint64_t it){
    // for debugging:
    // state.Umid[it].print_file("Umid_"+std::to_string(it)+"_"+std::to_string(global_id));
    // element_grid.Jac12.print_file("Jac12_"+std::to_string(it)+"_"+std::to_string(global_id));
    // element_grid.Jac21.print_file("Jac21_"+std::to_string(it)+"_"+std::to_string(global_id));
    // element_grid.rho12.print_file("rho12_"+std::to_string(it)+"_"+std::to_string(global_id));
    // element_grid.rho21.print_file("rho21_"+std::to_string(it)+"_"+std::to_string(global_id));
    
    // if I am the element with source, I do the following
    // source will not have been initialized in the process with no source element
    if(is_source(sim)){
        auto som = source.selement_id;
        auto invMsg12 = source.invMsg12;
        auto invMsg21 = source.invMsg21;
        auto Fit = source.Ft(it,0);

        if(global_id == som){
            logger->debug("DFDM::Element2D::source_injection in source element:{}", global_id);
            logger->debug("DFDM::Element2D::source_injection:: dU2dtt12 size: {}x{}, dU2dtt21 size: {}x{}, invMsg12 size: {}x{}, invMsg21 size: {}x{}", state.dU2dtt12.rows, state.dU2dtt12.cols, state.dU2dtt21.rows, state.dU2dtt21.cols, invMsg12.rows, invMsg12.cols, invMsg21.rows, invMsg21.cols);
            // state.dU2dtt12.print_file("./output/before_src_inj_state.dU2dtt12_"+std::to_string(it));
            // source.Ft.print_file("./output/src_Ft_"+std::to_string(it));
            state.dU2dtt12 = state.dU2dtt12 + invMsg12*Fit;
            state.dU2dtt21 = state.dU2dtt21 + invMsg21*Fit;
            // invMsg12.print_file("./output/invMsg12_"+std::to_string(it));
            // state.dU2dtt12.print_file("./output/src_inj_state.dU2dtt12_"+std::to_string(it));
            // state.dU2dtt21.print_file("./output/src_inj_state.dU2dtt21_"+std::to_string(it));
        }else{
            std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
            ": Conflict between source id and source element" << std::endl;
            exit(-1);
        }

    }
    
    logger->debug("DFDM::Element2D::source_injection state.dU2dtt12 size: {}x{}, element_grid.Jac12 size: {}x{}", 
                  state.dU2dtt12.rows, state.dU2dtt12.cols, element_grid.Jac12.rows, element_grid.Jac12.cols);

    auto temp12 = DFDM::matrix<double>::elementwise_div(state.dU2dtt12, element_grid.Jac12);
    logger->debug("DFDM::Element2D::source_injection temp12 size: {}x{}, element_grid.rho12 size: {}x{}", 
                  temp12.rows, temp12.cols, element_grid.rho12.rows, element_grid.rho12.cols);             
    state.dU2dtt12 = DFDM::matrix<double>::elementwise_div(temp12, element_grid.rho12);

    logger->debug("DFDM::Element2D::source_injection state.dU2dtt21 size: {}x{}, element_grid.Jac21 size: {}x{}", 
                  state.dU2dtt21.rows, state.dU2dtt21.cols, element_grid.Jac21.rows, element_grid.Jac21.cols); 
    auto temp21 = DFDM::matrix<double>::elementwise_div(state.dU2dtt21, element_grid.Jac21);
    logger->debug("DFDM::Element2D::source_injection temp21 size: {}x{}, element_grid.rho21 size: {}x{}", 
                  temp21.rows, temp21.cols, element_grid.rho21.rows, element_grid.rho21.cols);
    state.dU2dtt21 = DFDM::matrix<double>::elementwise_div(temp21, element_grid.rho21);
}

bool DFDM::Element2D::is_source(const DFDM::simulation& sim){
    return sim.source_elem_id == global_id;
}

bool DFDM::Element2D::is_receiver(const DFDM::simulation& sim){
    return std::find(sim.receiver_elem_ids.begin(), sim.receiver_elem_ids.end(), global_id) != sim.receiver_elem_ids.end();
}

DFDM::receiver& DFDM::Element2D::get_local_receiver(uint64_t global_id, std::vector<DFDM::receiver>& receivers) {
    for (auto& receiver : receivers) {
        if (receiver.relement_id == global_id) {
            return receiver;
        }
    }
    std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << ": Receiver not found for element ID!" << std::endl;
    exit(-1);
}


void DFDM::Element2D::receiver_record(std::vector<DFDM::receiver>& receivers, const DFDM::simulation& sim, uint64_t it){

    // if I am the element with receiver, I do the following
    // receiver will not have been initialized in the process with no receiver element
    if(is_receiver(sim)){   
        decltype(auto) receiver = get_local_receiver(global_id, receivers);
        auto Ub21 = DFDM::utils::tensorProduct2D(ops_x.invL22_t, ops_z.invL11_t, state.U21);
        auto Ub12 = DFDM::utils::tensorProduct2D(ops_x.invL11_t, ops_z.invL22_t, state.U12);

        auto Umid1 = DFDM::utils::tensorProduct2D(receiver.gbx2, receiver.gbz1, Ub21);
        auto Umid2 = DFDM::utils::tensorProduct2D(receiver.gbx1, receiver.gbz2, Ub12);

        auto Sigxx11 = DFDM::utils::tensorProduct2D(ops_x.invL11_t, ops_z.invL11_t, state.Sxx11);
        auto Sigxx22 = DFDM::utils::tensorProduct2D(ops_x.invL22_t, ops_z.invL22_t, state.Sxx22);

        auto Sigxx1 = DFDM::utils::tensorProduct2D(receiver.gbx1, receiver.gbz1, Sigxx11);
        auto Sigxx2 = DFDM::utils::tensorProduct2D(receiver.gbx2, receiver.gbz2, Sigxx22);

        auto Sigzz11 = DFDM::utils::tensorProduct2D(ops_x.invL11_t, ops_z.invL11_t, state.Szz11);
        auto Sigzz22 = DFDM::utils::tensorProduct2D(ops_x.invL22_t, ops_z.invL22_t, state.Szz22);

        auto Sigzz1 = DFDM::utils::tensorProduct2D(receiver.gbx1, receiver.gbz1, Sigzz11);
        auto Sigzz2 = DFDM::utils::tensorProduct2D(receiver.gbx2, receiver.gbz2, Sigzz22);


        if (Umid2.rows != 1 || Umid2.cols != 1) {
            std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << 
            ": Umid2 is not a 1x1 matrix" << std::endl;
            exit(-1);
        }else{
            receiver.ur.value_set(it, 1, Umid1.data_vec[0]);
            receiver.ur.value_set(it, 2, Umid2.data_vec[0]);

            receiver.sigr.value_set(it, 1, Sigxx1.data_vec[0]);
            receiver.sigr.value_set(it, 2, Sigxx2.data_vec[0]);

            receiver.sigr.value_set(it, 3, Sigzz1.data_vec[0]);
            receiver.sigr.value_set(it, 4, Sigzz2.data_vec[0]);
        }
    }
}

void DFDM::Element2D::wavefield_record(uint64_t it){

    auto Ub21 = DFDM::utils::tensorProduct2D(ops_x.invL22_t, ops_z.invL11_t, state.U21);
    auto Ub12 = DFDM::utils::tensorProduct2D(ops_x.invL11_t, ops_z.invL22_t, state.U12);

    auto Umid1 = DFDM::utils::tensorProduct2D(ops_x.B2_t, ops_z.B1_t, Ub21);
    auto Umid2 = DFDM::utils::tensorProduct2D(ops_x.B1_t, ops_z.B2_t, Ub12);

    state.Umid[it] = (Umid1 + Umid2) * 0.5;
}

void DFDM::Element2D::rotate_sij(const DFDM::matrix<double>& sxx, const DFDM::matrix<double>& szz, DFDM::matrix<double>& sxx_r, DFDM::matrix<double>& szz_r, const DFDM::matrix<double>& rot_mat){
    DFDM::matrix<double> sij(1, 2);
    sxx_r.resize(sxx.rows, sxx.cols);
    szz_r.resize(szz.rows, szz.cols);
    logger->debug("DFDM::Element2D::rotate_sij: sxx.rows: {}, sxx.cols: {}, szz.rows: {}, szz.cols: {}", sxx.rows, sxx.cols, szz.rows, szz.cols);
    for (uint32_t k = 0; k < sxx.cols; k++) {
        for (uint32_t i = 0; i < sxx.rows; i++) {
            // Assign the stress components to the sij matrix
            sij.value_set(0, 0, sxx(i, k));
            sij.value_set(0, 1, szz(i, k));
            // Perform the rotation using matrix multiplication
            sij = sij.matprod(rot_mat);
            // Assign the rotated stress components back to the original matrices
            sxx_r.value_set(i, k, sij(0, 0));
            szz_r.value_set(i, k, sij(0, 1));
        }
    }
}


// flips the boundary in 1D
DFDM::matrix<double> DFDM::Element2D::flip1D_boundary(DFDM::matrix<double>& face, uint32_t i1, uint32_t i2, uint32_t j1, uint32_t j2, bool orient){
    DFDM::matrix<double> R(i2-i1+1, j2-j1+1);
    if (orient == 1){
        for (uint32_t i = i1; i <= i2; i++) {
            for (uint32_t j = j1; j <= j2; j++) {
                R.value_set(i2-i, j2-j, face(i, j));
            }
        }
    }else{
        R = face;
    }
    return R;
}

DFDM::matrix<double> DFDM::Element2D::get_face(DFDM::matrix<double>& Umo, DFDM::matrix<double>& Upo, DFDM::matrix<double>& Uom, DFDM::matrix<double>& Uop, uint32_t jFace, bool orient){
    DFDM::matrix<double> Uout;
    switch (jFace){
        case 0:
            Uout = flip1D_boundary(Umo, 0, Umo.rows-1, 0, Umo.cols-1, orient);
            break;
        case 1:
            Uout = flip1D_boundary(Upo, 0, Upo.rows-1, 0, Upo.cols-1, orient);
            break;
        case 2:
            Uout = flip1D_boundary(Uom, 0, Uom.rows-1, 0, Uom.cols-1, orient);
            break;
        case 3:
            Uout = flip1D_boundary(Uop, 0, Uop.rows-1, 0, Uop.cols-1, orient);
            break;
        default:
            std::cerr << "Error: Invalid face index" << std::endl;
            exit(-1);
    }
    return Uout;
}

void DFDM::Element2D::gen_rotation_matrices(){
    rot_mat_mo.resize(2,2);
    rot_mat_po.resize(2,2);
    rot_mat_om.resize(2,2);
    rot_mat_op.resize(2,2);

    for(int jFace = 0; jFace < 4; jFace++){
        switch(jFace){
            case 0:{
                int iNbr_left = neighbors[jFace];
                int iFace_left = face_connections_nbrs[jFace];
                bool rot_left = rotation_nbr[jFace];
                int sign_left = rot_left == 1? -1 : 1;
                if(iNbr_left != -1){
                    switch(iFace_left){
                        case 0:
                            rot_mat_mo.value_set(0, 0, -1);
                            rot_mat_mo.value_set(1, 1, 1);
                            break;
                        case 1:
                            rot_mat_mo.value_set(0, 0, 1);
                            rot_mat_mo.value_set(1, 1, sign_left);
                            break;
                        case 2:
                            rot_mat_mo.value_set(0, 1, -1);
                            rot_mat_mo.value_set(1, 0, -1);
                            break;
                        case 3:
                            rot_mat_mo.value_set(0, 1, 1);
                            rot_mat_mo.value_set(1, 0, 1);
                            break;
                    }
                }
                break;
            }
            case 1:{
                int iNbr_right = neighbors[jFace];
                int iFace_right = face_connections_nbrs[jFace];
                bool rot_right = rotation_nbr[jFace];
                int sign_right = rot_right == 1? -1 : 1;
                if(iNbr_right != -1){
                    switch(iFace_right){
                        case 0:
                            rot_mat_po.value_set(0, 0, 1);
                            rot_mat_po.value_set(1, 1, sign_right);
                            break;
                        case 1:
                            rot_mat_po.value_set(0, 0, -1);
                            rot_mat_po.value_set(1, 1, 1);
                            break;
                        case 2:
                            rot_mat_po.value_set(0, 1, 1);
                            rot_mat_po.value_set(1, 0, 1);
                            break;
                        case 3:
                            rot_mat_po.value_set(0, 1, -1);
                            rot_mat_po.value_set(1, 0, -1);
                            break;
                    }
                }
                break;
            }
            case 2:{
                int iNbr_bottom = neighbors[jFace];
                int iFace_bottom = face_connections_nbrs[jFace];
                bool rot_bottom = rotation_nbr[jFace];
                int sign_bottom = rot_bottom == 1? -1 : 1;
                if(iNbr_bottom != -1){
                    switch(iFace_bottom){
                        case 0:
                            rot_mat_om.value_set(1, 0, -1);
                            rot_mat_om.value_set(0, 1, -1);
                            break;
                        case 1:
                            rot_mat_om.value_set(1, 0, 1);
                            rot_mat_om.value_set(0, 1, 1);
                            break;
                        case 2:
                            rot_mat_om.value_set(1, 1, -1);
                            rot_mat_om.value_set(0, 0, 1);
                            break;
                        case 3:
                            rot_mat_om.value_set(1, 1, 1);
                            rot_mat_om.value_set(0, 0, sign_bottom);
                            break;
                    }
                }
                break;
            }
            case 3:{
                int iNbr_top = neighbors[jFace];
                int iFace_top = face_connections_nbrs[jFace];
                bool rot_top = rotation_nbr[jFace];
                int sign_top = rot_top == 1? -1 : 1;
                if(iNbr_top != -1){
                    switch(iFace_top){
                        case 0:
                            rot_mat_op.value_set(1, 0, 1);
                            rot_mat_op.value_set(0, 1, 1);
                            break;
                        case 1:
                            rot_mat_op.value_set(1, 0, -1);
                            rot_mat_op.value_set(0, 1, -1);
                            break;
                        case 2:
                            rot_mat_op.value_set(1, 1, 1);
                            rot_mat_op.value_set(0, 0, sign_top);
                            break;
                        case 3:
                            rot_mat_op.value_set(1, 1, -1);
                            rot_mat_op.value_set(0, 0, sign_top);
                            break;
                    }
                }
                break;
            }
        }
    }

    // rot_mat_mo.print_file("rot_mat_mo_"+std::to_string(global_id));
    // rot_mat_po.print_file("rot_mat_po"+std::to_string(global_id));
    // rot_mat_om.print_file("rot_mat_om"+std::to_string(global_id));
    // rot_mat_op.print_file("rot_mat_op"+std::to_string(global_id));
}

bool DFDM::Element2D::is_local_element(uint64_t global_id, std::vector<DFDM::Element2D>& elements){
    for (auto& element : elements){
        if (element.global_id == global_id){
            return true;
        }
    }
    return false;
}

DFDM::Element2D& DFDM::Element2D::get_local_element(uint64_t global_id, std::vector<DFDM::Element2D>& elements){
    for (auto& element : elements) {
        if (element.global_id == global_id) {
            return element;
        }
    }
    std::cerr << "Error in " << __FILE__ << " at line " << __LINE__ << ": Element not found locally!" << std::endl;
    exit(-1);
}

void DFDM::Element2D::update_dt_nt(const std::vector<DFDM::Element2D>& elements, DFDM::simulation sim, std::shared_ptr<spdlog::logger> logger_){
    double CFL_condition = 0.1;
    // sim.delta_t = 9999999.9; // initialize to a large value
    // for(uint64_t iom = 0; iom < elements.size(); iom++){
    //     uint64_t Nx1 = elements[iom].Nx1;
    //     uint64_t Nz1 = elements[iom].Nz1;
    //     uint64_t ix = std::ceil((double)(Nx1+1)/2);
    //     uint64_t iz = std::ceil((double)(Nz1+1)/2);
    //     // note that grid_x and grid_z are matrices in the form of (Nx1+1) x (Nz1+1)
    //     double x1 = elements[iom].element_grid.x2d11(ix-1,iz-1); // fixing for matlab indexing 
    //     double z1 = elements[iom].element_grid.z2d11(ix-1,iz-1);

    //     double x2 = elements[iom].element_grid.x2d11(ix,iz);
    //     double z2 = elements[iom].element_grid.z2d11(ix,iz);

    //     double dxz = std::sqrt(std::pow(x2-x1,2) + std::pow(z2-z1,2));

    //     double tmp = CFL_condition*dxz/sim.Vmax;
    //     sim.delta_t = std::min(sim.delta_t, tmp);
    // }
    // delta_t = std::min(delta_t, 0.001); // make sure delta_t is not too small
    sim.time_steps = std::ceil(sim.duration/sim.delta_t); //assuming duration is defined somewhere
    // sim.delta_t = 0.1;
    // nt = ceil(duration/dt);

    // logger_->info("DFDM::Element::update_dt_nt: Time steps: {}, Delta size: {}", sim.time_steps, sim.delta_t);


}

void DFDM::Element2D::grid_refine(DFDM::simulation& sim){
    //p = order of b1, pt = order of b2, Vp = velocity, rho = density 
    // setting pt =2 for now
    // std::cout<< "mu in practice at refining:" << sim.model.mu[0] << std::endl;
    element_grid.refine(order_b1, 2, sim.model, gauss_order, ops_x, ops_z);
}