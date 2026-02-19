#include "operators.hpp"
#include "timing.hpp"
double DFDM::Operators::get_maxdiff(std::vector<double>& vec, std::vector<double>& yo){
    double max_val = 0;

    for(int i = 0; i < vec.size(); i++){
       double temp =  std::abs(yo[i]-vec[i]); // abs and std::abs behave differently in the presence of cmath
       if(temp > max_val){
        max_val = temp;
       }
    }
    return max_val;
}


DFDM::Operators::Operators(OP_TYPE op_dim_){
    op_dim = op_dim_;
}

void DFDM::Operators::init_operators(uint64_t N_, uint64_t p_, uint32_t gauss_order, uint64_t global_id, std::shared_ptr<spdlog::logger> logger_){
    logger = logger_;
    element_id = global_id;
    ord_gi = gauss_order;
    N1 = N_;
    N2 = N1 - 1;
    p1 = p_;
    p2 = p1 - 1;
    k1 = p1 + 1;
    k2 = k1 - 1;
    logger->debug("Operators::init_operator: operators initialized");
}

// bspline looks correct now.
void DFDM::Operators::gen_bspline(){ 
    B1.resize(N1, N1);
    B2.resize(N2, N1);
    B1_t.resize(N1, N1);
    B2_t.resize(N1, N2);
    logger->debug("Operators::gen_bspline: generating bsplines");
    std::vector<double> ps1(N1);
    double delta = 1.0 / (N1 - 1);
    std::iota(ps1.begin(), ps1.end(), 0);
    std::for_each(ps1.begin(), ps1.end(), [delta](double& val) { val *= delta; });

    logger->debug("Operators::gen_bspline: ps1 generated, ps1.size:{}", ps1.size());

    #pragma omp parallel for collapse(2) schedule(auto) 
    for(uint32_t i = 0; i < N1; i++){
        for(uint32_t j = 0; j < N1; j++){
            double xj = ps1[j];
            B1.value_set(i, j, DFDM::utils::bspline_c(t1, N1, i, k1, xj));
        }
    }
    logger->debug("Operators::gen_bspline: B1 generated");
    B1.transpose(B1_t);
    logger->debug("Operators::gen_bspline: B1_t generated");
// second bspline for dim
    logger->debug("Operators::gen_bspline: generating B2");
    #pragma omp parallel for collapse(2) schedule(auto) 
    for(uint32_t i = 0; i < N2; i++){
        for(uint32_t j = 0; j < N1; j++){
            double x0 = ps1[j];
            B2.value_set(i, j, DFDM::utils::bspline_c(t2, N2, i, k2, x0));
        }
    }
    logger->debug("Operators::gen_bspline: B2 generated");
    B2.transpose(B2_t);
    logger->debug("Operators::gen_bspline: B2_t generated");

    #ifdef ENABLE_TEST
    if(element_id == 2 || element_id == 7){
        if(op_dim == DFDM::OP_TYPE::X){
            B1.print_file("./tests/generated_mats/B1_X_"+std::to_string(element_id));
            B2.print_file("./tests/generated_mats/B2_X_"+std::to_string(element_id));
        }else if(op_dim == DFDM::OP_TYPE::Z){
            B1.print_file("./tests/generated_mats/B1_Z_"+std::to_string(element_id));
            B2.print_file("./tests/generated_mats/B2_Z_"+std::to_string(element_id));
        }
    }
    #endif

}

void DFDM::Operators::gaussian_int(uint32_t N, double a, double b, std::vector<double> &xints, std::vector<double> &wints){
    N= N-1;
    int N1 = N+1, N2=N+2;
    std::vector<double> xu, y_;
    int width = 1- (-1);
    double step = (double)width/(N);
    for(double i = -1; i <= 1; i += step){
        xu.push_back(i);
    }
    //initial guess
    for(int i = 0; i < N+1; i++){
        double guess =  cos((2*i + 1)*M_PI/(2*N+2))+(0.27/N1)*sin(M_PI*xu[i]*N/N2);
        y_.push_back(guess);
    }

    std::vector<std::vector<double>> L;
    std::vector<double> Lp(N1,0);
    for(int i=0; i < N2; i++){
        L.emplace_back(std::vector<double>(N1,0));
    }

    std::vector<double> yo(N1,2);
    while(get_maxdiff(y_,yo) > EPS){
        fill(L[0].begin(), L[0].end(), 1.0);

        L[1] = y_;
        for(int i = 1; i < N1; i++){
            for(uint32_t j = 0; j< L[i+1].size(); j++){
                L[i+1][j] = ((2*(i+1)-1)*y_[j]*L[i][j] - i*L[i-1][j]) / (i+1);
            }
        }

        for(int j = 0; j < Lp.size(); j++){
            Lp[j] = N2 * ((L[N1-1][j]) - (y_[j]*L[N2-1][j]))/(1-(y_[j]*y_[j]));
        }
        yo = y_;
        for(int j = 0; j < yo.size(); j++){
            y_[j] = yo[j] - L[N2-1][j]/Lp[j];
        }

    } // while end

    for(int j=0; j < xints.size(); j++){
        xints[j] = (a*(1-y_[j]) + b*(1+y_[j]))/2;
        double sq1=((double)N2/N1);
        double sq2 = y_[j]*y_[j];
        double sq3 = Lp[j]*Lp[j];
        wints[j] = (b-a)/((1-sq2)*(sq3))*(sq1*sq1);
    }
    std::reverse(xints.begin(), xints.end());
}

void DFDM::Operators::gen_knots(){
    logger->debug("Operators::gen_knots: generating knot vectors with N1:{}, k1:{}, N2:{}, k2:{}", N1, k1, N2, k2);
    t1 = DFDM::utils::get_knot_vector(N1, k1);
    t2 = DFDM::utils::get_knot_vector(N2, k2);
    logger->debug("Operators::gen_knots: knot vectors generated with sizes: t1.size:{}, t2.size:{}", t1.size(), t2.size());
}

void DFDM::Operators::gen_gauss_params(){
    logger->debug("Operators::gen_gauss_params: generating gauss points");
    uint32_t NB_intervals = N1 - p1; // refined dim - bspline order, changing this
    DFDM::matrix<double> intervals(NB_intervals, 2);
    logger->debug("Operators::gen_gauss_params: intervals generated NB_intervals:{}", NB_intervals);
    for(uint32_t i = 0; i < NB_intervals; i++){
        intervals.value_set(i, 0, t1[p1 + i]);
        intervals.value_set(i, 1, t1[p1 + 1 + i]);
    }
    logger->debug("Operators::gen_gauss_params: intervals set");
    logger->debug("Operators::gen_gauss_params: entering gaussian_int loop, org_gi:{}", ord_gi);

    for(uint32_t i = 0; i < NB_intervals; i++){
        std::vector<double> temp_x(ord_gi, 0), temp_w(ord_gi, 0);
        gaussian_int(ord_gi, intervals(i,0), intervals(i,1), temp_x, temp_w);
        xint.push_back(temp_x);
        wxint.push_back(temp_w);
    }
    logger->debug("Operators::gen_gauss_params: gauss points generated");

}

void DFDM::Operators::stiffness_matrix(uint32_t kind, DFDM::matrix<double>& kk){
    uint32_t NB_intervals = xint.size();
    uint32_t ord_gi = xint[0].size();
    uint32_t k = p1+1;

    if(kind == 0){
        kk.resize(N1, N1-1);
        #pragma omp parallel for collapse(2) schedule(auto)// parallelizing the loop
        for(uint32_t ib1 = 0; ib1 < 2*p1; ib1++){
            for(uint32_t jb1 = 0; jb1 < 3*p1; jb1++){
                for(uint32_t kd = 0; kd < NB_intervals; kd++){
                    for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                        double b1tmp1 = DFDM::utils::dbspline_c(t1, N1, ib1, k, xint[kd][lpt], 1);
                        double b1tmp2 = DFDM::utils::bspline_c(t2, N1-1, jb1, k-1, xint[kd][lpt]);
                        kk.value_set(ib1, jb1, kk(ib1, jb1) + b1tmp1*b1tmp2*wxint[kd][lpt]);
                    }
                }
            }
        }

        for(uint32_t i = 2*p1+1 -1; i < N1-p1*2; i++){// additional -1 for matlab to c++
            for(uint32_t j = i-p1 - 1; j < i+p1; j++){
                kk.value_set(i, j, kk(i-1, j-1)); 
            }
        }

        #pragma omp parallel for collapse(2) schedule(auto)// parallelizing the loop
        for(uint32_t ib1 = N1-p1*2+1 - 1; ib1 < N1; ib1++){
            for(uint32_t jb1 = N1-3*p1 - 1; jb1 < N1-1; jb1++){
                for(uint32_t kd = 0; kd < NB_intervals; kd++){
                    for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                        double b1tmp1 = DFDM::utils::dbspline_c(t1, N1, ib1, k, xint[kd][lpt], 1);
                        double b1tmp2 = DFDM::utils::bspline_c(t2, N1-1, jb1, k-1, xint[kd][lpt]);
                        kk.value_set(ib1, jb1, kk(ib1, jb1)+ b1tmp1*b1tmp2*wxint[kd][lpt]); // modified indexing
                    }
                }
            }
        }
    } else {
        kk.resize(N1-1, N1);
        #pragma omp parallel for collapse(2) schedule(auto)// parallelizing the loop
        for(uint32_t ib1 = 0; ib1 < 2*p1; ib1++){
            for(uint32_t jb1 = 0; jb1 < 3*p1; jb1++){
                for(uint32_t kd = 0; kd < NB_intervals; kd++){
                    for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                        double b1tmp1 = DFDM::utils::dbspline_c(t2, N1-1, ib1, k-1, xint[kd][lpt], 1);
                        double b1tmp2 = DFDM::utils::bspline_c(t1, N1, jb1, k, xint[kd][lpt]);
                        kk.value_set(ib1, jb1, kk(ib1, jb1) + b1tmp1*b1tmp2*wxint[kd][lpt]);
                    }
                }
            }
        }

        for(uint32_t i = 2*p1+1 - 1; i < N1-p1*2; i++){// indexing modified
            for(uint32_t j = i-p1 - 1; j < i+p1; j++){
                kk.value_set(i, j, kk(i-1, j-1)); 
            }
        }

        #pragma omp parallel for collapse(2) schedule(auto)// parallelizing the loop
        for(uint32_t ib1 = N1-p1*2+1 - 1; ib1 < N1-1; ib1++){//modified indexing
            for(uint32_t jb1 = N1-3*p1 - 1; jb1 < N1; jb1++){
                for(uint32_t kd = 0; kd < NB_intervals; kd++){
                    for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                        double b1tmp1 = DFDM::utils::dbspline_c(t2, N1-1, ib1, k-1, xint[kd][lpt], 1);
                        double b1tmp2 = DFDM::utils::bspline_c(t1, N1, jb1, k, xint[kd][lpt]);
                        kk.value_set(ib1 , jb1, kk(ib1, jb1 ) + b1tmp1*b1tmp2*wxint[kd][lpt]); 
                    }
                }
            }
        }
    }
}

void DFDM::Operators::mass_matrix(uint32_t option, DFDM::matrix<double>& mm){
    uint32_t N1_, N2_;
    std::vector<double> tvec1, tvec2;
    if(option == 0){
        N1_ = N1;
        N2_ = N1;
        tvec1 = t1;
        tvec2 = t1;
    }else{
        N1_ = N2;
        N2_ = N2;
        tvec1 = t2;      
        tvec2 = t2;
    }

    mm.resize(N1_, N2_);
    uint32_t NB_intervals = xint.size();
    uint32_t ord_gi = xint[0].size();
    uint32_t k;
    if (option == 0) {
        k = p1 + 1;
    } else {
        k = p1;
    }
    #pragma omp parallel for collapse(2) schedule(auto) 
    for (uint32_t ib1 = 0; ib1 < 2 * p1; ib1++) {
        for (uint32_t jb1 = 0; jb1 < 3 * p1; jb1++) {
            for (uint32_t kd = 0; kd < NB_intervals; kd++) {
                for (uint32_t lpt = 0; lpt < ord_gi; lpt++) {
                    double b1tmp1 = DFDM::utils::bspline_c(tvec1, N1_, ib1, k, xint[kd][lpt]);
                    double b1tmp2 = DFDM::utils::bspline_c(tvec2, N2_, jb1, k, xint[kd][lpt]);
                    mm.value_set(ib1, jb1, mm(ib1, jb1) + b1tmp1 * b1tmp2 * wxint[kd][lpt]);
                }
            }
        }
    }

    for (uint32_t i = 2 * p1; i < N2_ - p1 * 2; i++) { // modified indexing for matlab
        for (uint32_t j = i - p1 - 1; j < i + p1; j++) { // modifying indexing for matlab
            mm.value_set(i, j , mm(i - 1 , j -1)); // modifying indexing for matlab
        }
    }
    
    #pragma omp parallel for collapse(2) schedule(auto) 
    for (uint32_t ib1 = N1_ - p1 * 2 + 1 - 1; ib1 < N1_; ib1++) {
        for (uint32_t jb1 = N2_ - 3 * p1 + 1 - 1; jb1 < N2_; jb1++) {
            for (uint32_t kd = 0; kd < NB_intervals; kd++) {
                for (uint32_t lpt = 0; lpt < ord_gi; lpt++) {
                    double b1tmp1 = DFDM::utils::bspline_c(tvec1, N1_, ib1, k, xint[kd][lpt]);
                    double b1tmp2 = DFDM::utils::bspline_c(tvec2, N2_, jb1, k, xint[kd][lpt]);
                    mm.value_set(ib1, jb1, mm(ib1, jb1) + b1tmp1 * b1tmp2 * wxint[kd][lpt]);
                }
            }
        }
    }
}


DFDM::basis_param DFDM::Operators::setup_basis(uint32_t N, uint32_t p) {
    logger->debug("Operators::setup_basis: generating basis with N:{}, p:{}", N, p);
    DFDM::basis_param basis_;
    uint32_t k_ = p + 1;
    if (p > N - 1) {
        throw std::invalid_argument("Error in setup_basis, p should be <=n-1 ...");
    }
    basis_.ps.resize(N);
    for (uint32_t i = 0; i < N; i++) {
        basis_.ps[i] = static_cast<double>(i) / (N - 1);
    }
    basis_.t =  DFDM::utils::get_knot_vector(N, k_);

    basis_.nb = N;
    basis_.pb = p;

    return basis_;
}

DFDM::matrix<double> DFDM::Operators::inner_product2(const DFDM::basis_param& basis1, const DFDM::basis_param& basis2){
    logger->debug("Operators::inner_product2: computing inner product for basis1 and basis2 with N1:{}, p1:{}, N2:{}, p2:{}", basis1.nb, basis1.pb, basis2.nb, basis2.pb);
    uint32_t N1_ = basis1.nb;
    uint32_t pN1 = basis1.pb;
    const std::vector<double>& t1_ = basis1.t;

    uint32_t M2_ = basis2.nb;
    uint32_t pM2 = basis2.pb;
    const std::vector<double>& t2_ = basis2.t;

    uint32_t kN1 = pN1 + 1;
    uint32_t kM2 = pM2 + 1;
    uint32_t ord_gi_ = std::ceil((double)(pN1 + pM2 + 1) / 2); // order of gaussian interpolation

    logger->debug("Operators::inner_product2: generating nodes, begin:{}, end:{} with t size:{}", pN1 + 1 - 1,  N1_ + kN1 - pN1 - 1, t1_.size());
    logger->debug("Operators::inner_product2: generating nodes, begin:{}, end:{} with t size:{}", pM2 + 1 - 1,  M2_ + kM2 - pM2 - 1, t2_.size());
    std::vector<double> nodes1(t1_.begin() + pN1 + 1 - 1, t1_.begin() + N1_ + kN1 - pN1);
    std::vector<double> nodes2(t2_.begin() + pM2 + 1 - 1, t2_.begin() + M2_ + kM2 - pM2); // index corrected for Matlab
    logger->debug("Operators::inner_product2: nodes generated, ord_gi:{}", ord_gi_);

    std::set<double> node_set;
    node_set.insert(nodes1.begin(), nodes1.end());
    node_set.insert(nodes2.begin(), nodes2.end());
    std::vector<double> nodes12(node_set.begin(), node_set.end());
    
    logger->debug("Operators::inner_product2: nodes12 generated and sorted, nodes12 size:{}", nodes12.size());
    uint32_t NB_intervals12 = nodes12.size() - 1;

    std::vector<std::vector<double>> int12;
    std::vector<std::vector<double>> wint12;
    for (uint32_t kd = 0; kd < NB_intervals12; kd++) {
        std::vector<double> temp_x(ord_gi_,0), temp_w(ord_gi_,0);
        gaussian_int(ord_gi_, nodes12[kd], nodes12[kd + 1], temp_x, temp_w);
        int12.push_back(temp_x);
        wint12.push_back(temp_w);
    }
    logger->debug("Operators::inner_product2: gauss points generated");
    DFDM::matrix<double> T12(N1_, M2_);
    logger->debug("Operators::inner_product2: Entering Inner product loop, N1_:{}, M2_:{}, NB_intervals12:{}, ord_gi_:{}, total iterations:{}", N1_, M2_, NB_intervals12, ord_gi_, N1_* M2_* NB_intervals12* ord_gi_);
    auto loop_time_s = start_timing();

    
    #pragma omp parallel for collapse(2) schedule(auto) // parallelizing the loop
    for (uint32_t ib1 = 0; ib1 < N1_; ib1++) {
        for (uint32_t jb1 = 0; jb1 < M2_; jb1++) {
            for (uint32_t kd = 0; kd < NB_intervals12; kd++) {
                for (uint32_t lpt = 0; lpt < ord_gi_; lpt++) {
                    double b1tmp1 = DFDM::utils::bspline_c(t1_, N1_, ib1, kN1, int12[kd][lpt]);
                    double b1tmp2 = DFDM::utils::bspline_c(t2_, M2_, jb1, kM2, int12[kd][lpt]);
                    // #pragma omp critical
                    {
                        T12.value_set(ib1, jb1, T12(ib1, jb1) + b1tmp1 * b1tmp2 * wint12[kd][lpt]);
                    }
                }
            }
        }
    }

    auto loop_time_t = stop_timing_ms(loop_time_s);
    // logger->info("Operators::inner_product2: Inner product loop completed, time taken: {} ms", loop_time_t);
    logger->debug("Operators::inner_product2: Inner product loop completed, time taken: {} ms", loop_time_t);
    logger->debug("Operators::inner_product2: inner product matrix computed");

    return T12;
}




std::vector<DFDM::matrix<double>> DFDM::Operators::gen_connected_ops(int32_t neighbor_id, int32_t face_id, uint32_t nbr_Nx1, uint32_t nbr_Nz1, uint32_t nbr_px1, uint32_t nbr_pz1){
    std::vector<DFDM::matrix<double>> conn_mats;
    auto basis_N1 = setup_basis(N1, p1); // N and p values from the current element
    auto basis_N2 = setup_basis(N2, p2);
    logger->debug("Operators::gen_connected_ops: basis generated for current element, N1:{}, p1:{}, N2:{}, p2:{}", N1, p1, N2, p2);
    // this is an example on how to do, this is for the top neighbor (X ops)
    if (neighbor_id != -1) {
        if (face_id == 0 || face_id == 1) { // get z from z; so Dzz
            auto Mz1 = nbr_Nz1; 
            auto Mz2 = Mz1 - 1;
            auto pMz1 = nbr_pz1; 
            auto pMz2 = pMz1 - 1;
            logger->debug("Operators::gen_connected_ops: generating basis for neighbor element (neighbor:{}, face:{}), Mz1:{}, pMz1:{}, Mz2:{}, pMz2:{}", neighbor_id, face_id, Mz1, pMz1, Mz2, pMz2);
            
            auto basis_Mz1 = setup_basis(Mz1, pMz1);
            auto basis_Mz2 = setup_basis(Mz2, pMz2);
            conn_mats.push_back(inner_product2(basis_N2, basis_Mz1));
            conn_mats.push_back(inner_product2(basis_N1, basis_Mz2));
            conn_mats.push_back(inner_product2(basis_N1, basis_Mz1));
            conn_mats.push_back(inner_product2(basis_N2, basis_Mz2));
            logger->debug("Operators::gen_connected_ops: inner product matrices computed for neighbor element (neighbor:{}, face:{}), Mz1:{}, pMz1:{}, Mz2:{}, pMz2:{}", neighbor_id, face_id, Mz1, pMz1, Mz2, pMz2);
        } else if(face_id == 2 || face_id == 3){ // flip = 2 or 3, get z from x; so Dzx
            auto Mx1 = nbr_Nx1; 
            auto Mx2 = Mx1 - 1;
            auto pMx1 = nbr_px1; 
            auto pMx2 = pMx1 - 1;
            logger->debug("Operators::gen_connected_ops: generating basis for neighbor element (neighbor:{}, face:{}), Mx1:{}, pMx1:{}, Mx2:{}, pMx2:{}", neighbor_id, face_id, Mx1, pMx1, Mx2, pMx2);
            auto basis_Mx1 = setup_basis(Mx1, pMx1);
            auto basis_Mx2 = setup_basis(Mx2, pMx2);
            conn_mats.push_back(inner_product2(basis_N2, basis_Mx1));
            conn_mats.push_back(inner_product2(basis_N1, basis_Mx2));
            conn_mats.push_back(inner_product2(basis_N1, basis_Mx1));
            conn_mats.push_back(inner_product2(basis_N2, basis_Mx2));
            logger->debug("Operators::gen_connected_ops: inner product matrices computed for neighbor element (neighbor:{}, face:{}), Mx1:{}, pMx1:{}, Mx2:{}, pMx2:{}", neighbor_id, face_id, Mx1, pMx1, Mx2, pMx2);
        }else{
            logger->error("Invalid face ID when computing connected transform matrices. File: {}, Line: {}", __FILE__, __LINE__);
        }
    }else{
        logger->error("No neighbor found while computing connected transform matrices. File: {}, Line: {}", __FILE__, __LINE__);
    }

    return conn_mats;
}


void DFDM::Operators::compute_operators(){
    //KK12
    uint32_t option;
    option = 0;
    logger->debug("Operators::compute_operators: computing stiffness matrix with option:{}", option);
    stiffness_matrix(option, KK12);
    //KK21
    option = 1;
    logger->debug("Operators::compute_operators: computing stiffness matrix with option:{}", option);
    stiffness_matrix(option, KK21);

    //MM11
    option = 0;
    logger->debug("Operators::compute_operators: computing mass matrix with option:{}", option);
    mass_matrix(option, MM11);
    //MM22
    option = 1;
    logger->debug("Operators::compute_operators: computing mass matrix with option:{}", option);
    mass_matrix(option, MM22);

    logger->debug("Operators::compute_operators: stiffness and mass matrices computed");

    MM11.sqrtm_schur(L11, L11_t);
    L11.inverse(invL11);
    L11_t.inverse(invL11_t);
    logger->debug("Operators::compute_operators: cholesky factorization and inverse computed for M11");
    MM22.sqrtm_schur(L22, L22_t);
    L22.inverse(invL22);
    L22_t.inverse(invL22_t);

    // #ifdef DEBUG
        // // if(element_id == 0){
        //     if(op_dim == DFDM::OP_TYPE::X){
        //         KK12.print_file("./output/KK12_X_"+std::to_string(element_id));
        //         KK21.print_file("./output/KK21_X_"+std::to_string(element_id));
        //         MM11.print_file("./output/MM11_X_"+std::to_string(element_id));
        //         MM22.print_file("./output/MM22_X_"+std::to_string(element_id));                
        //         invL11.print_file("./output/invL11_X_"+std::to_string(element_id));
        //         invL22.print_file("./output/invL22_X_"+std::to_string(element_id));
        //     }else if(op_dim == DFDM::OP_TYPE::Z){
        //         KK12.print_file("./output/KK12_Z_"+std::to_string(element_id));
        //         KK21.print_file("./output/KK21_Z_"+std::to_string(element_id));
        //         MM11.print_file("./output/MM11_Z_"+std::to_string(element_id));
        //         MM22.print_file("./output/MM22_Z_"+std::to_string(element_id));                
        //         invL11.print_file("./output/invL11_Z_"+std::to_string(element_id));
        //         invL22.print_file("./output/invL22_Z_"+std::to_string(element_id));
        //     }
        // }
    // #endif
    logger->debug("Operators::compute_operators: cholesky factorization and inverse computed for M22");



}