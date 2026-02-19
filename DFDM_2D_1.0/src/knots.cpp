#include "knots.hpp"




void DFDM::Knots::init_knots(uint32_t N, uint32_t k){
    uint32_t p = k-1;
    for(uint32_t i = 0; i < N+k; i++){
        knot_t1.push_back(0);
    }

    for(uint32_t i = 0; i < (N-1)+(k-1); i++){
        knot_t2.push_back(0);
    }

    gen_knot(N, k, knot_t1);
    gen_knot(N-1, k-1, knot_t2);

    for(uint32_t i = p; i <  N+k - p-1; i++){
        intervals_start.push_back(knot_t1[i]);
    }

    for(uint32_t i = p + 1; i <  N+k - p; i++){
        intervals_end.push_back(knot_t1[i]);
    }
}

void DFDM::Knots::gen_knot(uint32_t N, uint32_t k, std::vector<double>& t_vec){
    double dx = (double)1.0/ (N + 1 - k);

    for(uint32_t i = k-1; i < N + 1; i++){
        t_vec[i] = (double)((i+1)-k)*dx; // i+1 to compensate for matlab index mismatch
    }
     
    for(uint32_t i = 0; i < k; i++){
        t_vec[i] = (double)(t_vec[k-1] - 20.0*EPS); //1D had 10*EPS
    } 

    for(uint32_t i = N; i < (N + k ); i++){
        t_vec[i] = (double)(t_vec[N] + 20.0*EPS);
    }
}