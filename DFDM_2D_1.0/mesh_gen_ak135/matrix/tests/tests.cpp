#include <iostream>
#include "matrix.hpp"

int main(){
    DFDM::matrix<double> m1(3,3);
    m1.value_set(0,0,4);
    m1.value_set(0,1,1);
    m1.value_set(0,2,2);
    m1.value_set(1,0,1);
    m1.value_set(1,1,9);
    m1.value_set(1,2,3);
    m1.value_set(2,0,2);
    m1.value_set(2,1,3);
    m1.value_set(2,2,16);
    // m1.value_set(0,0,1);
    // m1.value_set(0,1,2);
    // m1.value_set(0,2,3);
    // m1.value_set(1,0,4);
    // m1.value_set(1,1,5);
    // m1.value_set(1,2,6);
    // m1.value_set(2,0,7);
    // m1.value_set(2,1,8);
    // m1.value_set(2,2,9);
    DFDM::matrix<double> m1u(3,3);
    DFDM::matrix<double> m1ut(3,3);
    m1.sqrtm_schur(m1u,m1ut);

    // auto result_m = m1+m2;
    m1.print();
    std::cout << "m1u:"<< std::endl; 
    m1u.print();
    // std::cout << "miut:"<< std::endl; 
    // m1ut.print();

    auto prod = m1u.matprod(m1u);
    std::cout <<"prod:"<< std::endl; 
    prod.print();
    return 0;
}