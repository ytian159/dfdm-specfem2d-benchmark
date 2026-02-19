#pragma once

#include "matrix.hpp"
#include "spdlog/spdlog.h"
namespace DFDM{
    class State{
    public:
        std::shared_ptr<spdlog::logger> logger;
      // acceleration
        DFDM::matrix<double> dU2dtt12;
        DFDM::matrix<double> dU2dtt21;
      // displacement
        DFDM::matrix<double> U12_0;
        DFDM::matrix<double> U12_1;
        DFDM::matrix<double> U12;

        DFDM::matrix<double> U21_0;
        DFDM::matrix<double> U21_1;
        DFDM::matrix<double> U21;

        // Displacement boundary values left + right
        // U21
        DFDM::matrix<double> U21mo;     // zeros(1,Nz1);
        DFDM::matrix<double> U21po;     // zeros(1,Nz1);
        DFDM::matrix<double> U21mo_inn; // zeros(1,Nz1);
        DFDM::matrix<double> U21po_inn; // zeros(1,Nz1);
        DFDM::matrix<double> U21mo_out; // zeros(1,Nz1);
        DFDM::matrix<double> U21po_out; // zeros(1,Nz1);

        // U12
        DFDM::matrix<double> U12mo;     // zeros(1,Nz2);
        DFDM::matrix<double> U12po;     // zeros(1,Nz2);
        DFDM::matrix<double> U12mo_inn; // zeros(1,Nz2);
        DFDM::matrix<double> U12po_inn; // zeros(1,Nz2);
        DFDM::matrix<double> U12mo_out; // zeros(1,Nz2);
        DFDM::matrix<double> U12po_out; // zeros(1,Nz2);
        
        // Displacement boundary values bottom + upper
        // U12
        DFDM::matrix<double> U12om;     // zeros(Nx1,1);
        DFDM::matrix<double> U12op;     // zeros(Nx1,1);
        DFDM::matrix<double> U12om_inn; // zeros(Nx1,1);
        DFDM::matrix<double> U12op_inn; // zeros(Nx1,1);
        DFDM::matrix<double> U12om_out; // zeros(Nx1,1);
        DFDM::matrix<double> U12op_out; // zeros(Nx1,1);
        // U21
        DFDM::matrix<double> U21om;     // zeros(Nx2,1);
        DFDM::matrix<double> U21op;     // zeros(Nx2,1);
        DFDM::matrix<double> U21om_inn; // zeros(Nx2,1);
        DFDM::matrix<double> U21op_inn; // zeros(Nx2,1);
        DFDM::matrix<double> U21om_out; // zeros(Nx2,1);
        DFDM::matrix<double> U21op_out; // zeros(Nx2,1);


        // Stress 
        // S11
        DFDM::matrix<double> Sxx11;    // zeros(Nx1,Nz1); 
        DFDM::matrix<double> Szz11;    // zeros(Nx1,Nz1);
        // 22   
        DFDM::matrix<double> Sxx22;    // zeros(Nx2,Nz2);
        DFDM::matrix<double> Szz22;    // zeros(Nx2,Nz2);
        
        // S boundary values 
        // left + right
        DFDM::matrix<double> Sxx11mo;     // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11po;     // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11mo_inn; // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11po_inn; // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11mo_inn_r; // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11po_inn_r; // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11mo_out; // zeros(1,Nz1);
        DFDM::matrix<double> Sxx11po_out; // zeros(1,Nz1);

        DFDM::matrix<double> Sxx11om_inn; // zeros(Nx1,1);
        DFDM::matrix<double> Sxx11op_inn; // zeros(Nx1,1);
        DFDM::matrix<double> Sxx11om_inn_r; // zeros(Nx1,1);
        DFDM::matrix<double> Sxx11op_inn_r; // zeros(Nx1,1);

        // left + right
        DFDM::matrix<double> Sxx22mo;     // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22po;     // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22mo_inn; // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22po_inn; // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22mo_inn_r; // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22po_inn_r; // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22mo_out; // zeros(1,Nz2);
        DFDM::matrix<double> Sxx22po_out; // zeros(1,Nz2); 

        DFDM::matrix<double> Sxx22om_inn; // zeros(Nx2,1);
        DFDM::matrix<double> Sxx22op_inn; // zeros(Nx2,1);
        DFDM::matrix<double> Sxx22om_inn_r; // zeros(Nx2,1);
        DFDM::matrix<double> Sxx22op_inn_r; // zeros(Nx2,1);
        
        // S boundary values 
        // bottom + upper
        DFDM::matrix<double> Szz11om;     // zeros(Nx1,1);
        DFDM::matrix<double> Szz11op;     // zeros(Nx1,1);
        DFDM::matrix<double> Szz11om_inn; // zeros(Nx1,1);
        DFDM::matrix<double> Szz11op_inn; // zeros(Nx1,1);
        DFDM::matrix<double> Szz11om_out; // zeros(Nx1,1);
        DFDM::matrix<double> Szz11op_out; // zeros(Nx1,1);
        DFDM::matrix<double> Szz11om_inn_r; // zeros(Nx1,1);
        DFDM::matrix<double> Szz11op_inn_r; // zeros(Nx1,1);

        DFDM::matrix<double> Szz11mo_inn; // zeros(1,Nz1);
        DFDM::matrix<double> Szz11po_inn; // zeros(1,Nz1);
        DFDM::matrix<double> Szz11mo_inn_r; // zeros(1,Nz1);
        DFDM::matrix<double> Szz11po_inn_r; // zeros(1,Nz1);

        // lower + upper
        DFDM::matrix<double> Szz22om;     // zeros(Nx2,1);
        DFDM::matrix<double> Szz22op;     // zeros(Nx2,1);
        DFDM::matrix<double> Szz22om_inn; // zeros(Nx2,1);
        DFDM::matrix<double> Szz22op_inn; // zeros(Nx2,1);
        DFDM::matrix<double> Szz22om_out; // zeros(Nx2,1);
        DFDM::matrix<double> Szz22op_out; // zeros(Nx2,1);
        DFDM::matrix<double> Szz22om_inn_r; // zeros(Nx2,1);
        DFDM::matrix<double> Szz22op_inn_r; // zeros(Nx2,1);

        DFDM::matrix<double> Szz22mo_inn; // zeros(1,Nz2);
        DFDM::matrix<double> Szz22po_inn; // zeros(1,Nz2);
        DFDM::matrix<double> Szz22mo_inn_r; // zeros(1,Nz2);
        DFDM::matrix<double> Szz22po_inn_r; // zeros(1,Nz2);

        std::vector< DFDM::matrix<double>> Umid; //stores umid at every time step in a vector. zeros(Nx1,Nz1,floor(nt));
        
        State();
        void init(uint64_t Nx1, uint64_t Nz1, uint64_t nt, std::shared_ptr<spdlog::logger> logger_);
    };
}