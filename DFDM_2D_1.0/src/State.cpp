
#include "State.hpp"

DFDM::State::State(){
        //
}

// mo = 0, po = 1, om = 2, op = 3
void DFDM::State::init(uint64_t Nx1, uint64_t Nz1, uint64_t nt, std::shared_ptr<spdlog::logger> logger_){
        logger = logger_;
        uint64_t Nx2 = Nx1 - 1;
        uint64_t Nz2 = Nz1 - 1;

        dU2dtt12.resize(Nx1, Nz2);
        dU2dtt21.resize(Nx2, Nz1);

        U12_0.resize(Nx1, Nz2);
        U12_1.resize(Nx1, Nz2);
        U12.resize(Nx1, Nz2);

        U21_0.resize(Nx2, Nz1);
        U21_1.resize(Nx2, Nz1);
        U21.resize(Nx2, Nz1);

        U21mo.resize(1, Nz1);
        U21po.resize(1, Nz1);
        U21mo_inn.resize(1, Nz1);
        U21po_inn.resize(1, Nz1);
        U21mo_out.resize(1, Nz1);
        U21po_out.resize(1, Nz1);

        U12om.resize(Nx1, 1);
        U12op.resize(Nx1, 1);
        U12om_inn.resize(Nx1, 1);
        U12op_inn.resize(Nx1, 1);
        U12om_out.resize(Nx1, 1);
        U12op_out.resize(Nx1, 1);
        
        U12mo.resize(1, Nz2);
        U12po.resize(1, Nz2);
        U12mo_inn.resize(1, Nz2);
        U12po_inn.resize(1, Nz2);
        U12mo_out.resize(1, Nz2);
        U12po_out.resize(1, Nz2);


        U21om.resize(Nx2, 1);
        U21op.resize(Nx2, 1);
        U21om_inn.resize(Nx2, 1);
        U21op_inn.resize(Nx2, 1);
        U21om_out.resize(Nx2, 1);
        U21op_out.resize(Nx2, 1);

        Sxx11.resize(Nx1, Nz1);
        Szz11.resize(Nx1, Nz1);

        Sxx22.resize(Nx2, Nz2);
        Szz22.resize(Nx2, Nz2);

        Sxx11mo.resize(1, Nz1);
        Sxx11po.resize(1, Nz1);
        Sxx11mo_inn.resize(1, Nz1);
        Sxx11po_inn.resize(1, Nz1);
        Sxx11mo_inn_r.resize(1, Nz1);
        Sxx11po_inn_r.resize(1, Nz1);
        Sxx11mo_out.resize(1, Nz1);
        Sxx11po_out.resize(1, Nz1);

        Sxx11om_inn.resize(Nx1, 1);
        Sxx11op_inn.resize(Nx1, 1);

        Sxx22mo.resize(1, Nz2);
        Sxx22po.resize(1, Nz2);
        Sxx22mo_inn.resize(1, Nz2);
        Sxx22po_inn.resize(1, Nz2);
        Sxx22mo_inn_r.resize(1, Nz2);
        Sxx22po_inn_r.resize(1, Nz2);
        Sxx22mo_out.resize(1, Nz2);
        Sxx22po_out.resize(1, Nz2);

        Sxx22om_inn.resize(Nx2, 1);
        Sxx22op_inn.resize(Nx2, 1);

        Szz11om.resize(Nx1, 1);
        Szz11op.resize(Nx1, 1);
        Szz11om_inn.resize(Nx1, 1);
        Szz11op_inn.resize(Nx1, 1);
        Szz11om_out.resize(Nx1, 1);
        Szz11op_out.resize(Nx1, 1);

        Szz11mo_inn.resize(1, Nz1);
        Szz11po_inn.resize(1, Nz1);
        Szz11mo_inn_r.resize(1, Nz1);
        Szz11po_inn_r.resize(1, Nz1);

        Szz22om.resize(Nx2, 1);
        Szz22op.resize(Nx2, 1);
        Szz22om_inn.resize(Nx2, 1);
        Szz22op_inn.resize(Nx2, 1);
        Szz22om_out.resize(Nx2, 1);
        Szz22op_out.resize(Nx2, 1);

        Szz22mo_inn.resize(1, Nz2);
        Szz22po_inn.resize(1, Nz2);
        Szz22mo_inn_r.resize(1, Nz2);
        Szz22po_inn_r.resize(1, Nz2);

        for(uint64_t i = 0; i < nt; i++){
                Umid.push_back(DFDM::matrix<double>(Nx1, Nz1));
        }
        logger->debug("State::init all matrices initialized");

}