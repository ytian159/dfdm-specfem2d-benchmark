#include "matrix.hpp"
// #include "utils.hpp"
#include <vector>
#include "grid.hpp"
#include "spdlog/spdlog.h"
namespace DFDM{
    class receiver{
        public:
            std::shared_ptr<spdlog::logger> logger;
            uint32_t relement_id; //
            // uint64_t sloc; // (elem_size+1)/2  (middle of element)
            // double src_grid; //source value on grid x
            double xglobal;
            double zglobal;
            double xlocal;
            double zlocal;
            uint64_t Nx1;
            uint64_t Nx2;
            uint64_t Nz1;
            uint64_t Nz2;
            uint32_t px1;
            uint32_t px2;
            uint32_t pz1;
            uint32_t pz2;
            uint32_t kx1;
            uint32_t kx2;
            uint32_t kz1;
            uint32_t kz2;

            DFDM::matrix<double> ur; // for benchmarking
            DFDM::matrix<double> sigr; // for benchmarking
            DFDM::matrix<double> gbx1; // sgb1
            DFDM::matrix<double> gbx2; // sgb2
            DFDM::matrix<double> gbz1; // sgb3
            DFDM::matrix<double> gbz2; // sgb4

            DFDM::matrix<double> gbx_inv1;
            DFDM::matrix<double> gbx_inv2;
            DFDM::matrix<double> gbz_inv1;
            DFDM::matrix<double> gbz_inv2; 

            std::vector<double> tx1;
            std::vector<double> tx2;
            std::vector<double> tz1;
            std::vector<double> tz2;

            receiver();
            void receiver_init(uint32_t global_id, uint32_t Nx1, uint32_t Nz1, uint32_t order_b1, uint64_t time_steps, std::shared_ptr<spdlog::logger> logger_);
            void gen_receiver(double delta_t);
            void get_receiver_location(const DFDM::Grid& element_grid);
            
    };
}