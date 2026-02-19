#pragma once


#include <vector>
#include <cstdint>

#define EPS 2.2204e-16

namespace DFDM{
    class Knots{
        public:
            std::vector<double> knot_t1;
            std::vector<double> knot_t2;
            // integration intervals, initiated from knot vectors
            std::vector<double> intervals_start;
            std::vector<double> intervals_end;
            
            void init_knots(uint32_t N, uint32_t k);
        private:
            void gen_knot(uint32_t N, uint32_t k, std::vector<double>& t_vec);
    };
}