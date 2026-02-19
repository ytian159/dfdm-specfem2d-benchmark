#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <string>
#include "matrix.hpp"

namespace DFDM{
    class ElementData{
        public:
            uint32_t global_id;
            uint32_t local_id;
            uint32_t region_id;
            uint32_t domain_id; // this is basically the domain id or local rank id
            uint32_t owner_cpu;

            double inflat;
            double theta1, theta2;
            double r1,r2;
            double rx1,rx2;
            double rz1,rz2;
            double coeff;
            int npart;
            double rad0,rad1;

            uint32_t nax, naz, NAX, NAZ; // capital ones are sizes for refined matrix
            DFDM::matrix<double> xa; // (nax, naz)
            DFDM::matrix<double> za; // (nax, naz)
            int32_t neighbors[4]; // 
            int32_t neighbor_process[4];
            int32_t face_connections_nbrs[4]; // the connected faces in corresponding element
            int32_t face_connections_self[4]; // my face that connects the neighboring element
            bool rotation_nbr[4]; // 0 no rotation, 1 rotation
            bool rotation_self[4]; // 0 no rotation, 1 rotation
            uint32_t neighbor_Nx1[4];
            uint32_t neighbor_Nz1[4];
            uint32_t neighbor_px1[4];
            uint32_t neighbor_pz1[4];
            uint32_t my_cpu; // current CPU that owns this element

            ElementData();
            void update_neighbor_cpus();


    };
    
}