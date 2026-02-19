#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <string>
#include "matrix.hpp"

namespace Mesh{
    class Element{
        public:
            int32_t global_id;
            uint32_t domain_id; // this is basically the domain id or local rank id
            uint32_t region_id; // Which block (center, left, right, top or bottom) this element belongs to
            int32_t owner_cpu;

            double inflat;
            double theta1, theta2;
            double r1,r2;
            double rx1,rx2;
            double rz1,rz2;
            double rad0, rad1;
            double coeff;
            int npart;

            uint64_t nax, naz, NAX, NAZ; // capital ones are sizes for refined matrix
            DFDM::matrix<double> xa; // (nax, naz)
            DFDM::matrix<double> za; // (nax, naz)
            int32_t neighbors[4]; // 
            int32_t neighbor_process[4];
            int32_t face_connections_nbrs[4]; // the connected faces in corresponding element (neighbors face
            int32_t face_connections_self[4]; // which of my faces are connected to which neighbor ()
            bool rotation_nbr[4]; // 0 no rotation, 1 rotation, rotation of neighbor's face relative to my face
            bool rotation_self[4]; // 0 no rotation, 1 rotation, rotation of my face relative to neighbor's face
            uint32_t neighbor_Nx1[4];
            uint32_t neighbor_Nz1[4];
            uint32_t neighbor_px1[4];
            uint32_t neighbor_pz1[4];
            uint32_t my_cpu; // current CPU that owns this element

            Element();
            void update_neighbor_cpus();


    };
    
}