#include "ElementData.hpp"

DFDM::ElementData::ElementData(){
    global_id = 0;
    local_id = 0;
    region_id = 0;
    domain_id = 0;

    for (int i = 0; i < 4; ++i) {
        neighbors[i] = -1;
        face_connections_nbrs[i] = -1;
        face_connections_self[i] = -1;
        neighbor_Nx1[i] = 0;
        neighbor_Nz1[i] = 0;
        neighbor_px1[i] = 0;
        neighbor_pz1[i] = 0;
    }
}

void DFDM::ElementData::update_neighbor_cpus(){
    for (int i = 0; i < 4; ++i) {
        neighbor_process[i] = owner_cpu;
    }
}