#include "element.hpp"

Mesh::Element::Element(){
    global_id = -1;
    region_id = -1;
    domain_id = -1;

    for (int i = 0; i < 4; ++i) {
        neighbors[i] = -1;
        face_connections_nbrs[i] = -1;
        face_connections_self[i] = -1;
        rotation_nbr[i] = 0;
        rotation_self[i] = 0;
        neighbor_Nx1[i] = 0;
        neighbor_Nz1[i] = 0;
        neighbor_px1[i] = 0;
        neighbor_pz1[i] = 0;
    }
}

void Mesh::Element::update_neighbor_cpus(){
    for (int i = 0; i < 4; ++i) {
        neighbor_process[i] = owner_cpu;
    }
}