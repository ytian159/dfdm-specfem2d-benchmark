    
#include "utils.hpp"
#include <filesystem>
void Mesh::loadBalance(std::vector<Mesh::ProcessData> &cpus_data, std::vector<Mesh::ProcessData> &cpus_data_balanced){
    std::cout << "Load balancing..." << std::endl;
    std::vector<Mesh::Element> combined_elements;
    for (int j = 0; j < cpus_data.size(); j++) {
        for (int x = 0; x < cpus_data[j].element_list.size(); x++) {
            combined_elements.push_back(cpus_data[j].element_list[x]);
        }
    }
    std::cout << "Total elements in combined list: " << combined_elements.size() << std::endl;

    // sort the elements based on the size of the elements
    std::sort(combined_elements.begin(), combined_elements.end(), [](auto& a, auto& b) {
        return a.nax * a.naz > b.nax * b.naz;
    });
    std::cout << "LoadBalance: Sorted elements based on size" << std::endl;
    // distribute the elements to the cpus in a round robin fashion (is this the best?)
    for (int i = 0; i < combined_elements.size(); i++) {
        auto current_element = combined_elements[i];
        current_element.owner_cpu = i % cpus_data.size();
        current_element.update_neighbor_cpus();
        cpus_data_balanced[i % cpus_data.size()].element_list.push_back(current_element);
    }
}

void Mesh::findNeighbors(std::vector<Mesh::ProcessData> &cpus_data){
    std::cout << "Finding neighbors..." << std::endl;
    std::vector<std::reference_wrapper<Mesh::Element>> combined_elements;
    for (int j = 0; j < cpus_data.size(); j++) {
        for (int x = 0; x < cpus_data[j].element_list.size(); x++) {
            combined_elements.push_back(cpus_data[j].element_list[x]);
        }
    }
    std::cout << "Total elements in combined list: " << combined_elements.size() << std::endl;


    for (auto& element : combined_elements) {

        std::vector<double> x_corner_points;
        std::vector<double> z_corner_points;
        x_corner_points.push_back(element.get().xa(0, 0));
        x_corner_points.push_back(element.get().xa(element.get().nax - 1, 0));
        x_corner_points.push_back(element.get().xa(0, element.get().naz - 1));
        x_corner_points.push_back(element.get().xa(element.get().nax - 1, element.get().naz - 1));

        z_corner_points.push_back(element.get().za(0, 0));
        z_corner_points.push_back(element.get().za(element.get().nax - 1, 0));
        z_corner_points.push_back(element.get().za(0, element.get().naz - 1));
        z_corner_points.push_back(element.get().za(element.get().nax - 1, element.get().naz - 1));

        for (int iface = 0; iface < 4; iface++) {
            uint32_t indxA = -1, indxB = -1;
            if (element.get().neighbors[iface] == -1) {
                switch (iface) {
                    case 0:
                        indxA = 0;
                        indxB = 2;
                        break;
                    case 1:
                        indxA = 1;
                        indxB = 3;
                        break;
                    case 2:
                        indxA = 0;
                        indxB = 1;
                        break;
                    case 3:
                        indxA = 2;
                        indxB = 3;
                        break;
                    default:
                        std::cerr << "Error: Wrong face ID" << std::endl;
                        break;
                }

                double xmF = x_corner_points[indxA];
                double xpF = x_corner_points[indxB];
                double zmF = z_corner_points[indxA];
                double zpF = z_corner_points[indxB];

                double xm_mid = (xmF + xpF) / 2;
                double zm_mid = (zmF + zpF) / 2;

                for (auto& other_element : combined_elements) {
                    if (element.get().global_id != other_element.get().global_id) {
                        std::vector<double> x_corner_points_j;
                        std::vector<double> z_corner_points_j;
                        x_corner_points_j.push_back(other_element.get().xa(0, 0));
                        x_corner_points_j.push_back(other_element.get().xa(other_element.get().nax - 1, 0));
                        x_corner_points_j.push_back(other_element.get().xa(0, other_element.get().naz - 1));
                        x_corner_points_j.push_back(other_element.get().xa(other_element.get().nax - 1, other_element.get().naz - 1));

                        z_corner_points_j.push_back(other_element.get().za(0, 0));
                        z_corner_points_j.push_back(other_element.get().za(other_element.get().nax - 1, 0));
                        z_corner_points_j.push_back(other_element.get().za(0, other_element.get().naz - 1));
                        z_corner_points_j.push_back(other_element.get().za(other_element.get().nax - 1, other_element.get().naz - 1));
                        for (int jface = 0; jface < 4; jface++) {
                            uint32_t indxA_o = -1, indxB_o = -1;
                            if (other_element.get().neighbors[jface] == -1 && element.get().neighbors[iface] == -1) {
                                switch (jface) {
                                    case 0:
                                        indxA_o = 0;
                                        indxB_o = 2;
                                        break;
                                    case 1:
                                        indxA_o = 1;
                                        indxB_o = 3;
                                        break;
                                    case 2:
                                        indxA_o = 0;
                                        indxB_o = 1;
                                        break;
                                    case 3:
                                        indxA_o = 2;
                                        indxB_o = 3;
                                        break;
                                    default:
                                        std::cerr << "Error: Wrong face ID" << std::endl;
                                        break;
                                }

                                double xmFOther = x_corner_points_j[indxA_o];
                                double xpFOther = x_corner_points_j[indxB_o];
                                double zmFOther = z_corner_points_j[indxA_o];
                                double zpFOther = z_corner_points_j[indxB_o];

                                double xm_midOther = (xmFOther + xpFOther) / 2;
                                double zm_midOther = (zmFOther + zpFOther) / 2;
                                if((std::abs(xm_mid - xm_midOther) + std::abs(zm_mid - zm_midOther)) < 50){//<= std::numeric_limits<double>::epsilon()
                                    
                                    // std::cout << "curr element:"<<element.get().global_id << " neighbor id:"<<other_element.get().global_id<< " iface "<<iface<< " jface:" << jface << std::endl;
                                    // std::cout << "distance computed for neighbors:" <<(std::abs(xm_mid - xm_midOther) + std::abs(zm_mid - zm_midOther)) <<  std::endl; 
                                    // std::cout << "xm_mid:" << xm_mid << "  zm_mid:" << zm_mid << std::endl;
                                    // std::cout << "xm_midOther:" << xm_midOther << "  zm_midOther:" << zm_midOther << std::endl;

                                    other_element.get().neighbors[jface] = element.get().global_id;
                                    other_element.get().neighbor_process[jface] = element.get().owner_cpu;
                                    element.get().neighbors[iface] = other_element.get().global_id;
                                    element.get().neighbor_process[iface] = other_element.get().owner_cpu;
                                    element.get().neighbor_Nx1[iface] = other_element.get().NAX;
                                    element.get().neighbor_Nz1[iface] = other_element.get().NAZ;
                                    other_element.get().neighbor_Nx1[jface] = element.get().NAX;
                                    other_element.get().neighbor_Nz1[jface] = element.get().NAZ;
                                }
                            }  
                        }
                    }
                }
            }
        }
    }
}

void Mesh::findFaces(std::vector<Mesh::ProcessData> &cpus_data){

    std::vector<std::reference_wrapper<Mesh::Element>> combined_elements;
    for (int j = 0; j < cpus_data.size(); j++) {
        for (int x = 0; x < cpus_data[j].element_list.size(); x++) {
            combined_elements.push_back(cpus_data[j].element_list[x]);
        }
    }
    std::cout << "Total elements in combined list: " << combined_elements.size() << std::endl;


    for (auto& element : combined_elements) {
        std::vector<double> x_corner_points;
        std::vector<double> z_corner_points;
        x_corner_points.push_back(element.get().xa(0, 0));
        x_corner_points.push_back(element.get().xa(element.get().nax - 1, 0));
        x_corner_points.push_back(element.get().xa(0, element.get().naz - 1));
        x_corner_points.push_back(element.get().xa(element.get().nax - 1, element.get().naz - 1));

        z_corner_points.push_back(element.get().za(0, 0));
        z_corner_points.push_back(element.get().za(element.get().nax - 1, 0));
        z_corner_points.push_back(element.get().za(0, element.get().naz - 1));
        z_corner_points.push_back(element.get().za(element.get().nax - 1, element.get().naz - 1));
        int32_t nbr_element = -1;
        for (int iface = 0; iface < 4; iface++) {
            uint32_t indxA = -1, indxB = -1;
            double max_dist = std::numeric_limits<double>::max();
            switch (iface) {
                case 0:
                    nbr_element = element.get().neighbors[iface];
                    indxA = 0;
                    indxB = 2;
                    break;
                case 1:
                    nbr_element = element.get().neighbors[iface];
                    indxA = 1;
                    indxB = 3;
                    break;
                case 2:
                    nbr_element = element.get().neighbors[iface];
                    indxA = 0;
                    indxB = 1;
                    break;
                case 3:
                    nbr_element = element.get().neighbors[iface];
                    indxA = 2;
                    indxB = 3;
                    break;
                default:
                    std::cerr << "Error: Wrong face ID" << std::endl;
                    break;
            }

            if(nbr_element != -1){
                // finding the neighbor element in the combined list for extracting corners
                auto nbr_el = std::find_if(combined_elements.begin(), combined_elements.end(),
                                        [nbr_element](const auto& element) { return element.get().global_id == nbr_element; });
                if (nbr_el != combined_elements.end()) {
                    //corner points for the nbr_el
                    std::vector<double> x_corner_points_nbr;
                    std::vector<double> z_corner_points_nbr;
                    x_corner_points_nbr.push_back(nbr_el->get().xa(0, 0));
                    x_corner_points_nbr.push_back(nbr_el->get().xa(nbr_el->get().nax - 1, 0));
                    x_corner_points_nbr.push_back(nbr_el->get().xa(0, nbr_el->get().naz - 1));
                    x_corner_points_nbr.push_back(nbr_el->get().xa(nbr_el->get().nax - 1, nbr_el->get().naz - 1));

                    z_corner_points_nbr.push_back(nbr_el->get().za(0, 0));
                    z_corner_points_nbr.push_back(nbr_el->get().za(nbr_el->get().nax - 1, 0));
                    z_corner_points_nbr.push_back(nbr_el->get().za(0, nbr_el->get().naz - 1));
                    z_corner_points_nbr.push_back(nbr_el->get().za(nbr_el->get().nax - 1, nbr_el->get().naz - 1));
                    // iterating over faces of neighbor element
                    for (int jFace = 0; jFace < 4; ++jFace) {
                        int indx_nbr1, indx_nbr2;
                        switch (jFace) {
                            case 0:
                                indx_nbr1 = 0;
                                indx_nbr2 = 2;
                                break;
                            case 1:
                                indx_nbr1 = 1;
                                indx_nbr2 = 3;
                                break;
                            case 2:
                                indx_nbr1 = 0;
                                indx_nbr2 = 1;
                                break;
                            case 3:
                                indx_nbr1 = 2;
                                indx_nbr2 = 3;
                                break;
                        }

                        //corner points for the element we are comparing for
                        double xmF = x_corner_points[indxA];
                        double xpF = x_corner_points[indxB];
                        double zmF = z_corner_points[indxA];
                        double zpF = z_corner_points[indxB];
                        double xmFNbr = x_corner_points_nbr[indx_nbr1];
                        double xpFNbr = x_corner_points_nbr[indx_nbr2];
                        double zmFNbr = z_corner_points_nbr[indx_nbr1];
                        double zpFNbr = z_corner_points_nbr[indx_nbr2];

                        //estimate distance for non rotated face
                        double dist1 = std::pow(xmF-xmFNbr, 2) + std::pow(xpF-xpFNbr, 2) + std::pow(zmF-zmFNbr, 2) + std::pow(zpF-zpFNbr, 2);
                        if(dist1 < max_dist){
                            max_dist = dist1;
                            element.get().face_connections_nbrs[iface] = jFace;
                            element.get().rotation_nbr[iface] = 0;
                            element.get().face_connections_self[iface] = iface;
                            element.get().rotation_self[iface] = 0;
                        }
                        //check for rotation
                        double dist1_rot = std::pow(xmF-xpFNbr, 2) + std::pow(xpF-xmFNbr, 2) + std::pow(zmF-zpFNbr, 2) + std::pow(zpF-zmFNbr, 2);
                        if(dist1_rot < max_dist){
                            max_dist = dist1_rot;
                            element.get().face_connections_nbrs[iface] = jFace;
                            element.get().rotation_nbr[iface] = 1; // indicates rotation
                            element.get().face_connections_self[iface] = iface;
                            element.get().rotation_self[iface] = 1;
                        }
                    }
                }else {
                    std::cerr << "Error: Element with global ID " << nbr_element << " not found in findFace functions." << std::endl;
                }
            }
        }
    }
}

void Mesh::checkValidityNbrs(std::vector<Mesh::ProcessData> &cpus_data){
    std::cout << "Checking validity of neighbors..." << std::endl;
    std::vector<std::reference_wrapper<Mesh::Element>> combined_elements;
    for (int j = 0; j < cpus_data.size(); j++) {
        for (int x = 0; x < cpus_data[j].element_list.size(); x++) {
            combined_elements.push_back(cpus_data[j].element_list[x]);
        }
    }
    std::cout << "Total elements in combined list: " << combined_elements.size() << std::endl;

    for (auto& element : combined_elements) {
        for(int iface = 0; iface < 4; ++iface){
            if(element.get().neighbors[iface] != -1){
                if(element.get().face_connections_nbrs[iface] == -1){
                    std::cerr << "Error: Face connection not found for element with global ID " << element.get().global_id << " and face " << iface << std::endl;
                    exit(1);
                }
            }
        }
    }

    for(auto& element : combined_elements) {
        for(int iface = 0; iface < 4; ++iface){
            // Check if any two neighbors of the element are the same
            for (int iface1 = 0; iface1 < 4; ++iface1) {
                for (int iface2 = iface1 + 1; iface2 < 4; ++iface2) {
                    if (element.get().neighbors[iface1] != -1 && 
                        element.get().neighbors[iface2] != -1 && 
                        element.get().neighbors[iface1] == element.get().neighbors[iface2]) {
                            if (element.get().face_connections_nbrs[iface1] != -1 && 
                                element.get().face_connections_nbrs[iface2] != -1 && 
                                element.get().face_connections_nbrs[iface1] == element.get().face_connections_nbrs[iface2]) {
                                std::cerr << "Error: Duplicate neighbors with identical face connections found for element with global ID " 
                                          << element.get().global_id << " on faces " << iface1 << " and " << iface2 << std::endl;
                                exit(1);
                            }
                    }
                }
            }
        }
    }

}

void Mesh::printGrid(std::vector<Mesh::ProcessData> &cpus_data, std::string output_dir){

    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }
    std::cout << "Printing grid to directory:..." <<output_dir << std::endl;
    std::vector<std::reference_wrapper<Mesh::Element>> combined_elements;
    for (int j = 0; j < cpus_data.size(); j++) {
        for (int x = 0; x < cpus_data[j].element_list.size(); x++) {
            combined_elements.push_back(cpus_data[j].element_list[x]);
        }
    }
    std::cout << "Total elements in combined list: " << combined_elements.size() << std::endl;

    for (auto& element : combined_elements) {
        element.get().xa.print_file(output_dir+"/grid_x_"+std::to_string(element.get().global_id));
        element.get().za.print_file(output_dir+"/grid_z_"+std::to_string(element.get().global_id));
    }
}

double Mesh::interpolate(std::vector<double> param, std::vector<double> r_param, int num_layers, double r, bool lower_edge, bool upper_edge){
    double result = 0.0;
    
    for (int i = 0; i < num_layers; i++) {
        if (i==num_layers-1){
            result = param[i];
            break;
        }
        if (r > r_param[i] && r < r_param[i + 1]) {
            result = param[i] + ((r - r_param[i]) / (r_param[i + 1] - r_param[i])) * (param[i + 1] - param[i]);
            break;
        }
        if ((i<num_layers - 1) && r == r_param[i + 1] && r == r_param[i]) {
            if (lower_edge) {
                result = param[i];
                break;
            }
            else if (upper_edge) {
                result = param[i + 1];
                break;
            }
            else {
                result = (param[i] + param[i + 1])/2;
                break;
            }
        }
    }
            
    return result;
}