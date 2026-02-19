#pragma once

#include <vector>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <string>
#include "matrix.hpp"
#include "ProcessData.hpp"
#include "element.hpp"
#include "functional"

namespace Mesh{
    void loadBalance(std::vector<Mesh::ProcessData> &cpus_data, std::vector<Mesh::ProcessData> &cpus_data_balanced);
    void findNeighbors(std::vector<Mesh::ProcessData> &cpus_data);
    void findFaces(std::vector<Mesh::ProcessData> &cpus_data);
    void checkValidityNbrs(std::vector<Mesh::ProcessData> &cpus_data);
    void printGrid(std::vector<Mesh::ProcessData> &cpus_data, std::string output_dir);
    double interpolate(std::vector<double> param, std::vector<double> r_param, int num_layers, double r, bool lower_edge, bool upper_edge);
}