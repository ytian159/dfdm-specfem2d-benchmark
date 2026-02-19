#include "model.hpp"
#include <iostream>
#include "utils.hpp"



DFDM::model::model()
                            : vp({0.}), rho({2000.}), mu({0.}), num_layers(1), Model("uniform"),r({R_EARTH})
                            
                            {}

void DFDM::model::model_init(std::string Model_, std::string model_file){ 
    
    Model = Model_;

    if (Model.compare("ak135-f")==0){
        std::cout << "Model is ak135-f" << std::endl;

        std::string file_name = model_file;

        std::ifstream infile(file_name);
        if (!infile.is_open()) {
            throw std::runtime_error("Unable to open file " + file_name);
        }else{
        }

        std::string line;
        std::getline(infile, line); // skipping first two lines
        std::getline(infile, line);
        // std::getline(infile, line);
        infile >> num_layers;

        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line


        r.resize(num_layers);
        rho.resize(num_layers);
        vp.resize(num_layers);
        mu.resize(num_layers);
        char comma;
        for (int i = 0; i < num_layers; i++) {
            
            double dummy1;
            double dummy2;
            double dummy3;

            // std::getline(infile, line);
            infile >> r[i] >> comma >> rho[i] >> comma >> vp[i] >> comma >> dummy1 >> comma >> dummy2 >> comma >> dummy3;
            
            vp[i] *= 1000.0; // Convert km/s to m/s
            rho[i] *= 1000.0; // Convert g/cm^3 to kg/m^3
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line
            mu[i] = rho[i]*vp[i]*vp[i];
        }
        infile.close();
    }
}

double DFDM::model::get_value(std::string parameter, double x_global, double z_global){
    
    double radius = std::sqrt(x_global * x_global + z_global * z_global);
    if (Model.compare("uniform")==0){
        if (parameter.compare("vp")==0) {
            return vp[0];
        }
        else if (parameter.compare("rho")==0) {
            return rho[0];
        }
        else if (parameter.compare("mu")==0) {
            return mu[0];
        }
        else{
            return 0.0; // default case, should not happen
        }
    }
    else if (Model.compare("ak135-f")==0){

        std::vector<double> param;

        if (parameter.compare("vp")==0) {
            param = vp;
        }
        else if (parameter.compare("rho")==0) {
            param = rho;
        }
        else if (parameter.compare("mu")==0) {
            param = mu;
        }
        return Mesh::interpolate(param, r, num_layers, radius, false, false); // no boundary handling for now
    }
    else{
        return 0.0; // default case, should not happen
    }
}

double DFDM::model::get_max(std::string parameter, double rmin, double rmax){
    
    double out =  9999999;

    if (Model.compare("uniform")==0){
        if (parameter.compare("vp")==0) {
            return vp[0];
        }
        else if (parameter.compare("rho")==0) {
            return rho[0];
        }
        else if (parameter.compare("mu")==0) {
            return mu[0];
        }
        else{
            return 0.0; // default case, should not happen
        }
    }
    else if (Model.compare("ak135-f")==0){

        std::vector<double> param;

        if (parameter.compare("vp")==0) {
            param = vp;
        }
        else if (parameter.compare("rho")==0) {
            param = rho;
        }
        else if (parameter.compare("mu")==0) {
            param = mu;
        }

        for (int i = 0; i < num_layers; i++) {
            if (rmin < r[i] && rmax > r[i]) { // upper edge
                out = std::min(param[i],out);
            }
            else if (rmin == r[i]){ // lower edge
                if (rmin == r[i + 1]){
                    out = std::min(param[i+1],out);
                }
                else{
                    out = std::min(param[i],out);
                }
                
            }
            else if (rmax == r[i]){ // upper edge)
                if (i > 0 && rmax == r[i - 1]){
                    out = std::min(param[i-1],out);
                }
                else{
                    out = std::min(param[i],out);
                }
            }
        }
        if(out == 0){
            std::cout << "get_max:: rmin: " << rmin << ", rmax: " << rmax << ", max " << parameter << ": " << out << std::endl;
            std::cout << "Model:" << Model << ", num_layers: " << num_layers << std::endl;
            
            for(int i=0; i<num_layers; i++){
                std::cout << "Layer " << i << ", r: " << r[i] << ", " << parameter << ": " << param[i] << std::endl;
            }
            exit(-1);
        }
        return out;
    }
    else{
        std::cout << "ERROR::Uknown model type in get_max!" << std::endl;
        exit(-1);
        return 0.0; // default case, should not happen
    }
}