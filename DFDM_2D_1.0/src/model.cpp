#include "model.hpp"
#include <iostream>
#include "utils.hpp"
#include "spdlog/spdlog.h"



DFDM::model::model()
                            : vp({0.}), rho({2000.}), mu({0.}), num_layers(1), Model("uniform"),r({R_EARTH})
                            , logger(nullptr)
                            {}

void DFDM::model::model_init(std::string Model_, std::string model_file, std::shared_ptr<spdlog::logger> logger_){ 

    logger = logger_;

    Model = Model_;

    if (Model == "ak135-f"){

        std::ifstream infile(model_file);
        if (!infile.is_open()) {
            throw std::runtime_error("Unable to open file " + model_file);
        }else{
        }

        std::string line;
        std::getline(infile, line); // skipping first two lines
        std::getline(infile, line);
        // reading layers from third line of model file
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

            infile >> r[i] >> comma >> rho[i] >> comma >> vp[i] >> comma >> dummy1 >> comma >> dummy2 >> comma >> dummy3;
            
            vp[i] *= 1000.0; // Convert km/s to m/s
            rho[i] *= 1000.0; // Convert g/cm^3 to kg/m^3
            infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip to next line
            mu[i] = rho[i]*vp[i]*vp[i];
        }
        infile.close();
    }else if (Model == "uniform"){
        // Default uniform model already set in constructor
        logger->info("Using uniform model with vp: {}, rho: {}, mu: {}", vp[0], rho[0], mu[0]);
    }
}

double DFDM::model::get_value(std::string parameter, double x_global, double z_global){
    
    double radius = std::sqrt(x_global * x_global + z_global * z_global);
    // std::cout<< "mu in practice:" << mu[0] << std::endl;
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
        return DFDM::utils::interpolate(param, r, num_layers, radius, false, false); // no boundary handling for now
    }
    else{
        logger->error("Unknown model encountered: {}", Model);
        throw std::runtime_error("Unknown model: " + Model);
        return 0.0; // default case, should not happen
    }
}