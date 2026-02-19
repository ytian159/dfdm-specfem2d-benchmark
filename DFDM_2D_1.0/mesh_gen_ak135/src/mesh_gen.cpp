#include <iostream>
#include <cmath>
#include "element.hpp"
#include "ProcessData.hpp"
#include "toml.hpp"
#include <limits>
#include "utils.hpp"
#include "model.hpp"

int main (int argc, char **argv){
    auto config_file = argv[1];
    std::string output_dir(argv[2]);
    // std::string model_file(argv[3]);
    if(argc != 3){
        std::cout << "Usage: ./mesh_gen config_file output_dir" << std::endl;
        return -1;
    }

    double nmin = 18.0;
    double nunrefmin = 18.0;
    // double freq = 1;
    double R_EARTH = 6371.0*1000; // in meters



    auto data = toml::parse(config_file);
    // std::vector<double> rad{600*1000, 1221.5*1000, 1700*1000};
    uint32_t NCPUs = toml::find<int>(data, "NCPUs"); // multiple of 4 for keeping load balance as parts is going to be in power of 2.

    uint32_t npart = toml::find<int>(data, "nparts"); // needs to be multiple of 2, number of elements in each region
    double dxz = toml::find<double>(data, "dxz"); // distance between two points in an element. units: meters
    double inflat = toml::find<double>(data, "inflat"); // determines the curvature of the middle element's boundary 1 means a circle and 0 is flat
    std::string Model  = toml::find<std::string>(data, "Model");
    double freq = toml::find<double>(data, "frequency");
    double ppw = toml::find<double>(data, "ppw"); // points per wavelength?
    std::string model_file = toml::find<std::string>(data, "model_file");

    std::vector<double> rad = toml::find<std::vector<double>>(data, "rads");

    
    DFDM::model model;

    if (Model.compare("uniform")==0){
        std::cout << "Model is uniform" << std::endl;
        double vp  = toml::find<double>(data, "Vp");
        double vs = toml::find<double>(data, "Vs");
        double rho = toml::find<double>(data, "rho");
        
        model.Model = "uniform";
        model.vp = {vp};
        model.rho = {rho};
        model.mu = {rho*vp*vp};
        model.num_layers = 1;
        model.r = {*std::max_element(rad.begin(), rad.end())};

    }
    else if (Model.compare("ak135-f")==0){
        model.model_init("ak135-f", model_file);
    }


    if (NCPUs>( 4*(rad.size()-1) * npart + npart*npart)){
        std::cout << "NCPUs is too large, please reduce it to " << 4*(rad.size()-1) * npart + npart*npart << std::endl;
        return 0;
    }

    double fmax = freq*3;
    // double ppw = 3;
    uint32_t iparts = std::ceil((float)(npart * 4 )/NCPUs); // elements per cpu (excluding elements from cube), for total elements outside the cube multiply by 2 again, since 4 regions are symmetrical, each has a symmetric copy. Here we just compute for one and use them twice.
    std::vector<Mesh::ProcessData> cpus_data, balanced_cpus_data;
    for(uint32_t x = 0; x < NCPUs; x++){
        cpus_data.push_back(Mesh::ProcessData(x));
        balanced_cpus_data.push_back(Mesh::ProcessData(x));
    }

    double dtheta = 2*M_PI/(npart*4); // because 4 regions around the circle
    int32_t length_rad = rad.size();

    std::vector<double> allrad(2*length_rad+npart - 1, 0); // zero vector of size 2*length_rad+npart - 1;

    for(uint32_t i = 0; i < length_rad; i++){
        allrad[i] = -1* rad[(rad.size()-1) - i]; // allrad[0]=-max(rad), etc decreasing rad from 0 to length(rad)
        allrad[length_rad + npart + i - 1] = rad[i]; // allrad[-1]=max(rad), increasing rad. Hold in the middle: what is going on between length_rad and length_rad+npart?
    }

    int int_len = npart - 1;
    int fwd_dir_len = (int_len - 1) / 2;
    int startrad = rad.size() + fwd_dir_len + 1;

    for (int i = 0; i < fwd_dir_len; i++) { // closes hole in the middle
        double step = rad[0] * (1 << i) / (int_len - 1); // Using bitwise shift for power of 2
        allrad[startrad + i] = step; // center part of allrad: steps in powers of 2, in positive and negative direction
        allrad[startrad - i - 2] = -step;
    }
  
    std::cout<<"rad0:" << rad[0] << std::endl;
    std::cout << "allrad array: ";
    for (uint32_t i = 0; i < allrad.size(); i++) {
        std::cout << allrad[i] << " ";
    }
    std::cout << std::endl;
   
    auto length_allrad = allrad.size();
    int32_t ndomain = std::pow((length_allrad -1 ),2) - 4*std::pow((length_rad - 1),2); // this is the actual square domain, minus the 4 quarters.

    std::vector<Mesh::Element> full_domain;
    uint32_t region = 0;
    int32_t idomain = -1;

    double r1, r2, rx1, rx2, rz1, rz2, aR, dR; // aR = avg radius, dR is thickness
    int32_t nx, nz;
    double theta, theta1, theta2;
    double r, rabs;
    double x_tmp, z_tmp, zs, xs;
    double coeff, alpha;
    int icpu = -1;
    //region naming can be understood by following the overlapping squares in remarkable
    for(uint32_t jb = 0; jb < length_allrad - 1; jb++){ // elements in z direction
        for(uint32_t ib = 0; ib < length_allrad - 1; ib++){ // elements in x direction
            Mesh::Element current_element; 

            if(jb < length_rad - 2 && ib >= length_rad -1 && ib < length_rad + npart - 1)
                region = 1;
            else if(jb == length_rad - 2 && ib >= length_rad -1 && ib < length_rad + npart - 1)
                region = 2;
            else if(ib < length_rad - 2 && jb >= length_rad -1 && jb < length_rad + npart - 1)
                region = 3;
            else if(ib == length_rad - 2 && jb >= length_rad -1 && jb < length_rad + npart - 1)
                region = 4;
            else if(ib >= length_rad - 1 && ib < length_rad + npart -1 && jb >= length_rad - 1 && jb < length_rad + npart - 1)
                region = 5;  
            else if(ib == length_rad + npart - 1 && jb >= length_rad -1 && jb < length_rad + npart - 1)
                region = 6;   
            else if(ib > length_rad + npart - 1 && jb >= length_rad -1 && jb < length_rad + npart - 1)
                region = 7; 
            else if(jb == length_rad + npart - 1 && ib >= length_rad -1 && ib < length_rad + npart - 1)
                region = 8;   
            else if(jb > length_rad +npart - 1 && ib >= length_rad -1 && ib < length_rad + npart - 1)
                region = 9; 
            else 
                region = 0;

            switch (region){ // if nparts ==1 then element=region, if nparts=2 all outer regions in 2, central in 4 etc
                case 1:
                case 2:
                case 8:
                case 9:
                {
                    icpu = -1;
                    idomain++;
                    if(region == 1 || region == 2){
                        int iib = ((ib+1) - (length_rad) + 1);
                        icpu = std::ceil((float)iib/iparts) - 1;
                        if(icpu < 0 || icpu >= NCPUs){
                            std::cout << "ERROR (region 1, 2)! CPU can't be less than zero: "<< icpu << std::endl;
                            return -1;
                        }
                    }else if(region == 8 || region == 9){
                        int iib = ((ib+1) - (length_rad) + 1);
                        icpu = (std::ceil((float)iib/iparts) + (((float)NCPUs/4)*3)) - 1;
                        if(icpu < 0 || icpu >= NCPUs){
                            std::cout << "ERROR (region 8, 9)! CPU can't be less than zero"<< icpu << std::endl;
                            return -1;
                        }
                    }

                    current_element.global_id = idomain;
                    current_element.region_id = region;
                    current_element.domain_id = region;
                    r1 = allrad[jb];
                    r2 = allrad[jb + 1];
                    aR = std::abs(r2+r1)/2;
                    dR = std::abs(r2-r1);

                    if (region == 1 || region == 2){
                        theta1 = 5*(M_PI/4) + ((double)ib - (double)length_rad + 1)/npart * M_PI/2;
                        theta2 = 5*(M_PI/4) + ((double)ib - (double)length_rad + 2)/npart * M_PI/2;
                    }if(region == 8 || region == 9){
                        theta1 = 3*(M_PI/4) - ((double)ib - (double)length_rad + 1)/npart * M_PI/2;
                        theta2 = 3*(M_PI/4) - ((double)ib - (double)length_rad + 2)/npart * M_PI/2;
                    }

                    current_element.inflat=inflat;
                    current_element.r1 = r1;
                    current_element.r2 = r2;
                    current_element.theta1 = theta1;
                    current_element.theta2 = theta2;
                    current_element.npart = npart;
                    if (region == 1 || region == 2){
                        current_element.coeff = 1.0 - 2.0*(((double)ib - (double)length_rad + 1.0))/npart;
                    }if(region == 8 || region == 9){
                        current_element.coeff =-1.0 + 2.0*(((double)ib - (double)length_rad + 1.0))/npart;
                    }
                    current_element.rx1 = 0.0;
                    current_element.rx2 = 0.0;
                    current_element.rz1 = 0.0;
                    current_element.rz2 = 0.0;
                    current_element.rad0 = rad[0];
                    current_element.rad1 = rad[1];

                    nx = std::max(nunrefmin, std::ceil((aR*M_PI/2)/npart/dxz));
                    nz = std::max(nunrefmin, std::ceil(dR/dxz));

                    current_element.nax = nx;
                    current_element.naz = nz;

                    current_element.xa.resize(nx,nz);
                    current_element.za.resize(nx,nz);

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){

                            theta = (double)ii/(double)(nx - 1) * (theta2 - theta1) + theta1;
                            r = (double)jj/(double)(nz - 1) * (r2-r1) + r1;
                            rabs = std::abs(r);

                            x_tmp = rabs*std::cos(theta);
                            z_tmp = rabs*std::sin(theta);


                            zs = r/std::sqrt(2);
                            if (region == 1 || region == 2){
                                coeff = 1.0 - 2.0*(((double)ib - (double)length_rad + 1.0) + (double)ii/(nx - 1))/npart;
                            }if(region == 8 || region == 9){
                                coeff =-1.0 + 2.0*(((double)ib - (double)length_rad + 1.0) + (double)ii/(nx - 1))/npart;
                            }

                            xs = coeff*zs;

                            if(region == 2 || region == 8){
                                alpha = (1.0 - inflat) * (1.0 - (rabs - rad[0])/(rad[1] - rad[0]));
                                double value_x = alpha * xs + (1.0 - alpha) * x_tmp;
                                current_element.xa.value_set(ii,jj, value_x);
                                double value_z = alpha * zs + (1.0 - alpha) * z_tmp;
                                current_element.za.value_set(ii,jj, value_z);
                            }else{
                                current_element.xa.value_set(ii,jj, x_tmp);
                                current_element.za.value_set(ii,jj, z_tmp);   
                            }

                        }
                    }

                    
                    DFDM::matrix<double> rs(nx, nz); // initialize data vector
                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double x_squared = current_element.xa(ii,jj) * current_element.xa(ii,jj);
                            double z_squared = current_element.za(ii,jj) * current_element.za(ii,jj);
                            rs.value_set(ii,jj, std::sqrt(x_squared + z_squared));
                        }
                    }

                    double max_element = *std::max_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmax = max_element * (1 - 10 * std::numeric_limits<double>::epsilon());

                    double min_element = *std::min_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmin = min_element * (1 + 10 * std::numeric_limits<double>::epsilon());

                    double nnx = nmin;
                    double nnz = nmin;
                    
                    double vpmin=model.get_max("vp", rmin, rmax);
                    double vmin = vpmin;
                    double wavemin = vmin/fmax;

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double r0 = std::sqrt(std::pow(current_element.xa(ii,jj), 2) + std::pow(current_element.za(ii,jj), 2));
                            double r = std::min(std::max(rmin, r0), rmax); 

                            nnx = std::max(nnx, std::ceil((dtheta*r)/wavemin*ppw));
                            nnz = std::max(nnz, std::ceil((rmax-rmin)/wavemin*ppw));      

                        }
                    }

                    current_element.NAX = nnx;
                    current_element.NAZ = nnz;
                    current_element.owner_cpu = icpu;
                    current_element.update_neighbor_cpus();
                    cpus_data[icpu].element_list.push_back(current_element);
                    auto el_ = cpus_data[icpu].element_list[cpus_data[icpu].element_list.size() - 1];
                    break;
                }   
                case 3:
                case 4:
                case 6:
                case 7:
                {
                    idomain++;
                    icpu = -1;
                    if(region == 3 || region == 4){
                        int jjb = ((jb+1) - (length_rad) + 1);
                        icpu = (std::ceil((float)jjb/iparts) + (float)NCPUs/4) - 1;
                        
                        if(icpu < 0 || icpu >= NCPUs){
                            std::cout << "ERROR (region 3, 4)! CPU can't be less than zero"<< icpu << std::endl;
                            return -1;
                        }
                    }else if(region == 6 || region == 7){
                        int jjb = ((jb+1) - (length_rad) + 1);
                        icpu = (std::ceil((float)jjb/iparts) + (((float)NCPUs/4)*2)) - 1;
                        if(icpu < 0 || icpu >= NCPUs){
                            std::cout << "ERROR(region 6, 7)! CPU can't be less than zero"<< icpu << std::endl;
                            return -1;
                        }
                    }
                    current_element.global_id = idomain;
                    current_element.region_id = region;
                    current_element.domain_id = region;
                    r1 = allrad[ib];
                    r2 = allrad[ib + 1];
                    aR = std::abs(r2+r1)/2;
                    dR = std::abs(r2-r1);

                    if (region == 3 || region == 4){
                        theta1 = 5*(M_PI/4) - ((double)jb - (double)length_rad + 1)/npart * M_PI/2;
                        theta2 = 5*(M_PI/4) - ((double)jb - (double)length_rad + 2)/npart * M_PI/2;
                    }if(region == 6 || region == 7){
                        theta1 = 7*(M_PI/4) + ((double)jb - (double)length_rad + 1)/npart * M_PI/2;
                        theta2 = 7*(M_PI/4) + ((double)jb - (double)length_rad + 2)/npart * M_PI/2;
                    }

                    current_element.inflat=inflat;
                    current_element.r1 = r1;
                    current_element.r2 = r2;
                    current_element.theta1 = theta1;
                    current_element.theta2 = theta2;
                    current_element.npart = npart;
                    if (region == 6 || region == 7){
                        current_element.coeff = - 1 + 2*(((double)jb - (double)length_rad + 1))/npart;
                    }if(region == 3 || region == 4){
                        current_element.coeff = 1 - 2*(((double)jb - (double)length_rad + 1))/npart;
                    }
                    current_element.rx1 = 0.0;
                    current_element.rx2 = 0.0;
                    current_element.rz1 = 0.0;
                    current_element.rz2 = 0.0;
                    current_element.rad0 = rad[0];
                    current_element.rad1 = rad[1];

                    nx = std::max(nunrefmin, std::ceil(dR/dxz));
                    nz = std::max(nunrefmin, std::ceil((aR*M_PI/2)/npart/dxz));
                    current_element.nax = nx;
                    current_element.naz = nz;

                    current_element.xa.resize(nx,nz);
                    current_element.za.resize(nx,nz);

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){

                            theta = (double)jj/(double)(nz - 1) * (theta2 - theta1) + theta1;
                            r = (double)ii/((double)nx - 1) * (r2-r1) + r1;
                            rabs = std::abs(r);

                            x_tmp = rabs*std::cos(theta);
                            z_tmp = rabs*std::sin(theta);

                            xs = r/std::sqrt(2);
                            if (region == 6 || region == 7){
                                coeff = - 1 + 2*(((double)jb - (double)length_rad + 1) + (double)jj/(nz - 1))/npart;
                            }if(region == 3 || region == 4){
                                coeff = 1 - 2*(((double)jb - (double)length_rad + 1) + (double)jj/(nz - 1))/npart;
                            }

                            zs = coeff*xs;

                            if(region == 4 || region == 6){
                                alpha = (1 - inflat) * (1 - (rabs - rad[0])/(rad[1] - rad[0]));
                                double value_x = alpha * xs + (1 - alpha) * x_tmp;
                                current_element.xa.value_set(ii,jj, value_x);
                                double value_z = alpha * zs + (1 - alpha) * z_tmp;
                                current_element.za.value_set(ii,jj, value_z);
                            }else{
                                current_element.xa.value_set(ii,jj, x_tmp);
                                current_element.za.value_set(ii,jj, z_tmp);   
                            }

                        }
                    }

                    DFDM::matrix<double> rs(nx, nz); // initialize data vector
                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double x_squared = current_element.xa(ii,jj) * current_element.xa(ii,jj);
                            double z_squared = current_element.za(ii,jj) * current_element.za(ii,jj);
                            rs.value_set(ii,jj, std::sqrt(x_squared + z_squared));
                        }
                    }

                    double max_element = *std::max_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmax = max_element * (1 - 10 * std::numeric_limits<double>::epsilon());

                    double min_element = *std::min_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmin = min_element * (1 + 10 * std::numeric_limits<double>::epsilon());

                    double nnx = nmin;
                    double nnz = nmin;
                    
                    double vpmin=model.get_max("vp", rmin, rmax);

                    double vmin = vpmin;
                    double wavemin = vmin/fmax;

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double r0 = std::sqrt(std::pow(current_element.xa(ii,jj), 2) + std::pow(current_element.za(ii,jj), 2));
                            double r = std::min(std::max(rmin, r0), rmax);

                            nnz = std::max(nnz, std::ceil((dtheta*r)/wavemin*ppw));
                            nnx = std::max(nnx, std::ceil((rmax-rmin)/wavemin*ppw));

                        }
                    }

                    current_element.NAX = nnx;
                    current_element.NAZ = nnz;
                    // std::cout << "region 3,4,6,7, CPU:" << icpu<< std::endl;
                    current_element.owner_cpu = icpu;
                    current_element.update_neighbor_cpus();
                    cpus_data[icpu].element_list.push_back(current_element);

                    auto el_ = cpus_data[icpu].element_list[cpus_data[icpu].element_list.size() - 1]; 
                    break;  
                }
                case 5: // central region
                {
                    // std::cout << "region:" << region << std::endl;
                    idomain++; // same as element id
                    icpu = -1;
                    int iib = ((ib) - (length_rad-1) + 1);
                    int jjb = ((jb) - (length_rad-1) + 1);
                    int ij = (jjb - 1)*npart + iib;
                    int cube_ipart = std::ceil((float)(npart*npart)/NCPUs);
                    icpu = std::ceil((float)ij/cube_ipart) - 1; // cpu indies start from 0
                    if(icpu < 0 || icpu >= NCPUs){
                        std::cout << "ERROR (region 5)! CPU can't be less than zero: "<< icpu << std::endl;
                        return -1;
                    }
                    current_element.global_id = idomain;
                    current_element.region_id = region;
                    current_element.domain_id = region;
                    rx1 = allrad[ib]; // start and end radius of region
                    rx2 = allrad[ib + 1];
                    nx = std::max(nunrefmin, std::ceil((rx2 - rx1)/dxz));

                    rz1 = allrad[jb];
                    rz2 = allrad[jb + 1];
                    nz = std::max(nunrefmin, std::ceil((rz2 - rz1)/dxz));     

                    current_element.inflat=inflat;
                    current_element.rx1 = rx1;
                    current_element.rx2 = rx2;
                    current_element.rz1 = rz1;
                    current_element.rz2 = rz2;
                    current_element.rad0 = rad[0];
                    current_element.rad1 = rad[1];
                    current_element.npart = npart;
                    current_element.r1 = 0.0;
                    current_element.r2 = 0.0;
                    current_element.theta1 = 0.0;
                    current_element.theta2 = 0.0;
                    current_element.coeff = 0.0;

                    current_element.nax = nx; // size of unrefined element
                    current_element.naz = nz;

                    current_element.xa.resize(nx,nz); // array containing the spline values
                    current_element.za.resize(nx,nz);  

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            xs = (double)ii/((double)(nx-1.0)) * (rx2 - rx1) + rx1;
                            zs = (double)jj/((double)(nz-1.0)) * (rz2 - rz1) + rz1;
                            
                            if((zs <= -xs) && (zs <= xs)){
                                theta = (xs - zs)/std::abs(zs) *(M_PI/4) + 5*(M_PI/4);
                                r = std::abs(zs);
                            }else if((xs <= -zs) && (xs <= zs)) {
                                theta = (xs - zs)/std::abs(xs) *(M_PI/4) + 5*(M_PI/4);
                                r = std::abs(xs);
                            }else if((xs >= -zs) && (xs >= zs)) {
                                theta = (zs - xs)/std::abs(xs) *(M_PI/4) + (M_PI/4);
                                r = std::abs(xs);
                            }else if((zs >= -xs) && (zs >= xs)) {
                                theta = (zs - xs)/std::abs(zs) *(M_PI/4) + (M_PI/4);
                                r = std::abs(zs);
                            }

                            if(std::abs(r) < __DBL_EPSILON__){
                                r = 0.0;
                                theta = 0.0;
                            }
                            alpha = r/rad[0] * inflat;
                            double value_x = (1.0 - alpha) * xs/(std::sqrt(2)) + alpha * r * std::cos(theta);
                            current_element.xa.value_set(ii,jj, value_x);
                            double value_z = (1.0 - alpha) * zs/(std::sqrt(2)) + alpha * r * std::sin(theta);
                            current_element.za.value_set(ii,jj, value_z);                  
                        }
                    }
                    
                    
                    DFDM::matrix<double> rs(nx, nz); // initialize data vector
                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double x_squared = current_element.xa(ii,jj) * current_element.xa(ii,jj);
                            double z_squared = current_element.za(ii,jj) * current_element.za(ii,jj);
                            rs.value_set(ii,jj, std::sqrt(x_squared + z_squared));
                        }
                    }

                    double max_element = *std::max_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmax = max_element * (1 - 10 * std::numeric_limits<double>::epsilon());

                    double min_element = *std::min_element(rs.data_vec.begin(), rs.data_vec.end());
                    double rmin = min_element * (1 + 10 * std::numeric_limits<double>::epsilon());

                    double nnx = nmin;
                    double nnz = nmin;

                    double max_xa = *std::max_element(current_element.xa.data_vec.begin(), current_element.xa.data_vec.end());
                    double min_xa = *std::min_element(current_element.xa.data_vec.begin(), current_element.xa.data_vec.end());
                    double lx = max_xa - min_xa;

                    double max_za = *std::max_element(current_element.za.data_vec.begin(), current_element.za.data_vec.end());
                    double min_za = *std::min_element(current_element.za.data_vec.begin(), current_element.za.data_vec.end());
                    double lz = max_za - min_za;

                    
                    double vpmin=model.get_max("vp", rmin, rmax);

                    double vmin = vpmin;
                    double wavemin = vmin/fmax;

                    for(uint32_t jj = 0; jj < nz; jj++){
                        for(uint32_t ii = 0; ii < nx; ii++){
                            double r0 = std::sqrt(std::pow(current_element.xa(ii,jj), 2) + std::pow(current_element.za(ii,jj), 2));
                            double r = std::min(std::max(rmin, r0), rmax);

                            nnx = std::max(nnx, std::ceil(lx/wavemin*ppw));
                            nnz = std::max(nnz, std::ceil(lz/wavemin*ppw));
                            // if(nnx > nmin || nnz > nmin){
                            //     std::cout << "Error: nnx or nnz is zero!" << std::endl;
                            //     std::cout << "wavemin: " << wavemin << ", ppw:" << ppw << ", vpmin:"<< vpmin << ", rmin: " << rmin << ", rmax: " << rmax << std::endl;
                            //     exit(-1);
                            // }

                        }
                    }

                    current_element.NAX = nnx;
                    current_element.NAZ = nnz;

                    current_element.owner_cpu = icpu;
                    current_element.update_neighbor_cpus();
                    cpus_data[icpu].element_list.push_back(current_element);

                    auto el_ = cpus_data[icpu].element_list[cpus_data[icpu].element_list.size() - 1];     
                    break;
                }
                default:
                    break;
            }                           
        }

    }

    Mesh::loadBalance(cpus_data, balanced_cpus_data);
    // combining elements to find neighbors.
    Mesh::findNeighbors(balanced_cpus_data);
   // finding neighbor faces:
    Mesh::findFaces(balanced_cpus_data);
    // check validity of mesh:
    Mesh::checkValidityNbrs(balanced_cpus_data);
    std::string grid_output_dir = output_dir + "/grid";
    Mesh::printGrid(balanced_cpus_data, grid_output_dir);
    

    // Print out the size of balanced_cpus_data
    std::cout << "Size of balanced_cpus_data: " << balanced_cpus_data.size() << std::endl;

    // Print out the size of element_list within each element of balanced_cpus_data
    // and the process_id for each balanced_cpus_data element
    for (size_t i = 0; i < balanced_cpus_data.size(); ++i) {
        std::cout << "CPU " << i << " (process_id: " << balanced_cpus_data[i].process_id << "): "
                  << "Size of element_list: " << balanced_cpus_data[i].element_list.size() << std::endl;
    }

    // Print data for each process in a separate file
    for (size_t i = 0; i < balanced_cpus_data.size(); ++i) {
        std::string filename = output_dir + "/process_" + std::to_string(balanced_cpus_data[i].process_id) + "_data.txt";
        
        balanced_cpus_data[i].print_data(filename);
        std::cout << "Data for CPU " << i << " (process_id: " << balanced_cpus_data[i].process_id 
                  << ") written to "<< filename << std::endl;
    }

}