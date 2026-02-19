#include "grid.hpp"
#include "model.hpp"
#include <cmath>

void DFDM::Grid::grid_init(const DFDM::ElementData& el_data, const DFDM::simulation sim, std::shared_ptr<spdlog::logger> logger_){
    logger = logger_;
    element_id = el_data.global_id;
    Nxt = el_data.nax;
    Nzt = el_data.naz;
    Nx1 = el_data.NAX;
    Nz1 = el_data.NAZ;
    Nx2 = el_data.NAX - 1;
    Nz2 = el_data.NAZ - 1;

    x2d11t = el_data.xa;
    z2d11t = el_data.za;

    inflat = el_data.inflat;
    region_id = el_data.region_id;
    theta1 = el_data.theta1;
    theta2 = el_data.theta2;
    r1 = el_data.r1;
    r2 = el_data.r2;
    rx1 = el_data.rx1;
    rx2 = el_data.rx2;
    rz1 = el_data.rz1;
    rz2 = el_data.rz2;
    rad0 = el_data.rad0;
    rad1 = el_data.rad1;
    coeff = el_data.coeff;
    npart = el_data.npart;
    

}

void DFDM::Grid::refine(int p, int pt, DFDM::model model , uint32_t ord_gi, DFDM::Operators& xops, DFDM::Operators& zops){

    logger->debug("Grid::refine:: Refining element:{} ", element_id);
    int kpt = pt+1;

    double dx = 1.0 / (Nx1 - 1);
    double dz = 1.0 / (Nz1 - 1);
    std::vector<double> xps1(Nx1);
    std::vector<double> zps1(Nz1);
    //iota fills the vector starting with 0 and the incrementing by 1.
    std::iota(xps1.begin(), xps1.end(), 0);
    std::iota(zps1.begin(), zps1.end(), 0);

    std::for_each(xps1.begin(), xps1.end(), [dx](double& x){ x *= dx; });
    std::for_each(zps1.begin(), zps1.end(), [dz](double& z){ z *= dz; });
    logger->debug("Grid::refine:: xps1 and zps1 generated, xps1.size:{}, zps1.size:{}", xps1.size(), zps1.size());
    /*
    here assume that x2d11t = xa (from grid init), this needs to be added
    */
    // uint32_t Nxt = x2d11t.rows;
    // uint32_t Nzt = x2d11t.cols;
    
    double dxt = 1.0 / (Nxt - 1);
    double dzt = 1.0 / (Nzt - 1);
    std::vector<double> xpst1(Nxt);
    std::vector<double> zpst1(Nzt);
    std::iota(xpst1.begin(), xpst1.end(), 0);
    std::iota(zpst1.begin(), zpst1.end(), 0);

    std::for_each(xpst1.begin(), xpst1.end(), [dxt](double& x){ x *= dxt; });
    std::for_each(zpst1.begin(), zpst1.end(), [dzt](double& z){ z *= dzt; });
    logger->debug("Grid::refine:: xpst1 and zpst1 generated, xpst1.size:{}, zpst1.size:{}", xpst1.size(), zpst1.size());
    logger->debug("Grid::refine:: x2d11t.rows:{}, x2d11t.cols:{}", x2d11t.rows, x2d11t.cols);
    logger->debug("Grid::refine:: dxt:{}, dzt:{}", dxt, dzt);
    std::vector<double> txt1 = DFDM::utils::get_knot_vector_shape(Nxt, kpt);
    std::vector<double> tzt1 = DFDM::utils::get_knot_vector_shape(Nzt, kpt);
    logger->debug("Grid::refine:: txt1 and tzt1 generated, txt1.size:{}, tzt1.size:{}", txt1.size(), tzt1.size());
    //matrices are initialized to zero by default 
    DFDM::matrix<double> B1xp_mat(Nxt, Nxt);
    DFDM::matrix<double> B1zp_mat(Nzt, Nzt);

    
    for(uint32_t ibx = 0; ibx < Nxt; ibx++){
        for(uint32_t jpx = 0; jpx < Nxt; jpx++){
            double x0 = xpst1[jpx];
            B1xp_mat.value_set(jpx, ibx, DFDM::utils::bspline_c(txt1, Nxt, ibx, kpt, x0));
        }
    }

    for(uint32_t ibz = 0; ibz < Nzt; ibz++){
        for(uint32_t jpz = 0; jpz < Nzt; jpz++){
            double z0 = zpst1[jpz];
            B1zp_mat.value_set(jpz, ibz, DFDM::utils::bspline_c(tzt1, Nzt, ibz, kpt, z0));
        }
    }
    logger->debug("Grid::refine:: B1xp_mat and B1zp_mat generated");
    DFDM::matrix<double> invB1xp_mat;
    B1xp_mat.inverse(invB1xp_mat);

    DFDM::matrix<double> invB1zp_mat;
    B1zp_mat.inverse(invB1zp_mat);
    logger->debug("Grid::refine:: invB1xp_mat and invB1zp_mat generated");
    // The corresponding coeffs of the grids in x and z direction
    // Note that only need to invert one direction
    // x direction
    DFDM::matrix<double> coef11x_xp = invB1xp_mat.matprod(x2d11t);
    DFDM::matrix<double> coef11z_xp = invB1xp_mat.matprod(z2d11t);
    logger->debug("Grid::refine:: coef11x_xp and coef11z_xp generated");

    // z direction (need to take transpose for this direction)
    DFDM::matrix<double> x2d11tT = x2d11t.transpose_inplace();
    DFDM::matrix<double> coef11x_zp = invB1zp_mat.matprod(x2d11tT);
    coef11x_zp = coef11x_zp.transpose_inplace();
    
    logger->debug("Grid::refine:: coef11x_zp generated");
    DFDM::matrix<double> z2d11tT = z2d11t.transpose_inplace();
    DFDM::matrix<double> coef11z_zp = invB1zp_mat.matprod(z2d11tT);
    coef11z_zp = coef11z_zp.transpose_inplace();
    logger->debug("Grid::refine:: coef11z_zp generated");

    // Plotting matrices for x and z directions
    DFDM::matrix<double> B1xp_mat_plot(Nx1, Nxt);
    DFDM::matrix<double> B1zp_mat_plot(Nz1, Nzt);

    for(uint32_t jpx = 0; jpx < Nx1; jpx++){
        double x0 = xps1[jpx];
        for(uint32_t ib1x = 0; ib1x < Nxt; ib1x++){
            B1xp_mat_plot.value_set(jpx, ib1x, DFDM::utils::bspline_c(txt1, Nxt, ib1x, kpt, x0));
        }
    }

    for(uint32_t jpz = 0; jpz < Nz1; jpz++){
        double z0 = zps1[jpz];
        for(uint32_t ib1z = 0; ib1z < Nzt; ib1z++){
            B1zp_mat_plot.value_set(jpz, ib1z, DFDM::utils::bspline_c(tzt1, Nzt, ib1z, kpt, z0));
        }
    }
    logger->debug("Grid::refine:: B1xp_mat_plot and B1zp_mat_plot generated");
    // Calculating partial derivatives at the hat points in the reference element
    // x direction
    x2d11 = B1xp_mat_plot.matprod(coef11x_xp);
    auto x2d11T = x2d11.transpose_inplace();
    x2d11T = B1zp_mat_plot.matprod(x2d11T);
    x2d11 = x2d11T.transpose_inplace();
    logger->debug("Grid::refine:: x2d11 generated");
    // z direction
    z2d11 = B1xp_mat_plot.matprod(coef11z_zp);
    auto z2d11T = z2d11.transpose_inplace();
    z2d11T = B1zp_mat_plot.matprod(z2d11T);
    z2d11 = z2d11T.transpose_inplace();
    logger->debug("Grid::refine:: z2d11 generated");
    // x2d11.print_file("x2d11_"+std::to_string(element_id));
    // z2d11.print_file("z2d11_"+std::to_string(element_id));
    // skipping the plotting lines from 79 to 102.
    // looks like hat points can be generated along with operator matrices
    // and then passed here, operator matrices do not need refined points.
    std::vector<double> hatpointsx1, hatpointsx2, hatpointsz1, hatpointsz2;
    hat_points(hatpointsx1, hatpointsx2, xops, Nx1, p, ord_gi);
    logger->debug("Grid::refine:: hatpointsx1 and hatpointsx2 generated");
    hat_points(hatpointsz1, hatpointsz2, zops, Nz1, p, ord_gi);
    logger->debug("Grid::refine:: hatpointsz1 and hatpointsz2 generated");

    /************* ***************/

    
    // if(element_id == 0){
    // std::cout << " txt1: ";
    // for(uint32_t h = 0; h < txt1.size(); h++){
    //     std::cout << txt1[h] << "  ";
    // }
    // std::cout << std::endl;

    // std::cout << " tzt1: ";
    // for(uint32_t h = 0; h < tzt1.size(); h++){
    //     std::cout << tzt1[h] << "  ";
    // }
    // std::cout << std::endl;

    // std::cout << " hatpointsz1: ";
    // for(uint32_t h = 0; h < hatpointsz1.size(); h++){
    //     std::cout << hatpointsz1[h] << "  ";
    // }
    // std::cout << std::endl;

    // std::cout << " hatpointsz2: ";
    // for(uint32_t h = 0; h < hatpointsz2.size(); h++){
    //     std::cout << hatpointsz2[h] << "  ";
    // }
    // std::cout << std::endl;
    // }
    /************* ***************/
    // calculating the bspline basis values at hat points in x and z dir

    DFDM::matrix<double> B1xp_mat_hat1(Nx1, Nxt);
    DFDM::matrix<double> B1xp_mat_hat2(Nx2, Nxt);
    DFDM::matrix<double> B1zp_mat_hat1(Nz1, Nzt);
    DFDM::matrix<double> B1zp_mat_hat2(Nz2, Nzt);

    for(uint32_t ib1x = 0; ib1x < Nxt; ib1x++){
        for(uint32_t jpx = 0; jpx < Nx1; jpx++){
            double x0 = hatpointsx1[jpx];
            B1xp_mat_hat1.value_set(jpx, ib1x, DFDM::utils::bspline_c(txt1, Nxt, ib1x, kpt, x0));
        }
    }
    logger->debug("Grid::refine:: B1xp_mat_hat1 generated");
    for(uint32_t ib1x = 0; ib1x < Nxt; ib1x++){
        for(uint32_t jpx = 0; jpx < Nx2; jpx++){
            double x0 = hatpointsx2[jpx];
            B1xp_mat_hat2.value_set(jpx, ib1x, DFDM::utils::bspline_c(txt1, Nxt, ib1x, kpt, x0));
        }
    }
    logger->debug("Grid::refine:: B1xp_mat_hat2 generated");
    for(uint32_t ib1z = 0; ib1z < Nzt; ib1z++){
        for(uint32_t jpz = 0; jpz < Nz1; jpz++){
            double z0 = hatpointsz1[jpz];
            B1zp_mat_hat1.value_set(jpz, ib1z, DFDM::utils::bspline_c(tzt1, Nzt, ib1z, kpt, z0));
        }
    }
    logger->debug("Grid::refine:: B1zp_mat_hat1 generated");

    for(uint32_t ib1z = 0; ib1z < Nzt; ib1z++){
        for(uint32_t jpz = 0; jpz < Nz2; jpz++){
            double z0 = hatpointsz2[jpz];
            B1zp_mat_hat2.value_set(jpz, ib1z, DFDM::utils::bspline_c(tzt1, Nzt, ib1z, kpt, z0));
        }
    }
    logger->debug("Grid::refine:: B1zp_mat_hat2 generated");

    // calculating partial derivations at the hat points in the reference element
    DFDM::matrix<double> dB1xp_mat_hat1(Nx1, Nxt);
    DFDM::matrix<double> dB1xp_mat_hat2(Nx2, Nxt);
    DFDM::matrix<double> dB1zp_mat_hat1(Nz1, Nzt);
    DFDM::matrix<double> dB1zp_mat_hat2(Nz2, Nzt);

    for(uint32_t ib1x = 0; ib1x < Nxt; ib1x++){
        for(uint32_t jpx = 0; jpx < Nx1; jpx++){
            double x0 = hatpointsx1[jpx];
            dB1xp_mat_hat1.value_set(jpx, ib1x, DFDM::utils::dbspline_c(txt1, Nxt, ib1x, kpt, x0, 1));
        }
    }
    logger->debug("Grid::refine:: dB1xp_mat_hat1 generated");

    for(uint32_t ib1x = 0; ib1x < Nxt; ib1x++){
        for(uint32_t jpx = 0; jpx < Nx2; jpx++){
            double x0 = hatpointsx2[jpx];
            dB1xp_mat_hat2.value_set(jpx, ib1x, DFDM::utils::dbspline_c(txt1, Nxt, ib1x, kpt, x0, 1));
        }
    }
    logger->debug("Grid::refine:: dB1xp_mat_hat2 generated");

    for(uint32_t ib1z = 0; ib1z < Nzt; ib1z++){
        for(uint32_t jpz = 0; jpz < Nz1; jpz++){
            double z0 = hatpointsz1[jpz];
            dB1zp_mat_hat1.value_set(jpz, ib1z, DFDM::utils::dbspline_c(tzt1, Nzt, ib1z, kpt, z0, 1));
        }
    }
    logger->debug("Grid::refine:: dB1zp_mat_hat1 generated");

    for(uint32_t ib1z = 0; ib1z < Nzt; ib1z++){
        for(uint32_t jpz = 0; jpz < Nz2; jpz++){
            double z0 = hatpointsz2[jpz];
            dB1zp_mat_hat2.value_set(jpz, ib1z, DFDM::utils::dbspline_c(tzt1, Nzt, ib1z, kpt, z0, 1));
        }
    }
    logger->debug("Grid::refine:: dB1zp_mat_hat2 generated");

    /*
    matlab plotting code from 172 to 191, not converting to C++:
           % % using the coeffs obtained above to get the coordinates of hat points in the reference element
        % 11
        x2d11_hat  = pagemtimes(B1xp_mat_hat1,coef11x_xp);
        x2d11_hatT = permute(x2d11_hat,[2,1,3]);
        x2d11_hatT = pagemtimes(B1zp_mat_hat1, x2d11_hatT);
        x2d11_hat  = permute(x2d11_hatT,[2,1,3]);

        z2d11_hat  = pagemtimes(B1xp_mat_hat1,coef11z_zp);
        z2d11_hatT = permute(z2d11_hat,[2,1,3]);
        z2d11_hatT = pagemtimes(B1zp_mat_hat1, z2d11_hatT);
        z2d11_hat  = permute(z2d11_hatT,[2,1,3]);

        
        hold on
        scatter(x2d11_hat(1:Nx1*Nz1)/1e3,z2d11_hat(1:Nx1*Nz1)/1e3,'b.');
        % scatter(x2d12_hat(1:Nx1*Nz2)/1e3,z2d12_hat(1:Nx1*Nz2)/1e3,'g.');
        % scatter(x2d21_hat(1:Nx2*Nz1)/1e3,z2d21_hat(1:Nx2*Nz1)/1e3,'b.');
        % scatter(x2d22_hat(1:Nx2*Nz2)/1e3,z2d22_hat(1:Nx2*Nz2)/1e3,'m.');
        axis equal tight
        hold off
    */
    // calculating the spatial derivatives
    // x direction 
    // dxdxp11 
    // if(element_id == 0) coef11x_xp.print_file("coef11x_xp");
    // if(element_id == 0) coef11z_xp.print_file("coef11z_xp");
    // if(element_id == 0) coef11z_zp.print_file("coef11z_zp");
    // if(element_id == 0) coef11x_zp.print_file("coef11x_zp");
    auto dxdxp11_hat = dB1xp_mat_hat1.matprod(coef11x_xp);
    auto dxdxp11_hatT = dxdxp11_hat.transpose_inplace();
    auto dxdxp11_hatT2 = B1zp_mat_hat1.matprod(dxdxp11_hatT);
    dxdxp11_hat = dxdxp11_hatT2.transpose_inplace();
    // if(element_id == 0) dxdxp11_hat.print_file("dxdxp11_hat");
    logger->debug("Grid::refine:: dxdxp11_hat generated");
    // dxdxp12
    auto dxdxp12_hat = dB1xp_mat_hat1.matprod(coef11x_xp);
    auto dxdxp12_hatT = dxdxp12_hat.transpose_inplace();
    auto dxdxp12_hatT2 = B1zp_mat_hat2.matprod(dxdxp12_hatT);
    dxdxp12_hat = dxdxp12_hatT2.transpose_inplace();
    // if(element_id == 0) dxdxp12_hat.print_file("dxdxp12_hat");
    logger->debug("Grid::refine:: dxdxp12_hat generated");
    // dxdxp21
    auto dxdxp21_hat = dB1xp_mat_hat2.matprod(coef11x_xp);
    auto dxdxp21_hatT = dxdxp21_hat.transpose_inplace();
    auto dxdxp21_hatT2 = B1zp_mat_hat1.matprod(dxdxp21_hatT);
    dxdxp21_hat = dxdxp21_hatT2.transpose_inplace();
    // if(element_id == 0) dxdxp21_hat.print_file("dxdxp21_hat");
    logger->debug("Grid::refine:: dxdxp21_hat generated");
    // dxdxp22  
    auto dxdxp22_hat = dB1xp_mat_hat2.matprod(coef11x_xp);
    auto dxdxp22_hatT = dxdxp22_hat.transpose_inplace();
    auto dxdxp22_hatT2 = B1zp_mat_hat2.matprod(dxdxp22_hatT);
    dxdxp22_hat = dxdxp22_hatT2.transpose_inplace();
    // if(element_id == 0) dxdxp22_hat.print_file("dxdxp22_hat");
    logger->debug("Grid::refine:: dxdxp22_hat generated");
    // dzdxp11
    auto dzdxp11_hat = dB1xp_mat_hat1.matprod(coef11z_xp);
    auto dzdxp11_hatT = dzdxp11_hat.transpose_inplace();
    auto dzdxp11_hatT2 = B1zp_mat_hat1.matprod(dzdxp11_hatT);
    dzdxp11_hat = dzdxp11_hatT2.transpose_inplace();
    // if(element_id == 0) dzdxp11_hat.print_file("dzdxp11_hat");
    logger->debug("Grid::refine:: dzdxp11_hat generated");
    // dzdxp12
    auto dzdxp12_hat = dB1xp_mat_hat1.matprod(coef11z_xp);
    auto dzdxp12_hatT = dzdxp12_hat.transpose_inplace();
    auto dzdxp12_hatT2 = B1zp_mat_hat2.matprod(dzdxp12_hatT);
    dzdxp12_hat = dzdxp12_hatT2.transpose_inplace();
    // if(element_id == 0) dzdxp12_hat.print_file("dzdxp12_hat");
    logger->debug("Grid::refine:: dzdxp12_hat generated");
    // dzdxp21
    auto dzdxp21_hat = dB1xp_mat_hat2.matprod(coef11z_xp);
    auto dzdxp21_hatT = dzdxp21_hat.transpose_inplace();
    auto dzdxp21_hatT2 = B1zp_mat_hat1.matprod(dzdxp21_hatT);
    dzdxp21_hat = dzdxp21_hatT2.transpose_inplace();
    // if(element_id == 0) dzdxp21_hat.print_file("dzdxp21_hat");
    logger->debug("Grid::refine:: dzdxp21_hat generated");
    // dzdxp22
    auto dzdxp22_hat = dB1xp_mat_hat2.matprod(coef11z_xp);
    auto dzdxp22_hatT = dzdxp22_hat.transpose_inplace();
    auto dzdxp22_hatT2 = B1zp_mat_hat2.matprod(dzdxp22_hatT);
    dzdxp22_hat = dzdxp22_hatT2.transpose_inplace();
    // if(element_id == 0) dzdxp22_hat.print_file("dzdxp22_hat");
    logger->debug("Grid::refine:: dzdxp22_hat generated");
    // z direction
    // dxdzp11
    auto dxdzp11_hat = B1xp_mat_hat1.matprod(coef11x_zp);
    auto dxdzp11_hatT = dxdzp11_hat.transpose_inplace();
    auto dxdzp11_hatT2 = dB1zp_mat_hat1.matprod(dxdzp11_hatT);
    dxdzp11_hat = dxdzp11_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dxdzp11_hat generated");
    // dzdzp11
    auto dzdzp11_hat = B1xp_mat_hat1.matprod(coef11z_zp);
    auto dzdzp11_hatT = dzdzp11_hat.transpose_inplace();
    auto dzdzp11_hatT2 = dB1zp_mat_hat1.matprod(dzdzp11_hatT);
    dzdzp11_hat = dzdzp11_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dzdzp11_hat generated");
    // dxdzp21
    auto dxdzp21_hat = B1xp_mat_hat2.matprod(coef11x_zp);
    auto dxdzp21_hatT = dxdzp21_hat.transpose_inplace();
    auto dxdzp21_hatT2 = dB1zp_mat_hat1.matprod(dxdzp21_hatT);
    dxdzp21_hat = dxdzp21_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dxdzp21_hat generated");
    // dzdzp21
    auto dzdzp21_hat = B1xp_mat_hat2.matprod(coef11z_zp);
    auto dzdzp21_hatT = dzdzp21_hat.transpose_inplace();
    auto dzdzp21_hatT2 = dB1zp_mat_hat1.matprod(dzdzp21_hatT);
    dzdzp21_hat = dzdzp21_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dzdxp22_hat generated");
    // dxdzp12
    auto dxdzp12_hat = B1xp_mat_hat1.matprod(coef11x_zp);
    auto dxdzp12_hatT = dxdzp12_hat.transpose_inplace();
    auto dxdzp12_hatT2 = dB1zp_mat_hat2.matprod(dxdzp12_hatT);
    dxdzp12_hat = dxdzp12_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dxdzp12_hat generated");
    // dzdzp12
    auto dzdzp12_hat = B1xp_mat_hat1.matprod(coef11z_zp);
    auto dzdzp12_hatT = dzdzp12_hat.transpose_inplace();
    auto dzdzp12_hatT2 = dB1zp_mat_hat2.matprod(dzdzp12_hatT);
    dzdzp12_hat = dzdzp12_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dzdzp12_hat generated");
    // dxdzp22
    auto dxdzp22_hat = B1xp_mat_hat2.matprod(coef11x_zp);
    auto dxdzp22_hatT = dxdzp22_hat.transpose_inplace();
    auto dxdzp22_hatT2 = dB1zp_mat_hat2.matprod(dxdzp22_hatT);
    dxdzp22_hat = dxdzp22_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dxdzp22_hat generated");
    // dzdzp22
    auto dzdzp22_hat = B1xp_mat_hat2.matprod(coef11z_zp);
    auto dzdzp22_hatT = dzdzp22_hat.transpose_inplace();
    auto dzdzp22_hatT2 = dB1zp_mat_hat2.matprod(dzdzp22_hatT);
    dzdzp22_hat = dzdzp22_hatT2.transpose_inplace();
    logger->debug("Grid::refine:: dzdzp22_hat generated");


   // Jac
    auto Jac11_hat = DFDM::matrix<double>::elementwise_mult(dxdxp11_hat, dzdzp11_hat) - DFDM::matrix<double>::elementwise_mult(dzdxp11_hat, dxdzp11_hat);
    logger->debug("Grid::refine:: Jac11_hat generated");
    auto Jac12_hat = DFDM::matrix<double>::elementwise_mult(dxdxp12_hat, dzdzp12_hat) - DFDM::matrix<double>::elementwise_mult(dzdxp12_hat, dxdzp12_hat);
    logger->debug("Grid::refine:: Jac12_hat generated");
    auto Jac21_hat = DFDM::matrix<double>::elementwise_mult(dxdxp21_hat, dzdzp21_hat) - DFDM::matrix<double>::elementwise_mult(dzdxp21_hat, dxdzp21_hat);
    logger->debug("Grid::refine:: Jac21_hat generated");
    auto Jac22_hat = DFDM::matrix<double>::elementwise_mult(dxdxp22_hat, dzdzp22_hat) - DFDM::matrix<double>::elementwise_mult(dzdxp22_hat, dxdzp22_hat);
    logger->debug("Grid::refine:: Jac22_hat generated");

    // if(element_id == 0) Jac11_hat.print_file("Jac11_hat");
    // if(element_id == 0) Jac12_hat.print_file("Jac12_hat");
    // if(element_id == 0) Jac21_hat.print_file("Jac21_hat");
    // if(element_id == 0) Jac22_hat.print_file("Jac22_hat");
    //define a vactor and fill its data_vec with 1s, assuming same dims for all
    
    // 1/Jac
    DFDM::matrix<double> ones_Jac11_hat(Jac11_hat.rows, Jac11_hat.cols);
    std::fill(ones_Jac11_hat.data_vec.begin(), ones_Jac11_hat.data_vec.end(), 1);
    logger->debug("Grid::refine:: ones_Jac11_hat generated");
    auto tmp11_hat = DFDM::matrix<double>::elementwise_div(ones_Jac11_hat, Jac11_hat);
    
    DFDM::matrix<double> ones_Jac12_hat(Jac12_hat.rows, Jac12_hat.cols);
    std::fill(ones_Jac12_hat.data_vec.begin(), ones_Jac12_hat.data_vec.end(), 1);
    logger->debug("Grid::refine:: tmp11_hat generated");
    auto tmp12_hat = DFDM::matrix<double>::elementwise_div(ones_Jac12_hat, Jac12_hat);
    logger->debug("Grid::refine:: tmp12_hat generated");

    DFDM::matrix<double> ones_Jack21_hat(Jac21_hat.rows, Jac21_hat.cols);
    std::fill(ones_Jack21_hat.data_vec.begin(), ones_Jack21_hat.data_vec.end(), 1);
    auto tmp21_hat = DFDM::matrix<double>::elementwise_div(ones_Jack21_hat, Jac21_hat);
    logger->debug("Grid::refine:: tmp21_hat generated");

    DFDM::matrix<double> ones_Jac22_hat(Jac22_hat.rows, Jac22_hat.cols);
    std::fill(ones_Jac22_hat.data_vec.begin(), ones_Jac22_hat.data_vec.end(), 1);   
    auto tmp22_hat = DFDM::matrix<double>::elementwise_div(ones_Jac22_hat, Jac22_hat);
    logger->debug("Grid::refine:: tmp22_hat generated");

    // dxi_dx
    auto dxpdx11_hat = DFDM::matrix<double>::elementwise_mult(dzdzp11_hat, tmp11_hat);
    auto dzpdx11_hat = DFDM::matrix<double>::elementwise_mult(dzdxp11_hat, tmp11_hat) * -1 ;
    auto dxpdz11_hat = DFDM::matrix<double>::elementwise_mult(dxdzp11_hat, tmp11_hat) * -1;
    auto dzpdz11_hat = DFDM::matrix<double>::elementwise_mult(dxdxp11_hat, tmp11_hat);

    // if(element_id == 0) dxpdx11_hat.print_file("dxpdx11_hat");
    // if(element_id == 0) dzpdx11_hat.print_file("dzpdx11_hat");
    // if(element_id == 0) dxpdz11_hat.print_file("dxpdz11_hat");
    // if(element_id == 0) dzpdz11_hat.print_file("dzpdz11_hat");

    auto dxpdx12_hat = DFDM::matrix<double>::elementwise_mult(dzdzp12_hat, tmp12_hat);
    auto dzpdx12_hat = DFDM::matrix<double>::elementwise_mult(dzdxp12_hat, tmp12_hat) * -1;
    auto dxpdz12_hat = DFDM::matrix<double>::elementwise_mult(dxdzp12_hat, tmp12_hat) * -1;
    auto dzpdz12_hat = DFDM::matrix<double>::elementwise_mult(dxdxp12_hat, tmp12_hat);

    // if(element_id == 0) dxpdx12_hat.print_file("dxpdx12_hat");
    // if(element_id == 0) dzpdx12_hat.print_file("dzpdx12_hat");
    // if(element_id == 0) dxpdz12_hat.print_file("dxpdz12_hat");
    // if(element_id == 0) dzpdz12_hat.print_file("dzpdz12_hat");

    auto dxpdx21_hat = DFDM::matrix<double>::elementwise_mult(dzdzp21_hat, tmp21_hat);
    auto dzpdx21_hat = DFDM::matrix<double>::elementwise_mult(dzdxp21_hat, tmp21_hat) * -1;
    auto dxpdz21_hat = DFDM::matrix<double>::elementwise_mult(dxdzp21_hat, tmp21_hat) * -1;
    auto dzpdz21_hat = DFDM::matrix<double>::elementwise_mult(dxdxp21_hat, tmp21_hat);

    // if(element_id == 0) dxpdx21_hat.print_file("dxpdx21_hat");
    // if(element_id == 0) dzpdx21_hat.print_file("dzpdx21_hat");
    // if(element_id == 0) dxpdz21_hat.print_file("dxpdz21_hat");
    // if(element_id == 0) dzpdz21_hat.print_file("dzpdz21_hat");

    auto dxpdx22_hat = DFDM::matrix<double>::elementwise_mult(dzdzp22_hat, tmp22_hat);
    auto dzpdx22_hat = DFDM::matrix<double>::elementwise_mult(dzdxp22_hat, tmp22_hat) * -1;
    auto dxpdz22_hat = DFDM::matrix<double>::elementwise_mult(dxdzp22_hat, tmp22_hat) * -1;
    auto dzpdz22_hat = DFDM::matrix<double>::elementwise_mult(dxdxp22_hat, tmp22_hat);

    // if(element_id == 0) dxpdx22_hat.print_file("dxpdx22_hat");
    // if(element_id == 0) dzpdx22_hat.print_file("dzpdx22_hat");
    // if(element_id == 0) dxpdz22_hat.print_file("dxpdz22_hat");
    // if(element_id == 0) dzpdz22_hat.print_file("dzpdz22_hat");

    //setting up class variables
    px1 = p;
    pz1 = p;
    dxpdx11 = dxpdx11_hat;
    dzpdx11 = dzpdx11_hat;
    dxpdz11 = dxpdz11_hat;
    dzpdz11 = dzpdz11_hat;
    dxpdx12 = dxpdx12_hat;
    dzpdx12 = dzpdx12_hat;
    dxpdz12 = dxpdz12_hat;
    dzpdz12 = dzpdz12_hat;
    dxpdx21 = dxpdx21_hat;
    dzpdx21 = dzpdx21_hat;
    dxpdz21 = dxpdz21_hat;
    dzpdz21 = dzpdz21_hat;
    dxpdx22 = dxpdx22_hat;
    dzpdx22 = dzpdx22_hat;
    dxpdz22 = dxpdz22_hat;
    dzpdz22 = dzpdz22_hat;
    Jac12 = Jac12_hat;
    Jac21 = Jac21_hat;
    Jac11 = Jac11_hat;
    Jac22 = Jac22_hat;

    rho12.resize(Nx1, Nz2);
    rho21.resize(Nx2, Nz1);

    mu11.resize(Nx1, Nz1);
    mu22.resize(Nx2, Nz2);

    for (uint32_t i = 0; i < Nx1; i++){
        for(uint32_t j = 0; j < Nz1; j++){
            double hatpointsx11_global, hatpointsz11_global;
            local_to_global(hatpointsx1[i], hatpointsz1[j], hatpointsx11_global, hatpointsz11_global);
            mu11.value_set(i, j, model.get_value("mu", hatpointsx11_global, hatpointsz11_global));
        }
    }

    for (uint32_t i = 0; i < Nx2; i++){
        for(uint32_t j = 0; j < Nz2; j++){
            double hatpointsx22_global, hatpointsz22_global;
            local_to_global(hatpointsx2[i], hatpointsz2[j], hatpointsx22_global, hatpointsz22_global);
            mu22.value_set(i, j, model.get_value("mu", hatpointsx22_global, hatpointsz22_global));
        }
    }

    for (uint32_t i = 0; i < Nx1; i++){
        for(uint32_t j = 0; j < Nz2; j++){
            double hatpointsx12_global, hatpointsz12_global;
            local_to_global(hatpointsx1[i], hatpointsz2[j], hatpointsx12_global, hatpointsz12_global);
            rho12.value_set(i, j, model.get_value("rho", hatpointsx12_global, hatpointsz12_global));
        }
    }

    for (uint32_t i = 0; i < Nx2; i++){
        for(uint32_t j = 0; j < Nz1; j++){
            double hatpointsx21_global, hatpointsz21_global;
            local_to_global(hatpointsx2[i], hatpointsz1[j], hatpointsx21_global, hatpointsz21_global);
            rho21.value_set(i, j, model.get_value("rho", hatpointsx21_global, hatpointsz21_global));
        }
    }
}
void DFDM::Grid::hat_points(std::vector<double>& hatpoints1, std::vector<double>& hatpoints2, DFDM::Operators& ops, uint32_t N, uint32_t p, uint32_t ord_gi){
    // hatpoints are the center of mass of each spline
    // They are calculated for the local splines, not the orthogonalised splines. 

    auto NB_intervals = N - p;
    uint32_t k = p + 1;
    auto t1 = ops.t1;
    auto t2 = ops.t2;
    auto xint = ops.xint;
    auto wxint = ops.wxint;

    auto invl11 = ops.invL11;
    auto invl22 = ops.invL22;
    
    // if(ops.op_dim == DFDM::OP_TYPE::X){
    //     invl11.print_file("invl11_"+std::to_string(element_id));
    //     invl22.print_file("invl22_"+std::to_string(element_id));
    //     std::cout << "N:" << N <<std::endl;
    // }

    DFDM::matrix<double> intB1(N, 1);
    DFDM::matrix<double> intB2(N-1, 1);
    for(uint32_t ib1 = 0; ib1 < N; ib1++){
        for(uint32_t kd = 0; kd < NB_intervals; kd++){
            for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                double b1tmp1 = DFDM::utils::bspline_c(t1, N, ib1, k, xint[kd][lpt]);
                intB1.value_set(ib1, 0, intB1(ib1, 0) + b1tmp1 * wxint[kd][lpt]);
            }
        }
    }
    for(uint32_t ib2 = 0; ib2 < N-1; ib2++){
        for(uint32_t kd = 0; kd < NB_intervals; kd++){
            for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                double b2tmp1 = DFDM::utils::bspline_c(t2, N-1, ib2, k-1, xint[kd][lpt]);
                intB2.value_set(ib2, 0, intB2(ib2, 0) + b2tmp1 * wxint[kd][lpt]);
            }
        }
    }

    DFDM::matrix<double> intxB1(N, 1);
    DFDM::matrix<double> intxB2(N-1, 1);
    for(uint32_t ib1 = 0; ib1 < N; ib1++){
        for(uint32_t kd = 0; kd < NB_intervals; kd++){
            for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                double b1tmp1 = DFDM::utils::bspline_c(t1, N, ib1, k, xint[kd][lpt]);
                double tmp2 = xint[kd][lpt];
                intxB1.value_set(ib1, 0, intxB1(ib1,0)+ b1tmp1 * tmp2 * wxint[kd][lpt]);
            }
        }
    }
    for(uint32_t ib2 = 0; ib2 < N-1; ib2++){
        for(uint32_t kd = 0; kd < NB_intervals; kd++){
            for(uint32_t lpt = 0; lpt < ord_gi; lpt++){
                double b2tmp1 = DFDM::utils::bspline_c(t2, N-1, ib2, k-1, xint[kd][lpt]);
                double tmp2 = xint[kd][lpt];
                intxB2.value_set(ib2, 0, intxB2(ib2,0) + b2tmp1 * tmp2 * wxint[kd][lpt]);
            }
        }
    }

    //  if(ops.op_dim == DFDM::OP_TYPE::X){
    //     intB1.print_file("intB1_"+std::to_string(element_id));
    //     intB2.print_file("intB2_"+std::to_string(element_id));
    //     intxB1.print_file("intxB1_"+std::to_string(element_id));
    //     intxB2.print_file("intxB2_"+std::to_string(element_id));
    //  }

    auto numer1 = invl11.matprod(intxB1);
    auto denom1 = invl11.matprod(intB1);
    hatpoints1.resize(numer1.data_vec.size());
    for(uint32_t i = 0; i < numer1.data_vec.size(); i++){
        hatpoints1[i] = numer1.data_vec[i] / denom1.data_vec[i];
    }

    auto numer2 = invl22.matprod(intxB2);
    auto denom2 = invl22.matprod(intB2);
    hatpoints2.resize(numer2.data_vec.size());
    for(uint32_t i = 0; i < numer2.data_vec.size(); i++){
        hatpoints2[i] = numer2.data_vec[i] / denom2.data_vec[i];
    }

// below modifications is for when using cholesky
    // hatpoints1[hatpoints1.size() - 1] = 1;
    // hatpoints2[hatpoints2.size() - 1] = 1;
}

void DFDM::Grid::local_to_global(double x_local, double z_local, double& x_global, double& z_global){
    
    if (x_local > 1. || x_local < 0.0 || z_local > 1.0 || z_local < 0.0){
        throw std::runtime_error("x_local and z_local must be between 0 and 1");
    }

    if (region_id == 1 || region_id == 2 || region_id == 8 || region_id == 9){
        double aR = std::abs(r2+r1)/2;
        double dR = std::abs(r2-r1);

        double theta = x_local * (theta2 - theta1) + theta1;
        double r = z_local * (r2-r1) + r1;
        double rabs = std::abs(r);

        double x_tmp = rabs*std::cos(theta);
        double z_tmp = rabs*std::sin(theta);

        double zs = r/std::sqrt(2);

        double coeff_tmp;
        if (region_id == 1 || region_id == 2){
            coeff_tmp = coeff - 2.0*x_local/npart;
        }if(region_id == 8 || region_id == 9){
            coeff_tmp = coeff + 2.0*x_local/npart;
        }

        double xs = coeff_tmp*zs;

        if(region_id == 2 || region_id == 8){
            double alpha = (1.0 - inflat) * (1.0 - (rabs - rad0)/(rad1 - rad0));
            x_global = alpha * xs + (1.0 - alpha) * x_tmp;
            z_global = alpha * zs + (1.0 - alpha) * z_tmp;
        }else{
            x_global = x_tmp;
            z_global = z_tmp;
        }
    }
    else if (region_id == 3 || region_id == 4 || region_id == 6 || region_id == 7){
        double aR = std::abs(r2+r1)/2;
        double dR = std::abs(r2-r1);

        double theta = z_local * (theta2 - theta1) + theta1;
        double r = x_local * (r2-r1) + r1;
        double rabs = std::abs(r);

        double x_tmp = rabs*std::cos(theta);
        double z_tmp = rabs*std::sin(theta);

        double xs = r/std::sqrt(2);

        double coeff_tmp;
        if (region_id == 6 || region_id == 7){
            coeff_tmp = coeff + 2*z_local/npart;
        }if(region_id == 3 || region_id == 4){
            coeff_tmp = coeff - 2*z_local/npart;
        }

        double zs = coeff_tmp*xs;

        if(region_id == 4 || region_id == 6){
            double alpha = (1 - inflat) * (1 - (rabs - rad0)/(rad1 - rad0));
            x_global = alpha * xs + (1 - alpha) * x_tmp;
            z_global = alpha * zs + (1 - alpha) * z_tmp;
        }else{
            x_global = x_tmp;
            z_global = z_tmp;
        }
    }
    else if (region_id == 5){
        double xs = x_local * (rx2 - rx1) + rx1;
        double zs = z_local * (rz2 - rz1) + rz1;

        double r, theta;
        
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
        double alpha = r/rad0 * inflat;
        x_global = (1.0 - alpha) * xs/(std::sqrt(2)) + alpha * r * std::cos(theta);
        z_global = (1.0 - alpha) * zs/(std::sqrt(2)) + alpha * r * std::sin(theta);
    }
    else{
        throw std::runtime_error("Invalid region_id: " + std::to_string(region_id));
    }
}