#include "utils.hpp"

double DFDM::utils::max2d(DFDM::matrix<double>& in_){
    double max = 0;
        for(uint32_t i = 0; i < in_.rows; i++){
            for(uint32_t j = 0; j < in_.cols; j++){
                if (max < in_(i,j)){
                    max = in_(i,j);
                }
            }
        }
    return max;

}

double DFDM::utils::dbspline_c(const std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo, uint32_t d_order){
    double b = 0;
    double c1 = 0 ,c2 = 0;
    if(d_order == 0){
        b = bspline_c(t_,n,i,k,xo);
    }else if(d_order > 0 && d_order < k){
        if( std::abs(t_[i+k-1] - t_[i]) > EPS){
            c1 = (k-1)/(t_[i+k-1] - t_[i]); //
        }else{
            c1 = 0;
        }

        if(std::abs(t_[i+k] - t_[i+1] > EPS)){
            c2 = (k-1)/(t_[i+k]-t_[i+1]); // 
        }else{
            c2 = 0;
        }
        b = c1*dbspline_c(t_,n,i,k-1,xo,d_order-1)-c2*dbspline_c(t_,n,i+1,k-1,xo,d_order-1);
    }
    return b;
}

/*
t=kp-1;

t  = zeros(np+kp,1);
dx = 1/(np+1-kp);


xis = 0:1/(np-1):1;

for i = kp:np+1
    t(i) = (i-kp)*dx;
end

if mod(pt,2)==0
    for ii = pt+2:np
        id1 = ii - pt/2 - 1;
        id2 = ii - pt/2;
        t(ii) = ( xis(id1) + xis(id2) )/2;
    end 
else
    for ii = pt+2:np
        id1   = ii - (pt+1)/2;
        t(ii) = xis(id1);
    end
end

t(1:kp)       = t(kp)   - 10*eps;
t(np+1:np+kp) = t(np+1) + 10*eps;

end


*/

std::vector<double> DFDM::utils::get_knot_vector_shape(uint32_t np, uint32_t kp){
    uint32_t pt = kp - 1;
    std::vector<double> t(np + kp, 0.0);
    double dx = 1.0 / (np + 1 - kp);

    std::vector<double> xis(np, 0.0);
    for (uint32_t i = 0; i < np; i++) {
        xis[i] = (double)i / (np - 1);
    }

    for (uint32_t i = kp; i <= np + 1; i++) {
        t[i-1] = (i - kp) * dx; // to compensate for c++ indexing we use i-1
    }

    if (pt % 2 == 0) {
        for (uint32_t ii = pt + 2; ii <= np; ii++) {
            uint32_t id1 = ii - pt/2 - 1;
            uint32_t id2 = ii - pt/2;
            t[ii-1] = (xis[id1-1] + xis[id2-1]) / 2.0; // compensate for c++ index
        }
    } else {
        for (uint32_t ii = pt + 2; ii <= np; ii++) {
            uint32_t id1 = ii - (pt + 1) / 2;
            t[ii-1] = xis[id1-1];
        }
    }

    for (uint32_t i = 1; i <= kp; i++) {
        t[i-1] = t[kp-1] - 10 * EPS;
    }

    for (uint32_t i = np + 1; i <= np + kp; i++) {
        t[i-1] = t[np] + 10 * EPS;
    }

    return t;
}


std::vector<double> DFDM::utils::get_knot_vector(uint32_t n, uint32_t k){
    std::vector<double> t(n + k, 0.0);
    double dx = 1.0 / (n + 1 - k);

    for (uint32_t i = k; i <= n + 1; i++) {
        t[i-1] = (i - k) * dx; // to compensate for c++ indexing we use i-1
    }

    for (uint32_t i = 1; i <= k; i++) {
        t[i-1] = t[k-1] - 20 * EPS;
    }

    for (uint32_t i = n + 1; i <= n + k; i++) {
        t[i-1] = t[n] + 20 * EPS;
    }

    return t;
}


double DFDM::utils::bspline_c(const std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo){

    double b = 0;
    double c1 = 0 ,c2 = 0;
    if (k == 1){
        if (xo > t_[i] && !(xo > t_[i+1])){
            b = 1;
        }else{
            b = 0;
        }
    }else{
        if( std::abs(t_[i+k-1] - t_[i]) > EPS){
            c1 = (xo - t_[i])/(t_[i+k-1] - t_[i]);
        }else{
            c1 = 0;
        }

        if(t_[i+k] - t_[i+1] > EPS){
            c2 = (t_[i+k]-xo)/(t_[i+k]-t_[i+1]);
        }else{
            c2 = 0;
        }
        b = c1*bspline_c(t_,n,i,k-1,xo)+c2*bspline_c(t_,n,i+1,k-1,xo);
    }
    return b;
}

std::pair<std::vector<double>, std::vector<double>> DFDM::utils::lgwt(int N, double a, double b) {
    N = N - 1;
    int N1 = N + 1;
    int N2 = N + 2;
    std::vector<double> xu(N1);
    for (int i = 0; i < N1; i++) {
        xu[i] = -1.0 + (2.0 * i) / (N1);
    }

    std::vector<double> y(N1);
    for (int i = 0; i < N1; i++) {
        y[i] = std::cos((2 * (i + 1) - 1) * M_PI / (2 * N + 2)) + (0.27 / N1) * std::sin(M_PI * xu[i] * N / N2);
    }

    std::vector<std::vector<double>> L(N1, std::vector<double>(N2, 0.0));
    std::vector<std::vector<double>> Lp(N1, std::vector<double>(N2, 0.0));

    double y0 = 2.0;
    while (true) {
        double max_diff = 0.0;
        for (int i = 0; i < N1; i++) {
            max_diff = std::max(max_diff, std::abs(y[i] - y0));
        }
        if (max_diff < EPS) break;

        for (int i = 0; i < N1; i++) {
            L[i][0] = 1.0;
            L[i][1] = y[i];

            Lp[i][0] = 0.0;
            Lp[i][1] = 1.0;
        }

        for (int k = 2; k <= N1; k++) {
            for (int i = 0; i < N1; i++) {
                L[i][k] = ((2 * k - 1) * y[i] * L[i][k - 1] - (k - 1) * L[i][k - 2]) / k;
            }
        }

        for (int i = 0; i < N1; i++) {
            for (int k = 0; k < N2; k++) {
                Lp[i][k] = (N2 * (L[i][N1-1] - y[i] * L[i][N2-1])) / (1 - y[i] * y[i]);
            }
        }

        for (int i = 0; i < y.size(); i++) {
            y0 = y[i];
            y[i] = y0 - L[i][N2 - 1] / Lp[i][N1 - 1];
        }
    }

    std::vector<double> x(y.size());
    for (int i = 0; i < y.size(); i++) {
        x[i] = (a * (1 - y[i]) + b * (1 + y[i])) / 2.0;
    }

    std::vector<double> w(N1);
    for (int i = 0; i < y.size(); i++) {
        w[i] = (b - a) / ((1 - y[i] * y[i]) * Lp[i][N1 - 1] * Lp[i][N1 - 1]) * (N2 / N1) * (N2 / N1);
    }

    std::reverse(x.begin(), x.end());
    std::reverse(w.begin(), w.end());

    return std::make_pair(x, w);
}


DFDM::matrix<double> DFDM::utils::dFdxp(const DFDM::matrix<double>& KK, const DFDM::matrix<double>& U, const DFDM::matrix<double>& invL_t, const DFDM::matrix<double>& invL, const DFDM::matrix<double>& Umo, const DFDM::matrix<double>& Upo){
    // to obtain dUdx111 with U211
    // from orth base to bspline
    auto Ub = invL_t.matprod(U);
    auto dUdxp = (KK* (-1)).matprod(Ub); // derivative 
    // impose the boundary 
    auto row_mo = dUdxp.row_get(0); // get the first row
    std::transform(row_mo.begin(), row_mo.end(), Umo.data_vec.begin(), row_mo.begin(), [](double a, double b){ return a - b; }); // subtract Umo from the first row
    dUdxp.row_set(0, row_mo); // set the first row

    auto row_po = dUdxp.row_get(dUdxp.rows - 1); // get the last row
    std::transform(row_po.begin(), row_po.end(), Upo.data_vec.begin(), row_po.begin(), [](double a, double b){ return a + b; }); // add Upo to the last row
    dUdxp.row_set(dUdxp.rows - 1, row_po); // set the last row
    // back to orth base
    dUdxp = invL.matprod(dUdxp); // this can be replaced with fwd bwd method.

    return dUdxp;
}


DFDM::matrix<double> DFDM::utils::dFdzp(const DFDM::matrix<double>& KK, const DFDM::matrix<double>& U, const DFDM::matrix<double>& invL_t, const DFDM::matrix<double>& invL, const DFDM::matrix<double>& Uom, const DFDM::matrix<double>& Uop){
    // to obtain dUdzp11 with U12
    auto UT = U.transpose_inplace();
    auto Ub = invL_t.matprod(UT);
    auto dUdzp = (KK*-1).matprod(Ub); // derivative 
    dUdzp = dUdzp.transpose_inplace();
    // impose the boundary 
    auto col_om = dUdzp.col_get(0); // get the first column
    std::transform(col_om.begin(), col_om.end(), Uom.data_vec.begin(), col_om.begin(), [](double a, double b){ return a - b; }); // subtract U12om from the first column
    dUdzp.col_set(0, col_om); // set the first column
    
    auto col_op = dUdzp.col_get(dUdzp.cols - 1); // get the last column
    std::transform(col_op.begin(), col_op.end(), Uop.data_vec.begin(), col_op.begin(), [](double a, double b){ return a + b; }); // add U12op to the last column
    dUdzp.col_set(dUdzp.cols - 1, col_op); // set the last column
    //
    dUdzp = dUdzp.transpose_inplace(); // permute(dUdzp,[2,1])
    dUdzp = invL.matprod(dUdzp); // pagemtimes(invLz11,dUdzp11)
    dUdzp = dUdzp.transpose_inplace(); // permute(dUdzp11,[2,1])

    return dUdzp;
}


DFDM::matrix<double> DFDM::utils::tensorProduct2D(const DFDM::matrix<double>& Bx_mat, 
                                                  const DFDM::matrix<double>& By_mat, 
                                                  const DFDM::matrix<double>& coef2D){
    auto val2D = Bx_mat.matprod(coef2D);
    val2D = val2D.transpose_inplace();
    val2D = By_mat.matprod(val2D);
    val2D = val2D.transpose_inplace();

    return val2D;
}

double DFDM::utils::bilinear_interpolate(const DFDM::matrix<double>& grid2d11, uint64_t Nx1, uint64_t Nz1, double xlocal, double zlocal) {
    double index_x_db = (double)(Nx1 - 1) * xlocal;
    double index_z_db = (double)(Nz1 - 1) * zlocal;

    uint64_t index_x_floor = static_cast<uint64_t>(std::floor(index_x_db));
    uint64_t index_z_floor = static_cast<uint64_t>(std::floor(index_z_db));
    uint64_t index_x_ceil = static_cast<uint64_t>(std::ceil(index_x_db));
    uint64_t index_z_ceil = static_cast<uint64_t>(std::ceil(index_z_db));

    if (index_x_ceil >= Nx1) index_x_ceil = Nx1 - 1;
    if (index_z_ceil >= Nz1) index_z_ceil = Nz1 - 1;
    if (index_x_floor >= Nx1) index_x_floor = Nx1 - 1;
    if (index_z_floor >= Nz1) index_z_floor = Nz1 - 1;

    if (index_x_ceil < 0) index_x_ceil = 0;
    if (index_z_ceil < 0) index_z_ceil = 0;
    if (index_x_floor < 0) index_x_floor = 0;
    if (index_z_floor < 0) index_z_floor = 0;

    double remainder_x = index_x_db - index_x_floor;
    double remainder_z = index_z_db - index_z_floor;

    return grid2d11(index_x_floor, index_z_floor) * (1.0 - remainder_x) * (1.0 - remainder_z) +
           grid2d11(index_x_ceil, index_z_floor) * (remainder_x) * (1.0 - remainder_z) +
           grid2d11(index_x_floor, index_z_ceil) * (1.0 - remainder_x) * (remainder_z) +
           grid2d11(index_x_ceil, index_z_ceil) * (remainder_x) * (remainder_z);
}

double DFDM::utils::interpolate(std::vector<double> param, std::vector<double> r_param, int num_layers, double r, bool lower_edge, bool upper_edge){
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