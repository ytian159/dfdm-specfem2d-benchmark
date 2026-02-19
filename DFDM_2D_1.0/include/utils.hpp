#pragma once
#include <vector>
#include "matrix.hpp"
#include <cmath>
#include <algorithm>

#define EPS 2.2204e-16


namespace DFDM{
    namespace utils{
      struct double_int{
      double val;
      int rank;
    };

      double max2d(DFDM::matrix<double>& in_);
      std::vector<double> get_knot_vector_shape(uint32_t np, uint32_t kp);
      std::vector<double> get_knot_vector(uint32_t n, uint32_t k);
      double dbspline_c(const std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo, uint32_t d_order);
      double bspline_c(const std::vector<double>& t_, uint32_t n, uint32_t i, uint32_t k, double xo);
      double hat_pts(DFDM::matrix<double>& hatpoints1, DFDM::matrix<double>& hatpoints2, uint32_t N, uint32_t p);
      std::pair<std::vector<double>, std::vector<double>> lgwt(int N, double a, double b);
      DFDM::matrix<double> dFdxp(const DFDM::matrix<double>& KK, const DFDM::matrix<double>& U, 
                                const DFDM::matrix<double>& invL_t, const DFDM::matrix<double>& invL, 
                                const DFDM::matrix<double>& Umo, const DFDM::matrix<double>& Upo);
      DFDM::matrix<double> dFdzp(const DFDM::matrix<double>& KK, const DFDM::matrix<double>& U, 
                                 const DFDM::matrix<double>& invL_t, const DFDM::matrix<double>& invL, 
                                 const DFDM::matrix<double>& Uom, const DFDM::matrix<double>& Uop);
      double bilinear_interpolate(const DFDM::matrix<double>& grid2d11, uint64_t Nx1, uint64_t Nz1, double xlocal, double zlocal);
      DFDM::matrix<double> tensorProduct2D(const DFDM::matrix<double>& Bx_mat, 
                                                  const DFDM::matrix<double>& By_mat, 
                                                  const DFDM::matrix<double>& coef2D);
      double interpolate(std::vector<double> param, std::vector<double> r_param, int num_layers, double r, bool lower_edge, bool upper_edge);
    }
}

