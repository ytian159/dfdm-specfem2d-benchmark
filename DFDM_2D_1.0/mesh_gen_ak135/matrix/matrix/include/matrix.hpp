#pragma once

#include <iostream>
#include <vector>
#include <fstream>

#include "cblas.h"
#include "lapacke.h"
#include <algorithm>
#include <numeric>

namespace DFDM{
    template <class Type>
    class matrix{
        public:
            std::vector<Type> data_vec;
            uint32_t rows, cols;
            matrix();
            matrix(const matrix& mc);
            matrix(uint32_t rows_, uint32_t cols);
            void resize(uint32_t rows_, uint32_t cols);
            Type operator()(uint64_t row, uint64_t col) const;
            matrix& operator=(const matrix& matE);
            template<class D>
            matrix operator*(const D scalar)const;
            matrix operator+(const matrix& inm);
            matrix operator-(const matrix& inm);
            void value_set(uint32_t row, uint32_t col, Type in);
            void transpose(matrix& matT);
            matrix transpose_inplace()const;
            void reshape(uint32_t rows_, uint32_t cols_);
            void  inverse(matrix& matI);
            void cholesky(matrix& matl, matrix& matlt);
            void sqrtm_schur(matrix& matl, matrix& matlt);
            matrix vecprod(const matrix& vec);
            matrix matprod(const matrix& matI) const;
            void print();
            void print_file(std::string file_name) const;
            static void row_copy(const matrix& mat_src, uint32_t src_row_id, matrix& mat_tar, uint32_t tar_row_id);
            static void col_copy(const matrix& mat_src, uint32_t src_col_id, matrix& mat_tar, uint32_t tar_col_id);
            static matrix elementwise_mult(const matrix& ina, const matrix& inb);
            static matrix elementwise_add(const matrix& ina, const matrix& inb);
            static matrix elementwise_div(const matrix& ina, const matrix& inb);
            void row_set(uint32_t tgt_row, std::vector<Type> row_data);
            std::vector<Type> row_get(uint32_t row_id);
            void col_set(const uint32_t tgt_col, const std::vector<Type>& col_data);
            std::vector<Type> col_get(uint32_t col_id);
    };
}