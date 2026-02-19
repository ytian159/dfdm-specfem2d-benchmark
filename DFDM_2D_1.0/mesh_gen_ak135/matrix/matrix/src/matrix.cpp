#include "matrix.hpp"
template class DFDM::matrix<double>;
template DFDM::matrix<double> DFDM::matrix<double>::operator*(const double scalar)const;
template DFDM::matrix<double> DFDM::matrix<double>::operator*(const int32_t scalar)const;

template<class T>
DFDM::matrix<T>::matrix():rows(0), cols(0){
    data_vec.resize(rows*cols);
}

template<class T>
DFDM::matrix<T>::matrix(const DFDM::matrix<T>& mc){
    rows = mc.rows;
    cols = mc.cols;
    data_vec = mc.data_vec;
}

// initializes the matrix to zero
template<class T>
DFDM::matrix<T>::matrix(uint32_t rows_, uint32_t cols_): rows(rows_), cols(cols_){
    data_vec.resize(rows*cols,0);
}

//also sets everything to zero
template<class T>
void DFDM::matrix<T>::resize(uint32_t rows_, uint32_t cols_){
    rows = rows_;
    cols = cols_;
    data_vec.resize(rows*cols,0);
}
template<class T>
T DFDM::matrix<T>::operator()(uint64_t row, uint64_t col) const {
    if(row >= rows || col >= cols){
        std::cerr << "ERROR: Out of bounds access in Matrix operator(). Row index: " << row << ", Column index: " << col << std::endl;
        exit(-1);
    }
    return data_vec[row*cols + col];
}

template<class T>
DFDM::matrix<T>& DFDM::matrix<T>::operator=(const DFDM::matrix<T>& matE){
    rows = matE.rows;
    cols = matE.cols;
    data_vec = matE.data_vec;
    return *this;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_mult(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){
    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match." << std::endl;
        exit(1);
    }
    DFDM::matrix<T> prod(ina.rows, ina.cols);
    std::vector<uint32_t> indx(ina.rows*ina.cols);
    std::iota(indx.begin(), indx.end(), 1);
    std::for_each(std::begin(indx), std::end(indx), [&](uint32_t idx){
        prod.data_vec[idx - 1] = ina.data_vec[idx - 1] * inb.data_vec[idx - 1];
    });
    
    return prod;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_add(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){
    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match." << std::endl;
        exit(1);
    }
    DFDM::matrix<T> prod(ina.rows, ina.cols);
    std::vector<uint32_t> indx(ina.rows*ina.cols);
    std::iota(indx.begin(), indx.end(), 1);
    std::for_each(std::begin(indx), std::end(indx), [&](uint32_t idx){
        prod.data_vec[idx - 1] = ina.data_vec[idx - 1] + inb.data_vec[idx - 1];
    });
    
    return prod;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_div(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){
    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match in elementwise_div." << std::endl;
        exit(1);
    }
    DFDM::matrix<T> prod(ina.rows, ina.cols);
    std::vector<uint32_t> indx(ina.rows*ina.cols);
    std::iota(indx.begin(), indx.end(), 1);
    std::for_each(std::begin(indx), std::end(indx), [&](uint32_t idx){
        prod.data_vec[idx - 1] = ina.data_vec[idx - 1] / inb.data_vec[idx - 1];
    });
    
    return prod;
}

template<class T>
void DFDM::matrix<T>::value_set(uint32_t row, uint32_t col, T in){
    if(row >= rows || col >= cols){
        std::cerr << "ERROR: Out of bounds access in Matrix value_set. Row index: " << row << ", Column index: " << col 
               << ".Actual matrix dimensions: rows = " << rows << ", cols = " << cols << std::endl; 
        exit(-1);
    }
    data_vec[row*cols + col] = in;
}

template<class T>
void DFDM::matrix<T>::transpose(DFDM::matrix<T>& matT){
    matT.resize(this->cols, this->rows);
    cblas_domatcopy(CblasRowMajor, CblasTrans, this->rows , this->cols, 1.0, this->data_vec.data(), this->cols, matT.data_vec.data(), matT.cols);
}

//does not modify the original matrix, copies result in new matrix
template<class T>
DFDM::matrix<T> DFDM::matrix<T>::transpose_inplace()const{
    DFDM::matrix<T> matT(this->cols, this->rows);
    cblas_domatcopy(CblasRowMajor, CblasTrans, this->rows , this->cols, 1.0, this->data_vec.data(), this->cols, matT.data_vec.data(), matT.cols);
    return matT;
}

template<class T>
void DFDM::matrix<T>::reshape(uint32_t rows_, uint32_t cols_){
    if(rows_ * cols_ != rows * cols){
        std::cerr << "ERROR: Total number of elements must remain the same after reshape." << std::endl;
        exit(-1);
    }

    // std::vector<T> temp_data_vec = data_vec;
    // for (uint32_t i = 0; i < rows * cols; i++) {
    //     uint32_t new_row = i / cols_;
    //     uint32_t new_col = i % cols_;
    //     data_vec[new_row * cols_ + new_col] = temp_data_vec[i];
    // }

    rows = rows_;
    cols = cols_;
}



template<class T>
void DFDM::matrix<T>::inverse(DFDM::matrix<T>& matI){
    matI.rows = this->rows;
    matI.cols = this->cols;
    matI.data_vec = this->data_vec;
 
    std::vector<int32_t> ipiv(std::min(matI.rows,matI.cols));
    auto inf =  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, matI.rows, matI.cols, matI.data_vec.data(), matI.cols, ipiv.data());
    if (inf > 0){
        std::cerr<<"ERROR: error computing matrix inverse, factor of U is singular!!"<<"\n";
        exit(-1);
    }else if(inf !=0){
        std::cerr<<"ERROR: error computing matrix inverse, illegal arguement: "<<inf<<"\n";
        exit(-1);      
    }else{
        inf = LAPACKE_dgetri(LAPACK_ROW_MAJOR, matI.rows, matI.data_vec.data(), matI.cols, ipiv.data());
    }
}

template<class T>
void DFDM::matrix<T>::cholesky(DFDM::matrix<T>& matl, DFDM::matrix<T>& matlt){
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform Cholesky decomposition." << std::endl;
        exit(1);
    }
    for (int i = 0; i < this->rows; i++) {
        for(int j = 0; j < this->cols; j++){
            if(std::abs(this->data_vec[i*this->cols + j] - this->data_vec[j*this->cols + i]) > 1e-8){
                std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not symmetric, cannot perform Cholesky decomposition." << std::endl;
                exit(1);
            }
        }
    }

    matl.rows = this->rows;
    matl.cols = this->cols;
    matl.data_vec = this->data_vec;
 
    auto res = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', matl.rows, matl.data_vec.data(), matl.cols);
    if(res != 0){
        std::cerr<<"Cholesky Failed with error:" << res << std::endl;
        exit(-1);
    }

    for (int i = 0; i < matl.rows; i++) {
        for(int j = 0; j < matl.cols; j++){
            if( j > i)
                matl.data_vec[i*matl.cols + j] = 0;
        }
    }

    matl.transpose(matlt);
}

template<class T>
void DFDM::matrix<T>::sqrtm_schur(DFDM::matrix<T>& matl, DFDM::matrix<T>& matlt){
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform Schur decomposition." << std::endl;
        exit(1);
    }
    
    DFDM::matrix<double> U(this->rows,this->rows);
    DFDM::matrix<double> R(this->rows,this->rows);
    std::vector<double> wr(this->rows), wi(this->rows);

    matl.rows = this->rows;
    matl.cols = this->cols;
    matl.data_vec = this->data_vec;
    
        // Perform Schur decomposition (LAPACK function dgees)
    int sdim;
    auto info =LAPACKE_dgees(LAPACK_ROW_MAJOR, 'V', 'N', NULL, this->rows, matl.data_vec.data(), U.rows, &sdim, wr.data(), wi.data(), U.data_vec.data(), this->rows);

    if (info != 0) {
        std::cerr << "Error in Schur decomposition: " << info << std::endl;
        return;
    }
    // Compute the square root of the upper triangular matrix matl and store in R
    int n = this->rows;
    for (int i = 0; i < this->rows; ++i) {
        R.data_vec[i * n + i] = std::sqrt(matl.data_vec[i * n + i]);
        for (int j = i + 1; j < this->rows; ++j) {
            double sum = 0;
            for (int k = i; k < j; ++k) {
                sum += R.data_vec[i * n + k] * R.data_vec[k * n + j];
            }
            R.data_vec[i * n + j] = (matl.data_vec[i * n + j] - sum) / (R.data_vec[i * n + i] + R.data_vec[j * n + j]);
        }
    }
    // Transform back to the original basis
    // sqrtA = U * R * U^T

    std::vector<double> temp(n * n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, U.data_vec.data(), n, R.data_vec.data(), n, 0.0, temp.data(), n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, temp.data(), n, U.data_vec.data(), n, 0.0, matl.data_vec.data(), n);
    matl.transpose(matlt);
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::vecprod(const DFDM::matrix<T>& vec){
    DFDM::matrix<T> prod(this->rows, vec.cols);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, this->rows, this->cols, 1.0, this->data_vec.data(), this->cols, vec.data_vec.data(), 1, 0.0, prod.data_vec.data(), 1);

    return prod;
}
template<class T>
DFDM::matrix<T> DFDM::matrix<T>::matprod(const DFDM::matrix<T>& matI) const {
    DFDM::matrix<T> prod(this->rows, matI.cols);
    if(this->cols != matI.rows){
        std::cerr<<"Prod failed, mat dim mismatch"<< std::endl;
        exit(-1);
    }
    // cblas_dgemv(CblasRowMajor, CblasNoTrans, this->rows, this->cols, 1.0, this->data_vec.data(), this->cols, vec.data_vec.data(), 1, 0.0, prod.data_vec.data(), 1);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, this->rows, matI.cols, this->cols, 1.0, this->data_vec.data(), this->cols, matI.data_vec.data(), matI.cols, 0.0, prod.data_vec.data(), prod.cols);

    return prod;
}
template<class T>
template<class D>
DFDM::matrix<T> DFDM::matrix<T>::operator*(const D scalar)const{
    DFDM::matrix<T> sprod(rows,cols);
    sprod.data_vec = data_vec;
    cblas_dscal(rows*cols, scalar, sprod.data_vec.data(), 1);
    return sprod;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::operator+(const matrix& inm){
    if(inm.cols != cols || inm.rows != rows){
        std::cerr << "Dims Don't Match for Matrix/Vector Sum"<<std::endl;
        exit(-1);
    }
    DFDM::matrix<T> ssum(rows,cols);
    ssum.data_vec = data_vec;
    cblas_daxpy(rows*cols, 1.0, inm.data_vec.data(), 1,ssum.data_vec.data() , 1);
    return ssum;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::operator-(const matrix& inm){
    if(inm.cols != cols || inm.rows != rows){
        std::cerr << "Dims Don't Match for Matrix/Vector Sub"<<std::endl;
        exit(-1);
    }
    DFDM::matrix<T> ssum(rows,cols);
    ssum.data_vec = data_vec;
    cblas_daxpy(rows*cols, -1.0, inm.data_vec.data(), 1,ssum.data_vec.data() , 1);
    return ssum;
}

template<class T>
void DFDM::matrix<T>::row_copy(const matrix& mat_src, uint32_t src_row_id, matrix& mat_tar, uint32_t tar_row_id){
    if(mat_src.cols != mat_tar.cols){
        std::cerr << "Dimension mismatch in row copy" << std::endl;
        exit(-1);
    }
    
    if(src_row_id >= mat_src.rows || tar_row_id >= mat_tar.rows){
        std::cerr << " ERROR:Rows to be copied cannot be larger than total number of rows in target or source matrix!" << std::endl;
        exit(-1);
    }
    auto src_it = mat_src.data_vec.begin();
    std::advance(src_it, src_row_id * mat_src.cols);

    auto tar_it = mat_tar.data_vec.begin();
    std::advance(tar_it, tar_row_id * mat_tar.cols);

    std::copy(src_it, src_it + mat_src.cols, tar_it);
}

template<class T>
void DFDM::matrix<T>::col_copy(const matrix& mat_src, uint32_t src_col_id, matrix& mat_tar, uint32_t tar_col_id){
    if(mat_src.rows != mat_tar.rows){
        std::cerr << "Dimension mismatch in col copy" << std::endl;
        exit(-1);
    }
    if(src_col_id >= mat_src.cols || tar_col_id >= mat_tar.cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }

    auto src_it = mat_src.data_vec.begin();
    std::advance(src_it, src_col_id);

    auto tar_it = mat_tar.data_vec.begin();
    std::advance(tar_it, tar_col_id);
    if(mat_src.cols == 1 && mat_tar.cols == 1){ // simple vector copy
        std::copy(src_it, src_it + mat_src.rows, tar_it);
    }else{
        for(uint32_t i = 0; i < mat_src.rows; i++){
            *tar_it = *src_it;
            std::advance(src_it, mat_src.cols);
            std::advance(tar_it, mat_tar.cols);
        }
    }
}


template<class T>
void DFDM::matrix<T>::row_set(uint32_t tgt_row, std::vector<T> row_data){
    if(row_data.size() != cols){
        std::cerr << "Dimension mismatch in row set, size of source matrix is not equal to target matrix" << std::endl;
        exit(-1);
    }else if(tgt_row >= rows){
        std::cerr << " ERROR:Rows to be copied cannot be larger than total number of rows in target or source matrix!" << std::endl;
        exit(-1);
    }

    std::copy(row_data.begin(), row_data.end(), data_vec.begin() + tgt_row * cols);
}

template<class T>
void DFDM::matrix<T>::col_set(const uint32_t tgt_col, const std::vector<T>& col_data){
    if(col_data.size() != rows){
        std::cerr << "Dimension mismatch in col_set, size of source matrix is not equal to target matrix" << ", matrix_rows:" << rows << " col length:"<< col_data.size()<< std::endl;
        exit(-1);
    }
    if(tgt_col >= cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }

    for (uint32_t i = 0; i < rows; i++) {
        data_vec[i * cols + tgt_col] = col_data[i];
    }
}

template<class T>
std::vector<T> DFDM::matrix<T>::col_get(uint32_t col_id){
    if(col_id >= cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }
    std::vector<T> col_data(rows);
    for (uint32_t i = 0; i < rows; i++) {
        col_data[i] = data_vec[i * cols + col_id];
    }
    return col_data;
}

template<class T>
std::vector<T> DFDM::matrix<T>::row_get(uint32_t row_id){
    if(row_id >= rows){
        std::cerr << " ERROR:Rows to be copied cannot be larger than total number of rows in target or source matrix!" << std::endl;
        exit(-1);
    }
    std::vector<T> row_vec(cols);
    std::copy(data_vec.begin() + row_id * cols, data_vec.begin() + (row_id + 1) * cols, row_vec.begin());
    return row_vec;
}

template<class T>
void DFDM::matrix<T>::print(){
    for(uint32_t i = 0; i < rows; i++){
        for(uint32_t j = 0; j < cols; j++){
            std::cout<< data_vec[i*cols + j] <<"\t";
        }
        std::cout << "\n";
    }
}

template<class T>
void DFDM::matrix<T>::print_file(std::string file_name) const{
    std::ofstream file_out(file_name);
    // file_out << "rows:" << rows << "\tcols:" << cols << std::endl;
    for(uint32_t i = 0; i < rows; i++){
        for(uint32_t j = 0; j < cols; j++){
            file_out<< data_vec[i*cols + j];
            if(j != cols - 1){
                file_out<<",";
            }
        }
        file_out << "\n";
    }
    file_out.close();
}