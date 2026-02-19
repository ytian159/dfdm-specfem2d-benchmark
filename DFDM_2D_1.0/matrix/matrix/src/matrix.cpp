#include "matrix.hpp"
#include <type_traits>
template class DFDM::matrix<double>;
template DFDM::matrix<double> DFDM::matrix<double>::operator*(const double scalar)const;
template DFDM::matrix<double> DFDM::matrix<double>::operator*(const int32_t scalar)const;

template<class T>
DFDM::matrix<T>::matrix():rows(0), cols(0), banded(false), kd(-1){
    data_vec.resize(rows*cols);
}

template<class T>
DFDM::matrix<T>::matrix(const DFDM::matrix<T>& mc){
    rows = mc.rows;
    cols = mc.cols;
    banded = mc.banded;
    data_vec = mc.data_vec;
    kd = mc.kd;
}

// initializes the matrix to zero
template<class T>
DFDM::matrix<T>::matrix(uint32_t rows_, uint32_t cols_): rows(rows_), cols(cols_){
    data_vec.resize(rows*cols,0);
    banded = false;
}

//also sets everything to zero
template<class T>
void DFDM::matrix<T>::resize(uint32_t rows_, uint32_t cols_){
    if(banded == true){
        std::cerr << "ERROR: Cannot resize a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }
    rows = rows_;
    cols = cols_;
    data_vec.resize(rows*cols,0);
    banded = false;
}

template<class T>
T DFDM::matrix<T>::operator()(uint64_t row, uint64_t col) const {
    if (banded) {
        std::cerr << "ERROR: Attempting to access elements of a banded matrix using operator(). Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }
    
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
    banded = matE.banded;
    kd = matE.kd;
    return *this;
}


template<class T>
void DFDM::matrix<T>::triangular_solve(const matrix& matL, matrix& matX, matrix& matB, const char* uplo){
    if(matL.banded == true){
        std::cerr << "ERROR: Cannot perform triangular solve on a banded matrix. Please convert it back to a full matrix first or use banded_triangular_solve." << std::endl;
        exit(-1);
    }

    if(matL.rows != matL.cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform triangular solve." << std::endl;
        exit(1);
    }
    if(matL.rows != matB.rows ){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix dimensions do not match for triangular solve." << std::endl;
        exit(1);
    }
    // Compare the first character pointed to by uplo
    if((*uplo != 'L') && (*uplo != 'U')){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Invalid value for uplo parameter. Use 'L' for lower triangular or 'U' for upper triangular. Current value:" << uplo << std::endl;
        exit(1);
    }

    if(*uplo == 'U')
        cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, matB.rows, matB.cols, 1.0, matL.data_vec.data(), matL.cols, matB.data_vec.data(), matB.cols);
    else
        cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, matB.rows, matB.cols, 1.0, matL.data_vec.data(), matL.cols, matB.data_vec.data(), matB.cols);

    matX = matB;
}

template<class T>
void DFDM::matrix<T>::banded_triangular_solve(const matrix& matL, matrix& matX, matrix& matB, const char* uplo){


    if(matL.banded == false){
        std::cerr << "ERROR: Banded triangular solve requires the matrix to be in banded format. Please convert it to banded format first." << std::endl;
        exit(-1);
    }

    if(matL.rows != matL.cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform triangular solve." << std::endl;
        exit(1);
    }
    if(matL.rows != matB.rows){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix dimensions do not match for triangular solve." << std::endl;
        exit(1);
    }
    // Compare the first character pointed to by uplo
    if((*uplo != 'L') && (*uplo != 'U')){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Invalid value for uplo parameter. Use 'L' for lower triangular or 'U' for upper triangular. Current value:" << uplo << std::endl;
        exit(1);
    }


    /*
(	int 	matrix_layout,
char 	uplo,
char 	trans,
char 	diag,
lapack_int 	n,
lapack_int 	kd,
lapack_int 	nrhs,
const double * 	ab,
lapack_int 	ldab,
double * 	b,
lapack_int 	ldb 
)	
    */
    if(*uplo == 'U')
        LAPACKE_dtbtrs(LAPACK_ROW_MAJOR, *uplo, 'N', 'N', matL.rows, matL.kd, matB.cols, matL.data_vec.data() , matL.cols, matB.data_vec.data(), matB.cols);
    else
        LAPACKE_dtbtrs(LAPACK_ROW_MAJOR, *uplo, 'N', 'N', matL.rows, matL.kd, matB.cols, matL.data_vec.data() , matL.cols, matB.data_vec.data(), matB.cols);

    matX = matB;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::toBanded(int kd, const char* uplo){

    if(banded == true){
        std::cerr << "ERROR: Matrix is already in banded format." << std::endl;
        exit(-1);
    }

    if(rows != cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot convert to banded format." << std::endl;
        exit(1);
    }

    int rows_banded = kd + 1;
    DFDM::matrix<T> banded_mat(rows, rows);
    banded_mat.data_vec.resize(rows_banded * rows, 0.0); // Initialize with zeros

    // kd: number of sub/super diagonals
    // uplo: 'U' for upper, 'L' for lower triangular
    int n = rows;
    if (*uplo == 'U') {
        // Upper banded: AB(kd+1+i-j, j) = A(i, j) for max(0, j-kd) <= i <= j (from BLAS docs)
        for (int j = 0; j < n; ++j) {
            int i_start = std::max(0, j - kd);
            for (int i = i_start; i <= j; ++i) {
                int ab_row = kd + i - j;
                banded_mat.data_vec[ab_row * n + j] = this->data_vec[i * n + j];
            }
        }
    } else if (*uplo == 'L') {
        // Lower banded: AB(i-j, j) = A(i, j) for j <= i <= min(n-1, j+kd) (from BLAS docs)
        for (int j = 0; j < n; ++j) {
            int i_end = std::min(n - 1, j + kd);
            for (int i = j; i <= i_end; ++i) {
                int ab_row = i - j;
                banded_mat.data_vec[ab_row * n + j] = this->data_vec[i * n + j];
            }
        }
    } else {
        std::cerr << "Invalid value for uplo in toBanded: " << *uplo << std::endl;
        exit(1);
    }

    banded_mat.banded = true;
    banded_mat.kd = kd;

    return banded_mat;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_mult(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){
    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match." << std::endl;
        exit(1);
    }

    DFDM::matrix<T> prod(ina.rows, ina.cols);
    
    // Use size_t to prevent overflow on very large matrices (>4 billion elements)
    const size_t total_elements = static_cast<size_t>(ina.rows) * ina.cols;

    // OpenMP directive to parallelize the loop
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < total_elements; ++i) {
        prod.data_vec[i] = ina.data_vec[i] * inb.data_vec[i];
    }
    
    return prod;
}
// template<class T>
// DFDM::matrix<T> DFDM::matrix<T>::elementwise_mult(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){

//     if(ina.banded == true || inb.banded == true){
//         std::cerr << "Error: Element-wise multiplication not supported for banded matrices." << std::endl;
//         exit(1);
//     }

//     if(ina.rows != inb.rows || ina.cols != inb.cols){
//         std::cerr << "Error: Matrix dimensions do not match." << std::endl;
//         exit(1);
//     }
//     DFDM::matrix<T> prod(ina.rows, ina.cols);
//     std::vector<uint32_t> indx(ina.rows*ina.cols);
//     std::iota(indx.begin(), indx.end(), 1);
//     std::for_each(std::begin(indx), std::end(indx), [&](uint32_t idx){
//         prod.data_vec[idx - 1] = ina.data_vec[idx - 1] * inb.data_vec[idx - 1];
//     });
    
//     return prod;
// }

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_add(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){

    if(ina.banded == true || inb.banded == true){
        std::cerr << "Error: Element-wise addition not supported for banded matrices." << std::endl;
        exit(1);
    }

    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match." << std::endl;
        exit(1);
    }
    DFDM::matrix<T> prod(ina.rows, ina.cols);
    const size_t total_elements = static_cast<size_t>(ina.rows) * ina.cols;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < total_elements; ++i) {
        prod.data_vec[i] = ina.data_vec[i] + inb.data_vec[i];
    }
    
    return prod;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::elementwise_div(const DFDM::matrix<T>& ina, const DFDM::matrix<T>& inb){

    if(ina.banded == true || inb.banded == true){
        std::cerr << "Error: Element-wise division not supported for banded matrices." << std::endl;
        exit(1);
    }

    if(ina.rows != inb.rows || ina.cols != inb.cols){
        std::cerr << "Error: Matrix dimensions do not match in elementwise_div." << std::endl;
        exit(1);
    }
    DFDM::matrix<T> prod(ina.rows, ina.cols);
    const size_t total_elements = static_cast<size_t>(ina.rows) * ina.cols;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < total_elements; ++i) {
        prod.data_vec[i] = ina.data_vec[i] / inb.data_vec[i];
    }
    
    return prod;
}

template<class T>
void DFDM::matrix<T>::value_set(uint32_t row, uint32_t col, T in){
    if(banded == true){
        std::cerr << "ERROR: Attempting to set elements of a banded matrix using value_set. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(row >= rows || col >= cols){
        std::cerr << "ERROR: Out of bounds access in Matrix value_set. Row index: " << row << ", Column index: " << col 
               << ".Actual matrix dimensions: rows = " << rows << ", cols = " << cols << std::endl; 
        exit(-1);
    }
    data_vec[row*cols + col] = in;
}

template<class T>
void DFDM::matrix<T>::transpose(DFDM::matrix<T>& matT){
    if(banded == true){
        std::cerr << "ERROR: Cannot transpose a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }
    matT.resize(this->cols, this->rows);
    cblas_domatcopy(CblasRowMajor, CblasTrans, this->rows , this->cols, 1.0, this->data_vec.data(), this->cols, matT.data_vec.data(), matT.cols);
}

//does not modify the original matrix, copies result in new matrix
template<class T>
DFDM::matrix<T> DFDM::matrix<T>::transpose_inplace()const{
    if(banded == true){
        std::cerr << "ERROR: Cannot transpose a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    DFDM::matrix<T> matT(this->cols, this->rows);
    cblas_domatcopy(CblasRowMajor, CblasTrans, this->rows , this->cols, 1.0, this->data_vec.data(), this->cols, matT.data_vec.data(), matT.cols);
    return matT;
}

template<class T>
void DFDM::matrix<T>::reshape(uint32_t rows_, uint32_t cols_){
    if(banded == true){
        std::cerr << "ERROR: Cannot reshape a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(rows_ * cols_ != rows * cols){
        std::cerr << "ERROR: Total number of elements must remain the same after reshape." << std::endl;
        exit(-1);
    }

    rows = rows_;
    cols = cols_;
}



template<class T>
void DFDM::matrix<T>::inverse(DFDM::matrix<T>& matI){
    if(this->banded == true){
        std::cerr << "ERROR: Cannot compute inverse of a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

    if(this->banded == true){
        std::cerr << "ERROR: Cannot compute Cholesky decomposition of a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

    if(this->banded == true){
        std::cerr << "ERROR: Cannot compute square root of a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

// python code for reference
// def factorize(M_in, p):
//     M = M_in.copy()  # Avoid modifying the input matrix
//     n = M.shape[0]
//     A = np.eye(n)

//     for j in range(p, 0, -1):
//         if j % 2 == 0:
//             M = M[::-1, ::-1]  # Reverse rows and columns

//         L = np.eye(n)

//         for i in range(j, n):
//             num = M[i, i - j]
//             den = M[i - 1, i - j]  # i - 1 corresponds to i in MATLAB (1-based)
//             L[i, i - 1] = num / den if den != 0 else 0  # Avoid division by zero
//             M[i, :] = M[i, :] - L[i, i - 1] * M[i - 1, :]

//         M = M @ np.linalg.inv(L.T)

//         M[np.where(M<1e-14)]=0

//         plt.figure()
//         plt.pcolor(M,norm=LogNorm())
//         plt.title('M for now')
//         plt.colorbar()

//         if j % 2 == 0:
//             M = M[::-1, ::-1]
//             L = L[::-1, ::-1]

//         A = A @ L

//     A = A @ sqrtm(M)
//     A[np.where(A<1e-14)]=0
//     return A

template<class T>
void DFDM::matrix<T>::alt_factorisation(DFDM::matrix<T>& matl, DFDM::matrix<T>& matlt, uint64_t p){
    // Matrix factorisation recommended by Lyu et al 2024 Appendix D
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform matrix decomposition." << std::endl;
        exit(1);
    }

    int n = this->rows;
    
    // DFDM::matrix<double> M(this->rows,this->rows);
    matl.resize(n,n);
    matlt.resize(n,n);

    DFDM::matrix<double> A(n,n);
    

    DFDM::matrix<double> revert_rows_and_columns(n,n);
    for (int i = 0; i<n; i++){
        revert_rows_and_columns.value_set(i,n-i-1,1.0);
    }

    auto M = *this; // I need a deep copy

    for (int i = 0; i<n; i++){
        A.value_set(i,i,1.0);
    }

    for (int j=p; j>0; j--){
        DFDM::matrix<double> L(n,n);
        DFDM::matrix<double> LT_inv;
        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));

        }

        for (int i = 0; i<n; i++){
            L.value_set(i,i,1.0);
        }

        for (int i = j; i<n; i++){
            auto num = M(i,i-j);
            auto den = M(i-1,i-j);
            auto val_ = (den != 0) ? (num / den) : 0; // Avoid division by zero
            L.value_set(i,i-1,val_);
            for (int k = 0; k<n; k++){
                M.value_set(i,k,M(i,k) - (L(i,i-1)*M(i-1,k)));
            }
        }

        auto MT=M.transpose_inplace();

        for (int i = j; i<n; i++){
            for (int k = 0; k<n; k++){
                MT.value_set(i,k,MT(i,k) - (L(i,i-1)*MT(i-1,k)));
            }
        }

        M=MT.transpose_inplace(); // probably not necessary, since M is symetric

        // set the elements outside of the band to zero
        // might be necessary in the future, if accuracy problems arise
        for (int i = 0; i<n; i++){
            for (int k = i+j; k<n; k++){
                M.value_set(i,k,0.0);
            }

            for (int k = 0; k<i-j; k++){
                M.value_set(i,k,0.0);
            }
        }
    
        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
            L = revert_rows_and_columns.matprod(L.matprod(revert_rows_and_columns));
        }

        A = A.matprod(L);
    }
    M.sqrtm_schur(matl,matlt);
    A = A.matprod(matl);

    matl=A;
    matlt=A.transpose_inplace();
    
}


// Python code for reference
// def fact_lower_diag_twist(j, M, isplit):
//     M = M.copy()  # avoid modifying the original M
//     n = M.shape[0]
    
//     L = np.eye(n)

//     # First loop
//     for i in range(j - 1, isplit + j - 2):  # MATLAB is 1-based, Python is 0-based
//         factor = M[i + 1, i - j + 1] / M[i, i - j + 1]
//         L[i + 1, :] = L[i + 1, :] - L[i, :] * factor
//         M[i + 1, :] = M[i + 1, :] - M[i, :] * factor

//     # Reverse M and L vertically and horizontally
//     M = M[::-1, ::-1]
//     L = L[::-1, ::-1]

//     # Second loop
//     for i in range(j - 1, n - isplit):
//         factor = M[i + 1, i - j + 1] / M[i, i - j + 1]
//         L[i + 1, :] = L[i + 1, :] - L[i, :] * factor
//         M[i + 1, :] = M[i + 1, :] - M[i, :] * factor

//     # Reverse M and L back
//     M = M[::-1, ::-1]
//     L = L[::-1, ::-1]

//     # Invert L
//     L = np.linalg.inv(L)

//     return L

template<class T>
void DFDM::matrix<T>::fact_lower_diag_twist(int j, DFDM::matrix<T>& L, int isplit){
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform Schur decomposition." << std::endl;
        exit(1);
    }

    auto M = *this; // I need a deep copy

    int n = this->rows;
    
    DFDM::matrix<double> L_inv(n,n);

    DFDM::matrix<double> revert_rows_and_columns(n,n);
    for (int i = 0; i<n; i++){
        revert_rows_and_columns.value_set(i,n-i-1,1.0);
    }

    for (int i = 0; i<n; i++){
        L_inv.value_set(i,i,1.0);
    }

    for (int i=j-1; i<isplit + j - 2; i++){
        double factor = M(i+1,i-j+1) / M(i,i-j+1);
        for (int k = 0; k<n; k++){
            L_inv.value_set(i+1,k,L_inv(i+1,k) - L_inv(i,k) * factor);
            M.value_set(i+1,k,M(i+1,k) - M(i,k) * factor);
        }
    }
    
    M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
    L_inv = revert_rows_and_columns.matprod(L_inv.matprod(revert_rows_and_columns));

    for (int i=j-1; i<n - isplit; i++){
        double factor = M(i+1,i-j+1) / M(i,i-j+1);

        for (int k = 0; k<n; k++){
            L_inv.value_set(i+1,k,L_inv(i+1,k) - L_inv(i,k) * factor);
            M.value_set(i+1,k,M(i+1,k) - M(i,k) * factor);
        }
    }

    M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
    L_inv = revert_rows_and_columns.matprod(L_inv.matprod(revert_rows_and_columns));

    L_inv.inverse(L);
    
}

// Python code for reference
// def fact_mass_twist(M_in, p):
//     M = M_in.copy()  # avoid modifying input
//     n = M.shape[0]
//     A = np.eye(n)

//     for j in range(p, 0, -1):  # MATLAB p:-1:1

//         if j % 2 == 0:
//             M = M[::-1, ::-1]

//         itwist = min(n - p, n // 2 + 1)

//         L = fact_lower_diag_twist(j, M, itwist)

//         # Apply similarity transformation
//         Linv = np.linalg.inv(L)
//         M = Linv @ M @ Linv.T

//         if j % 2 == 0:
//             M = M[::-1, ::-1]
//             L = L[::-1, ::-1]

//         A = A @ L

//     # Replace diagonal with sqrt(abs(diagonal))
//     for i in range(n):
//         M[i, i] = np.sqrt(abs(M[i, i]))

//     A = A @ M

//     return A

template<class T>
void DFDM::matrix<T>::fact_mass_twist(DFDM::matrix<T>& matl, DFDM::matrix<T>& matlt, uint64_t p){
    // Other possible matrix factorisation, seems unstable
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform Schur decomposition." << std::endl;
        exit(1);
    }

    auto M = *this; // I need a deep copy
    DFDM::matrix<T> A(this->rows,this->cols);
    int n = this->rows;
    
    // DFDM::matrix<double> M(this->rows,this->rows);
    matl.resize(this->rows,this->cols);
    matlt.resize(this->rows,this->cols);

    DFDM::matrix<double> revert_rows_and_columns(this->rows,this->cols);
    for (int i = 0; i<n; i++){
        revert_rows_and_columns.value_set(i,n-i-1,1.0);
    }

    for (int i = 0; i<n; i++){
        A.value_set(i,i,1.0);
    }

    for (int j=p; j>0; j--){
        DFDM::matrix<T> L(n,n);
        DFDM::matrix<T> L_inv;
        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
        }

        int itwist = std::min((int)(n - p), n / 2 + 1);

        M.fact_lower_diag_twist(j, L, itwist);
        L.inverse(L_inv);

        M = L_inv.matprod(M.matprod(L_inv.transpose_inplace()));

        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
            L = revert_rows_and_columns.matprod(L.matprod(revert_rows_and_columns));
        }

        A = A.matprod(L);
    }

    for (int i = 0; i<n; i++){
        M.value_set(i,i,std::sqrt(std::abs(M(i,i))));
    }
    A = A.matprod(M);

    matl=A;
    matlt=A.transpose_inplace();
    
}

// Python code for reference
// def fact_mass_twist_1(M_in, p):
//     M = M_in.copy()  # Avoid modifying original matrix
//     n = M.shape[0]
    
//     A = np.eye(n)

//     for j in range(p, 0, -1):  # Equivalent to MATLAB's p:-1:1

//         if j % 2 == 0:
//             M = M[::-1, ::-1]

//         itwist = 1
//         if j == 1:
//             itwist = min(n - p, n // 2 + 1)

//         L = fact_lower_diag_twist(j, M, itwist)

//         Linv = np.linalg.inv(L)
//         M = Linv @ M @ Linv.T

//         if j % 2 == 0:
//             M = M[::-1, ::-1]
//             L = L[::-1, ::-1]

//         A = A @ L

//     # Adjust diagonal of M
//     for i in range(n):
//         M[i, i] = np.sqrt(abs(M[i, i]))

//     A = A @ M

//     return A

template<class T>
void DFDM::matrix<T>::fact_mass_twist_1(DFDM::matrix<T>& matl, DFDM::matrix<T>& matlt, uint64_t p){
    // Other possible matrix factorisation, seems unstable
    if(this->rows != this->cols){
        std::cerr << "DFDM::Matrix:: Error at line " << __LINE__ << " in function " << __func__ << " in file " << __FILE__ << ": Matrix is not square, cannot perform Schur decomposition." << std::endl;
        exit(1);
    }

    auto M = *this; // I need a deep copy
    
    int n = this->rows;

    DFDM::matrix<T> A(n,n);
    
    matl.resize(n,n);
    matlt.resize(n,n);

    DFDM::matrix<double> revert_rows_and_columns(n,n);
    for (int i = 0; i<n; i++){
        revert_rows_and_columns.value_set(i,n-i-1,1.0);
    }

    for (int i = 0; i<n; i++){
        A.value_set(i,i,1.0);
    }

    for (int j=p; j>0; j--){
        DFDM::matrix<T> L(n,n);
        DFDM::matrix<T> L_inv;
        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
        }

        int itwist = 1;
        if (j == 1){
            itwist = std::min((int)(n - p), n / 2 + 1);
        }

        M.fact_lower_diag_twist(j, L, itwist);
        L.inverse(L_inv);

        M = L_inv.matprod(M.matprod(L_inv.transpose_inplace()));

        if (j % 2 ==0){
            M = revert_rows_and_columns.matprod(M.matprod(revert_rows_and_columns));
            L = revert_rows_and_columns.matprod(L.matprod(revert_rows_and_columns));
        }

        A = A.matprod(L);
    }

    for (int i = 0; i<n; i++){
        M.value_set(i,i,std::sqrt(std::abs(M(i,i))));
    }
    A = A.matprod(M);

    matl=A;
    matlt=A.transpose_inplace();
    
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::vecprod(const DFDM::matrix<T>& vec){
    if(this->banded == true){
        std::cerr << "ERROR: Cannot perform matrix-vector product on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    DFDM::matrix<T> prod(this->rows, vec.cols);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, this->rows, this->cols, 1.0, this->data_vec.data(), this->cols, vec.data_vec.data(), 1, 0.0, prod.data_vec.data(), 1);

    return prod;
}
template<class T>
DFDM::matrix<T> DFDM::matrix<T>::matprod(const DFDM::matrix<T>& matI) const {
    if(this->banded == true || matI.banded == true){
        std::cerr << "ERROR: Cannot perform matrix-matrix product on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

    if(banded == true){
        std::cerr << "ERROR: Cannot perform scalar multiplication on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    DFDM::matrix<T> sprod(*this);
    cblas_dscal(rows*cols, scalar, sprod.data_vec.data(), 1);
    return sprod;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::operator+(const matrix& inm){

    if(banded == true || inm.banded == true){
        std::cerr << "ERROR: Cannot perform addition on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(inm.cols != cols || inm.rows != rows){
        std::cerr << "Dims Don't Match for Matrix/Vector Sum"<<std::endl;
        exit(-1);
    }
    DFDM::matrix<T> ssum(*this);
    cblas_daxpy(rows*cols, 1.0, inm.data_vec.data(), 1,ssum.data_vec.data() , 1);
    return ssum;
}

template<class T>
DFDM::matrix<T> DFDM::matrix<T>::operator-(const matrix& inm){
    if(banded == true || inm.banded == true){
        std::cerr << "ERROR: Cannot perform subtraction on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(inm.cols != cols || inm.rows != rows){
        std::cerr << "Dims Don't Match for Matrix/Vector Sub"<<std::endl;
        exit(-1);
    }
    DFDM::matrix<T> ssum(*this);
    cblas_daxpy(rows*cols, -1.0, inm.data_vec.data(), 1,ssum.data_vec.data() , 1);
    return ssum;
}

template<class T>
void DFDM::matrix<T>::row_copy(const matrix& mat_src, uint32_t src_row_id, matrix& mat_tar, uint32_t tar_row_id){
    if(mat_src.banded == true || mat_tar.banded == true){
        std::cerr << "ERROR: Cannot perform row copy on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

    if(mat_src.banded == true || mat_tar.banded == true){
        std::cerr << "ERROR: Cannot perform col copy on a banded matrix. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(mat_src.rows != mat_tar.rows){
        std::cerr << "Dimension mismatch in col copy" << std::endl;
        exit(-1);
    }
    if(src_col_id >= mat_src.cols || tar_col_id >= mat_tar.cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }

    if constexpr (std::is_same<T,double>::value) {
        const double* src_ptr = mat_src.data_vec.data();
        double* tar_ptr = mat_tar.data_vec.data();
        cblas_dcopy(mat_src.rows,
                    src_ptr + src_col_id,
                    mat_src.cols,
                    tar_ptr + tar_col_id,
                    mat_tar.cols);
    } else {
        auto src_it = mat_src.data_vec.begin();
        std::advance(src_it, src_col_id);

        auto tar_it = mat_tar.data_vec.begin();
        std::advance(tar_it, tar_col_id);
        if(mat_src.cols == 1 && mat_tar.cols == 1){
            std::copy(src_it, src_it + mat_src.rows, tar_it);
        }else{
            for(uint32_t i = 0; i < mat_src.rows; i++){
                *tar_it = *src_it;
                std::advance(src_it, mat_src.cols);
                std::advance(tar_it, mat_tar.cols);
            }
        }
    }
}


template<class T>
void DFDM::matrix<T>::row_set(uint32_t tgt_row, std::vector<T> row_data){

    if(banded == true){
        std::cerr << "ERROR: Cannot set elements of a banded matrix using row_set. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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

    if(banded == true){
        std::cerr << "ERROR: Cannot set elements of a banded matrix using col_set. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(col_data.size() != rows){
        std::cerr << "Dimension mismatch in col_set, size of source matrix is not equal to target matrix" << ", matrix_rows:" << rows << " col length:"<< col_data.size()<< std::endl;
        exit(-1);
    }
    if(tgt_col >= cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }

    if constexpr (std::is_same<T,double>::value) {
        const double* src_ptr = col_data.data();
        double* dst_ptr = data_vec.data();
        cblas_dcopy(rows,
                    src_ptr,
                    1,
                    dst_ptr + tgt_col,
                    cols);
    } else {
        for (uint32_t i = 0; i < rows; i++) {
            data_vec[i * cols + tgt_col] = col_data[i];
        }
    }
}

template<class T>
std::vector<T> DFDM::matrix<T>::col_get(uint32_t col_id){

    if(banded == true){
        std::cerr << "ERROR: Cannot get elements of a banded matrix using col_get. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

    if(col_id >= cols){
        std::cerr << " ERROR:Cols to be copied cannot be larger than total number of cols in target or source matrix!" << std::endl;
        exit(-1);
    }
    std::vector<T> col_data(rows);
    if constexpr (std::is_same<T,double>::value) {
        const double* src_ptr = data_vec.data();
        double* dst_ptr = col_data.data();
        cblas_dcopy(rows,
                    src_ptr + col_id,
                    cols,
                    dst_ptr,
                    1);
    } else {
        for (uint32_t i = 0; i < rows; i++) {
            col_data[i] = data_vec[i * cols + col_id];
        }
    }
    return col_data;
}

template<class T>
std::vector<T> DFDM::matrix<T>::row_get(uint32_t row_id){
    if(banded == true){
        std::cerr << "ERROR: Cannot get elements of a banded matrix using row_get. Please convert it back to a full matrix first." << std::endl;
        exit(-1);
    }

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
    if(banded == true){
        std::cout << "This is a banded matrix, printing in banded matrix format:" << std::endl;

        for(uint32_t i = 0; i < (kd+1); i++){
            for(uint32_t j = 0; j < cols; j++){
                std::cout<< data_vec[i*cols + j] <<"\t";
             }
            std::cout << "\n";
        }

    }else{
        for(uint32_t i = 0; i < rows; i++){
            for(uint32_t j = 0; j < cols; j++){
                std::cout<< data_vec[i*cols + j] <<"\t";
            }
            std::cout << "\n";
        }
    }


}

template<class T>
void DFDM::matrix<T>::print_file(std::string file_name) const{
    std::ofstream file_out(file_name);
    
    if(banded == true){
        std::cout << "This is a banded matrix, printing in banded matrix format:" << std::endl;
    }

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