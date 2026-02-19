#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "matrix.hpp"

// Minimal test helpers (double-only)
static int g_failures = 0;

inline bool near(double a, double b, double atol=1e-10, double rtol=1e-10) {
	double diff = std::fabs(a - b);
	double thresh = atol + rtol * std::max(std::fabs(a), std::fabs(b));
	return diff <= thresh;
}

inline void assert_true(bool cond, const char* msg) {
	if (!cond) {
		std::cerr << "[FAIL] " << msg << std::endl;
		++g_failures;
	}
}

inline void assert_near(double a, double b, const char* msg, double tol=1e-10) {
	if (!near(a, b, tol, tol)) {
		std::cerr << "[FAIL] " << msg << ": got=" << a << ", expected=" << b << std::endl;
		++g_failures;
	}
}

static void assert_matrix_near(const DFDM::matrix<double>& A, const DFDM::matrix<double>& B, const char* name, double tol=1e-10) {
	assert_true(A.rows == B.rows && A.cols == B.cols, "matrix dims mismatch in assert_matrix_near");
	for (uint32_t i = 0; i < A.rows; ++i) {
		for (uint32_t j = 0; j < A.cols; ++j) {
			const double va = A.data_vec[i*A.cols + j];
			const double vb = B.data_vec[i*B.cols + j];
			if (!near(va, vb, tol, tol)) {
				std::cerr << "[FAIL] " << name << " at (" << i << "," << j << ") got=" << va << ", expected=" << vb << std::endl;
				++g_failures;
				return;
			}
		}
	}
}

// Test basic math: scalar multiply, add, subtract
static void test_basic_arith() {
	// Construct a small 2x3 matrix with explicit values
	DFDM::matrix<double> A(2,3);
	A.value_set(0,0,1.0); A.value_set(0,1,2.0); A.value_set(0,2,3.0);
	A.value_set(1,0,4.0); A.value_set(1,1,5.0); A.value_set(1,2,6.0);

	// Scalar multiply
	auto B = A * 2.0;
	DFDM::matrix<double> Bexp(2,3);
	for (uint32_t i=0;i<2;i++) 
        for (uint32_t j=0;j<3;j++)
		    Bexp.value_set(i,j, A.data_vec[i*A.cols + j] * 2.0);
	assert_matrix_near(B, Bexp, "scalar multiply");

	// Add and subtract
	auto C = A + B;
	auto D = C - A;
	assert_matrix_near(D, B, "add/subtract consistency");
}

// Test matrix-matrix product against a hand-computed result
static void test_matprod() {
	DFDM::matrix<double> A(2,3);
	// A = [1 2 3; 4 5 6]
	A.value_set(0,0,1.0); A.value_set(0,1,2.0); A.value_set(0,2,3.0);
	A.value_set(1,0,4.0); A.value_set(1,1,5.0); A.value_set(1,2,6.0);

	DFDM::matrix<double> E(3,2);
	// E = [7 8; 9 10; 11 12]
	E.value_set(0,0,7.0);  E.value_set(0,1,8.0);
	E.value_set(1,0,9.0);  E.value_set(1,1,10.0);
	E.value_set(2,0,11.0); E.value_set(2,1,12.0);

	auto P = A.matprod(E);

	DFDM::matrix<double> Pexp(2,2);
	// Pexp = A*E = [58 64; 139 154]
	Pexp.value_set(0,0,58.0);  Pexp.value_set(0,1,64.0);
	Pexp.value_set(1,0,139.0); Pexp.value_set(1,1,154.0);

	assert_matrix_near(P, Pexp, "matprod accuracy");
}

// Test matrix-vector product (vec is represented as an n x 1 matrix)
static void test_vecprod() {
	DFDM::matrix<double> A(2,2);
	// A = [2 1; 0 3]
	A.value_set(0,0,2.0); A.value_set(0,1,1.0);
	A.value_set(1,0,0.0); A.value_set(1,1,3.0);

	DFDM::matrix<double> x(2,1);
	x.value_set(0,0,4.0);
	x.value_set(1,0,5.0);

	auto y = A.vecprod(x);
	// y = [2*4 + 1*5; 0*4 + 3*5] = [13; 15]
	assert_true(y.rows == 2 && y.cols == 1, "vecprod dims");
	assert_near(y.data_vec[0], 13.0, "vecprod y0");
	assert_near(y.data_vec[1], 15.0, "vecprod y1");
}

// Test transpose operations
static void test_transpose() {
	DFDM::matrix<double> A(2,3);
	// Fill with row-major increasing values 1..6
	double v=1.0;
	for (uint32_t i=0;i<2;i++) for (uint32_t j=0;j<3;j++) A.value_set(i,j,v++);

	DFDM::matrix<double> AT;
	A.transpose(AT);
	auto AT2 = A.transpose_inplace();
	assert_true(AT.rows == 3 && AT.cols == 2, "transpose dims");
	assert_matrix_near(AT, AT2, "transpose vs transpose_inplace");

	// Check entry mapping
	for(uint32_t i=0;i<2;i++) for(uint32_t j=0;j<3;j++) {
		assert_near(AT.data_vec[j*AT.cols + i], A.data_vec[i*A.cols + j], "transpose value");
	}
}

// Test elementwise operations (add, mul, div)
static void test_elementwise() {
	DFDM::matrix<double> A(2,2), B(2,2);
	// A = [1 2; 3 4], B = [5 6; 7 8]
	double avals[4] = {1,2,3,4};
	double bvals[4] = {5,6,7,8};
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
		A.value_set(i,j, avals[i*2+j]);
		B.value_set(i,j, bvals[i*2+j]);
	}

	auto Cadd = DFDM::matrix<double>::elementwise_add(A,B);
	auto Cmul = DFDM::matrix<double>::elementwise_mult(A,B);
	auto Cdiv = DFDM::matrix<double>::elementwise_div(B,A);

	for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
		double a = avals[i*2+j];
		double b = bvals[i*2+j];
		assert_near(Cadd.data_vec[i*2+j], a+b, "elem add");
		assert_near(Cmul.data_vec[i*2+j], a*b, "elem mul");
		assert_near(Cdiv.data_vec[i*2+j], b/a, "elem div");
	}
}

// Test row/col get/set and copy
static void test_row_col_ops() {
	DFDM::matrix<double> A(3,3);
	// A = diag(1,2,3)
	for (int i=0;i<3;i++) A.value_set(i,i, (double)(i+1));

	// col_get/col_set
	auto col1 = A.col_get(1);
	assert_true(col1.size()==3, "col_get size");
	assert_near(col1[0], 0.0, "col1[0]");
	assert_near(col1[1], 2.0, "col1[1]");
	assert_near(col1[2], 0.0, "col1[2]");

	std::vector<double> newcol = {10,20,30};
	A.col_set(2, newcol);
	auto col2 = A.col_get(2);
	for (int i=0;i<3;i++) assert_near(col2[i], newcol[i], "col_set accuracy");

	// row_get/row_set
	auto row0 = A.row_get(0);
    assert_true(row0.size()==3, "row_get size");
    assert_near(row0[0], 1.0, "row0[0]");
    assert_near(row0[1], 0.0, "row0[1]");
    assert_near(row0[2], 10.0, "row0[2]");



	std::vector<double> newrow = {7,8,9};
	A.row_set(1, newrow);
	auto row1 = A.row_get(1);
	for (int j=0;j<3;j++) assert_near(row1[j], newrow[j], "row_set accuracy");

	// row_copy / col_copy
	DFDM::matrix<double> B(3,3);
	DFDM::matrix<double>::row_copy(A, 1, B, 0);
	auto b0 = B.row_get(0);
	for (int j=0;j<3;j++) assert_near(b0[j], newrow[j], "row_copy accuracy");

	DFDM::matrix<double>::col_copy(A, 2, B, 1);
	auto bcol1 = B.col_get(1);
	auto target = A.col_get(2);
	for (int i=0;i<3;i++) assert_near(bcol1[i], target[i], "col_copy accuracy");
}

// Test triangular solve (lower)
static void test_trsm_lower() {
	// L is 3x3 lower-triangular with 1s on diagonal
	DFDM::matrix<double> L(3,3);
	L.value_set(0,0,1.0);
	L.value_set(1,0,2.0); L.value_set(1,1,1.0);
	L.value_set(2,0,3.0); L.value_set(2,1,4.0); L.value_set(2,2,1.0);

	// Solve L X = B where B is 3x2
	DFDM::matrix<double> Bm(3,2);
	Bm.value_set(0,0,1.0); Bm.value_set(0,1,2.0);
	Bm.value_set(1,0,3.0); Bm.value_set(1,1,4.0);
	Bm.value_set(2,0,5.0); Bm.value_set(2,1,6.0);

	DFDM::matrix<double> X;
	const char uplo = 'L';
	DFDM::matrix<double>::triangular_solve(L, X, Bm, &uplo);

	// Verify L*X ≈ original B (since X overwrote Bm inside call and returned into X)
	auto BX = L.matprod(X);
	DFDM::matrix<double> Borig(3,2);
	Borig.value_set(0,0,1.0); Borig.value_set(0,1,2.0);
	Borig.value_set(1,0,3.0); Borig.value_set(1,1,4.0);
	Borig.value_set(2,0,5.0); Borig.value_set(2,1,6.0);
	assert_matrix_near(BX, Borig, "triangular_solve lower");
}

// Test banded triangular solve (lower)
static void test_banded_trsm_lower() {
	// Use a strictly lower-triangular matrix with kd >= (n-1) so full lower is in band
	DFDM::matrix<double> L(3,3);
	L.value_set(0,0,2.0);
	L.value_set(1,0,1.0); L.value_set(1,1,3.0);
	L.value_set(2,0,0.5); L.value_set(2,1,2.0); L.value_set(2,2,4.0);

	int kd = 2; // capture all sub-diagonals for 3x3
	const char uplo = 'L';
	auto Lb = L.toBanded(kd, &uplo);

	DFDM::matrix<double> Bm(3,1);
	Bm.value_set(0,0,1.0);
	Bm.value_set(1,0,2.0);
	Bm.value_set(2,0,3.0);

	DFDM::matrix<double> X;
	DFDM::matrix<double>::banded_triangular_solve(Lb, X, Bm, &uplo);
	auto BX = L.matprod(X);

	DFDM::matrix<double> Borig(3,1);
	Borig.value_set(0,0,1.0); Borig.value_set(1,0,2.0); Borig.value_set(2,0,3.0);
	assert_matrix_near(BX, Borig, "banded triangular_solve lower");
}

// Test inverse: A * inv(A) ≈ I
static void test_inverse() {
	DFDM::matrix<double> A(3,3);
	// A is well-conditioned small matrix
	double vals[9] = {4,1,2, 0,3, -1, 2,0,1};
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) A.value_set(i,j, vals[i*3+j]);

	DFDM::matrix<double> Ainv;
	A.inverse(Ainv);
	auto I = A.matprod(Ainv);

	// Compare to identity
	for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
		double expect = (i==j) ? 1.0 : 0.0;
		assert_near(I.data_vec[i*3+j], expect, "inverse identity check", 1e-8);
	}
}

// Test Cholesky: A = L * L^T
static void test_cholesky() {
	DFDM::matrix<double> A(2,2);
	// SPD matrix
	// [4 2; 2 3]
	A.value_set(0,0,4.0); A.value_set(0,1,2.0);
	A.value_set(1,0,2.0); A.value_set(1,1,3.0);

	DFDM::matrix<double> L, LT;
	A.cholesky(L, LT);

	auto R = L.matprod(LT);
	assert_matrix_near(R, A, "cholesky reconstruction", 1e-8);
}

// Test sqrtm_schur on a diagonal positive matrix (sqrt is trivial)
static void test_sqrtm_schur() {
	DFDM::matrix<double> A(2,2);
	// diag(4, 9)
	A.value_set(0,0,4.0); A.value_set(0,1,0.0);
	A.value_set(1,0,0.0); A.value_set(1,1,9.0);

	DFDM::matrix<double> S, ST;
	A.sqrtm_schur(S, ST);
	// Expect diag(2,3)
	DFDM::matrix<double> Sexp(2,2);
	Sexp.value_set(0,0,2.0); Sexp.value_set(0,1,0.0);
	Sexp.value_set(1,0,0.0); Sexp.value_set(1,1,3.0);

	assert_matrix_near(S, Sexp, "sqrtm_schur diagonal");
	auto back = S.matprod(ST);
	assert_matrix_near(back, A, "sqrtm_schur reconstruction");
}

// Test toBanded mapping for lower/upper storage integrity
static void test_toBanded_layout() {
	// 4x4 with distinct entries to verify placement
	DFDM::matrix<double> A(4,4);
	double v = 1.0;
	for (uint32_t i=0;i<4;i++) for (uint32_t j=0;j<4;j++) A.value_set(i,j,v++);

	// Lower band with kd=1
	const char Luplo = 'L';
	auto Lb = A.toBanded(1, &Luplo);
	// Check a few positions: AB(i-j, j) = A(i,j) for j <= i <= min(n-1, j+kd)
	int n = 4;
	for (int j=0;j<n;j++) {
		int i_end = std::min(n-1, j+1);
		for (int i=j;i<=i_end;i++) {
			int ab_row = i - j;
			double ab = Lb.data_vec[ab_row * n + j];
			double a = A.data_vec[i*4 + j];
			assert_near(ab, a, "toBanded L mapping");
		}
	}

	// Upper band with kd=1
	const char Uuplo = 'U';
	auto Ub = A.toBanded(1, &Uuplo);
	// AB(kd+1+i-j, j) = A(i,j) for max(0,j-kd) <= i <= j
	for (int j=0;j<n;j++) {
		int i_start = std::max(0, j-1);
		for (int i=i_start;i<=j;i++) {
			int ab_row = 1 + i - j;
			double ab = Ub.data_vec[ab_row * n + j];
			double a = A.data_vec[i*4 + j];
			assert_near(ab, a, "toBanded U mapping");
		}
	}
}

int main() {
	std::cout << "Running DFDM matrix accuracy tests (double)..." << std::endl;

	test_basic_arith();
	test_matprod();
	test_vecprod();
	test_transpose();
	test_elementwise();
	test_row_col_ops();
	test_trsm_lower();
	test_banded_trsm_lower();
	test_inverse();
	test_cholesky();
	test_sqrtm_schur();
	test_toBanded_layout();

	if (g_failures == 0) {
		std::cout << "All tests passed." << std::endl;
		return EXIT_SUCCESS;
	} else {
		std::cerr << g_failures << " test(s) failed." << std::endl;
		return EXIT_FAILURE;
	}
}
