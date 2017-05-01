//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Configuration.h"

// Always test the fallback code, regardless of MKL definition
#undef HAVE_MKL

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x) std::cout << x; exit(0);

#include "AFQMC/Matrix/SparseMatrix.h"
#include "AFQMC/Matrix/tests/sparse_matrix_helpers.h"

// Include the templates directly so all the needed types get instantiated
//  and so the undef of HAVE_MKL has an effect
#include "AFQMC/Numerics/SparseMatrixOperations.cpp"

#include <stdio.h>
#include <string>
#include <complex>

using std::string;
using std::complex;
using std::cout;


namespace qmcplusplus
{
template<typename T> void test_matrix_mult_1_1()
{
  SparseMatrix<T> A(1,1);

  std::vector<s1D<T>> I;
  I.push_back(s1D<T>(0, 1.0));
  A.initFroms1D(I,true);

  T *B = new T[1];
  T *C = new T[1];
  B[0] = 1.0;
  C[0] = 0.0;
  int nrows = 1;
  SparseMatrixOperators::product_SpMatV(nrows, A, B, C);

  REQUIRE(realPart(C[0]) == Approx(1.0));
}

TEST_CASE("sparse_matrix_mult", "[sparse_matrix]")
{
  test_matrix_mult_1_1<complex<double>>();
}

TEST_CASE("sparse_matrix_mult_2x2", "[sparse_matrix]")
{

  SparseMatrix<complex<double>> A(2,2);

  std::vector< s2D<complex<double>> > I;
  I.push_back(s2D<complex<double>>(0,0, 2.0));
  I.push_back(s2D<complex<double>>(0,1, 4.0));
  // [[2.0 4.0]
  //  [0.0 0.0]]
  A.initFroms2D(I,true);


  complex<double> *B = new complex<double>[2];
  complex<double> *C = new complex<double>[2];
  B[0] = 4.0;
  B[1] = 3.0;
  C[0] = 0.0;
  C[1] = 0.0;
  int nrows = 2;
  SparseMatrixOperators::product_SpMatV(nrows, A, B, C);

  REQUIRE(C[0].real() == Approx(20.0));
  REQUIRE(C[0].imag() == Approx(0.0));

  REQUIRE(C[1].real() == Approx(0.0));
  REQUIRE(C[1].imag() == Approx(0.0));
}

TEST_CASE("sparse_matrix_axpy_1x1", "[sparse_matrix]")
{

  SparseMatrix<complex<double>> A(1,1);

  std::vector< s1D<complex<double>> > I;
  I.push_back(s1D<complex<double>>(0, 1.0));
  A.initFroms1D(I,true);

  //output_matrix(A);

  complex<double> *B = new complex<double>[1];
  complex<double> *C = new complex<double>[1];
  B[0] = 1.0;
  C[0] = 3.0;
  int nrows = 1;
  complex<double> alpha = 2.0;
  complex<double> beta = 2.0;
  SparseMatrixOperators::product_SpMatV(nrows, nrows, alpha, A.values(), A.column_data(), A.row_index(), B, beta, C);

  REQUIRE(C[0].real() == Approx(8.0));
  REQUIRE(C[0].imag() == Approx(0.0));
}

TEST_CASE("sparse_matrix_zaxpy_2x2", "[sparse_matrix]")
{

  SparseMatrix<complex<double>> A(2,2);

  std::vector< s2D<complex<double>> > I;
  I.push_back(s2D<complex<double>>(0,0, 2.0));
  I.push_back(s2D<complex<double>>(0,1, 4.0));
  // [[2.0 4.0]
  //  [0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  complex<double> *B = new complex<double>[2];
  complex<double> *C = new complex<double>[2];
  B[0] = 4.0;
  B[1] = 3.0;
  C[0] = 2.0;
  C[1] = 1.3;
  int nrows = 2;
  complex<double> alpha = 2.0;
  complex<double> beta = 2.0;
  SparseMatrixOperators::product_SpMatV(nrows, nrows, alpha, A.values(), A.column_data(), A.row_index(), B, beta, C);

  REQUIRE(C[0].real() == Approx(44.0));
  REQUIRE(C[0].imag() == Approx(0.0));

  REQUIRE(C[1].real() == Approx(2.6));
  REQUIRE(C[1].imag() == Approx(0.0));
}

TEST_CASE("sparse_matrix_daxpy_2x2", "[sparse_matrix]")
{

  SparseMatrix<double> A(2,2);

  std::vector< s2D<double> > I;
  I.push_back(s2D<double>(0,0, 2.0));
  I.push_back(s2D<double>(0,1, 4.0));
  // [[2.0 4.0]
  //  [0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  double *B = new double[2];
  double *C = new double[2];
  B[0] = 4.0;
  B[1] = 3.0;
  C[0] = 2.0;
  C[1] = 1.3;
  int nrows = 2;
  double alpha = 2.0;
  double beta = 2.0;
  SparseMatrixOperators::product_SpMatV(nrows, nrows, alpha, A.values(), A.column_data(), A.row_index(), B, beta, C);

  REQUIRE(C[0] == Approx(44.0));
  REQUIRE(C[1] == Approx(2.6));
}

TEST_CASE("sparse_matrix_zaxpy_T_2x2", "[sparse_matrix]")
{

  SparseMatrix<complex<double>> A(2,2);

  std::vector< s2D<complex<double>> > I;
  I.push_back(s2D<complex<double>>(0,0, 2.0));
  I.push_back(s2D<complex<double>>(0,1, 4.0));
  // [[2.0 4.0]
  //  [0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  complex<double> *B = new complex<double>[2];
  complex<double> *C = new complex<double>[2];
  B[0] = 4.0;
  B[1] = 3.0;
  C[0] = 2.0;
  C[1] = 1.3;
  int nrows = 2;
  complex<double> alpha = 2.0;
  complex<double> beta = 2.0;
  SparseMatrixOperators::product_SpMatTV(nrows, nrows, alpha, A, B, beta, C);

  for (int i = 0; i < nrows; i++) {
    cout << C[i] << std::endl;
  }

  REQUIRE(C[0].real() == Approx(20.0));
  REQUIRE(C[0].imag() == Approx(0.0));

  REQUIRE(C[1].real() == Approx(34.6));
  REQUIRE(C[1].imag() == Approx(0.0));

  C[0] = 2.0;
  C[1] = 1.3;
  SparseMatrixOperators::product_SpMatTV(nrows, nrows, alpha, A.values(), A.column_data(), A.row_index(), B, beta, C);

  REQUIRE(C[0].real() == Approx(20.0));
  REQUIRE(C[0].imag() == Approx(0.0));

  REQUIRE(C[1].real() == Approx(34.6));
  REQUIRE(C[1].imag() == Approx(0.0));

}

TEST_CASE("sparse_matrix_mm_2x2", "[sparse_matrix]")
{

  SparseMatrix<complex<double>> A(2,2);

  std::vector< s2D<complex<double>> > I;
  I.push_back(s2D<complex<double>>(0,0, 2.0));
  I.push_back(s2D<complex<double>>(0,1, 4.0));
  // [[2.0 4.0]
  //  [0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  complex<double> *B = new complex<double>[4];
  complex<double> *C = new complex<double>[4];
  B[0] = 4.0;
  B[1] = 3.0;
  B[2] = 1.0;
  B[3] = 5.0;

  C[0] = 2.0;
  C[1] = 1.3;
  C[2] = 1.3;
  C[3] = 1.3;
  int nrows = 2;
  complex<double> alpha = 2.0;
  complex<double> beta = 2.0;
  SparseMatrixOperators::product_SpMatM(nrows, nrows, nrows, alpha, A, B, nrows, beta, C, nrows);

  //std::cout << C[0] << std::endl;
  //std::cout << C[1] << std::endl;
  //std::cout << C[2] << std::endl;
  //std::cout << C[3] << std::endl;

  REQUIRE(C[0].real() == Approx(28.0));
  REQUIRE(C[0].imag() == Approx(0.0));

  REQUIRE(C[1].real() == Approx(54.6));
  REQUIRE(C[1].imag() == Approx(0.0));

  REQUIRE(C[2].real() == Approx(2.6));
  REQUIRE(C[2].imag() == Approx(0.0));

  REQUIRE(C[3].real() == Approx(2.6));
  REQUIRE(C[3].imag() == Approx(0.0));
}

TEST_CASE("sparse_matrix_real_mm_2x2", "[sparse_matrix]")
{

  const int N = 2; // nrows of A
  const int M = 3; // ncols of C
  const int K = 4; // ncols of A
  // A is MxK
  // B is KxN
  // C is MxN
  SparseMatrix<double> A(M,K);

  std::vector< s2D<double> > I;
  I.push_back(s2D<double>(0,0, 2.0));
  I.push_back(s2D<double>(0,1, 4.0));
  I.push_back(s2D<double>(1,2, 5.0));
  I.push_back(s2D<double>(2,1, 1.0));
  // [[2.0 4.0 0.0 0.0]
  //  [0.0 0.0 5.0 0.0]
  //  [1.0 0.0 0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  double *B = new double[K*N];
  double *C = new double[M*N];
  B[0] = 4.0;
  B[1] = 3.0;
  B[2] = 1.0;
  B[3] = 5.0;
  B[4] = 0.0;
  B[5] = 0.0;
  B[6] = 0.0;
  B[7] = 0.0;

  C[0] = 2.0;
  C[1] = 1.3;
  C[2] = 1.3;
  C[3] = 1.3;
  C[4] = 0.0;
  C[5] = 0.0;
  double alpha = 2.0;
  double beta = 2.0;
  int ldb = N;
  int ldc = N;
  SparseMatrixOperators::product_SpMatM(M, N, K, alpha, A, B, ldb, beta, C, ldc);

  //for (int i = 0; i < M*N; i++)
  //{
  //  std::cout << C[i] << std::endl;
  //}

  REQUIRE(C[0] == Approx(28.0));
  REQUIRE(C[1] == Approx(54.6));
  REQUIRE(C[2] == Approx(2.6));
  REQUIRE(C[3] == Approx(2.6));
  REQUIRE(C[4] == Approx(2.0));
  REQUIRE(C[5] == Approx(10.0));

}


TEST_CASE("sparse_matrix_float_mm_2x2", "[sparse_matrix]")
{

  const int N = 2; // nrows of A
  const int M = 3; // ncols of C
  const int K = 4; // ncols of A
  // A is MxK
  // B is KxN
  // C is MxN
  SparseMatrix<float> A(M,K);

  std::vector< s2D<double> > I;
  I.push_back(s2D<double>(0,0, 2.0));
  I.push_back(s2D<double>(0,1, 4.0));
  I.push_back(s2D<double>(0,3, 1.0));
  I.push_back(s2D<double>(1,2, 5.0));
  I.push_back(s2D<double>(2,1, 1.0));
  // [[2.0 4.0 0.0 1.0]
  //  [0.0 0.0 5.0 0.0]
  //  [1.0 0.0 0.0 0.0]]
  A.initFroms2D(I,true);

  //output_matrix(A);

  float *B = new float[K*N];
  float *C = new float[M*N];
  for (int i = 0; i < K*N; i++) B[i] = 0;
  for (int i = 0; i < M*N; i++) C[i] = 0;
  B[0] = 4.0;
  B[1] = 3.0;
  B[2] = 1.0;
  B[3] = 5.0;
  B[6] = 1.0;

  C[0] = 2.0;
  C[1] = 1.3;
  C[2] = 1.3;
  C[3] = 1.3;
  float alpha = 2.0;
  float beta = 2.0;
  int ldb = N;
  int ldc = N;
  SparseMatrixOperators::product_SpMatM(M, N, K, alpha, A, B, ldb, beta, C, ldc);

  //for (int i = 0; i < M*N; i++)
  //{
  //  std::cout << C[i] << std::endl;
  //}

  REQUIRE(C[0] == Approx(30.0));
  REQUIRE(C[1] == Approx(54.6));
  REQUIRE(C[2] == Approx(2.6));
  REQUIRE(C[3] == Approx(2.6));
  REQUIRE(C[4] == Approx(2.0));
  REQUIRE(C[5] == Approx(10.0));

}

#include "sparse_mult_cases.cpp"

}

