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
#include "AFQMC/Matrix/SparseMatrix.h"

#include "AFQMC/Matrix/tests/sparse_matrix_helpers.h"
#include <stdio.h>
#include <string>
#include <complex>


using std::string;
using std::complex;

namespace qmcplusplus
{

template<typename T> void test_init_one_D()
{
  SparseMatrix<T> M(1,1);

  // s1D is tuple<uint32_t, T>
  std::vector< s1D<T> > I;
  I.push_back(s1D<T>(0, 1.0));
  M.initFroms1D(I,true);

  int idx = M.find_element(0, 0);
  REQUIRE(idx == 0);
  const T &val = M(0,0);
  REQUIRE(realPart(val) == Approx(1.0));
}

TEST_CASE("sparse_matrix_init_one_D", "[sparse_matrix]")
{
  test_init_one_D<double>();
  test_init_one_D<complex<double>>();
}

template<typename T> void test_init_one_D_4()
{
  SparseMatrix<T> M(1,4);

  // s1D is tuple<uint32_t, T>
  std::vector< s1D<T> > I;
  I.push_back(s1D<T>(1, 4.0));
  I.push_back(s1D<T>(2, 2.0));
  M.initFroms1D(I,true);

  //output_matrix(M);

  int idx = M.find_element(0, 0);
  REQUIRE(idx == -1);
  idx = M.find_element(0, 1);
  REQUIRE(idx == 0);
  idx = M.find_element(0, 2);
  REQUIRE(idx == 1);
  idx = M.find_element(0, 3);
  REQUIRE(idx == -1);

  const T &val0 = M(0,0);
  REQUIRE(realPart(val0) == Approx(0.0));
  const T &val1 = M(0,1);
  REQUIRE(realPart(val1) == Approx(4.0));
  const T &val2 = M(0,2);
  REQUIRE(realPart(val2) == Approx(2.0));
  const T &val3 = M(0,3);
  REQUIRE(realPart(val3) == Approx(0.0));
}

TEST_CASE("sparse_matrix_init_one_D_4", "[sparse_matrix]")
{
  test_init_one_D_4<double>();
  test_init_one_D_4<complex<double> >();
}


template<typename T> void test_init_two_D_4()
{
  SparseMatrix<T> M(1,4);

  // s2D is tuple<uint32_t, T>
  std::vector< s2D<T> > I;
  I.push_back(s2D<T>(0, 0, 4.0));
  M.initFroms2D(I,true);

  //output_matrix(M);

  int idx = M.find_element(0, 0);
  REQUIRE(idx == 0);

  T val0 = M(0,0);
  REQUIRE(realPart(val0) == Approx(4.0));
}

TEST_CASE("sparse_matrix_init_two_D_4", "[sparse_matrix]")
{
  test_init_two_D_4<double>();
  test_init_two_D_4<complex<double>>();
}

template<typename T> void test_matrix_init_two_D_22()
{
  SparseMatrix<T> M(2,2);

  // s1D is tuple<uint32_t, uint32_t, T>
  std::vector< s2D<T> > I;
  I.push_back(s2D<T>(0, 1, 4.0));
  I.push_back(s2D<T>(1, 1, 2.0));
  M.initFroms2D(I,true);

  //output_matrix(M);

  int idx = M.find_element(0, 0);
  REQUIRE(idx == -1);
  idx = M.find_element(0, 1);
  REQUIRE(idx == 0);
  idx = M.find_element(1, 0);
  REQUIRE(idx == -1);
  idx = M.find_element(1, 1);
  REQUIRE(idx == 1);

  T &val0 = M(0,0);
  REQUIRE(realPart(val0) == Approx(0.0));
  T &val1 = M(0,1);
  REQUIRE(realPart(val1) == Approx(4.0));
  T &val2 = M(1,0);
  REQUIRE(realPart(val2) == Approx(0.0));
  T &val3 = M(1,1);
  REQUIRE(realPart(val3) == Approx(2.0));
}

TEST_CASE("sparse_matrix_init_two_D_2x2", "[sparse_matrix]")
{
  test_matrix_init_two_D_22<double>();
  test_matrix_init_two_D_22<complex<double>>();
}


template<typename T> void test_matrix_invariant()
{
  for (int m = 1; m < 4; m++) {
    for (int n = 1; n < 4; n++) {
      std::vector< s2D<T> > I;
      T *A = new T[m*n];
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          // Should use random selection of zero elements
          if ((i+j)%3 == 0)
          //if (false) // all non-zero (dense)
          {
            A[i*n+j] = 0;
          } else {
            T val = 1.0*(i+j+1);
            I.push_back(s2D<T>(i, j, val));
            A[i*n+j] = val;
          }
        }
      }
      SparseMatrix<T> M(m,n);
      M.initFroms2D(I,true);
#if 0
      printf("m = %d  n = %d\n",m,n);
      output_data(A, m*n);
      printf("CSR:\n");
      output_matrix(M);
      printf("\n");
#endif

      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          T tmp_M = M(i,j);
          T tmp_A = A[i*n+j];
          REQUIRE(tmp_M == tmp_A);
        }
      }
    }
  }
}

TEST_CASE("sparse_matrix_invariant", "[sparse_matrix]")
{
  test_matrix_invariant<double>();
  test_matrix_invariant<complex<double>>();
}

}

