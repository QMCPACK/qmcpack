//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by: 
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

// Always test the fallback code, regardless of MKL definition
#undef HAVE_MKL
#define MKL_INT         int
#define MKL_Complex8    std::complex<float> 
#define MKL_Complex16   std::complex<double>

#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x; exit(0);}

#include<iostream>
#include<vector>

#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"

#include<boost/multi_array.hpp>

using std::vector;
using boost::extents;
using boost::indices;
using range_t = boost::multi_array_types::index_range;
using boost::multi_array;
using boost::multi_array_ref;
using csr_matrix = ma::sparse::csr_matrix<double,int,int>;

namespace qmcplusplus
{

// tests dispatching through ma_operations
void test_sparse_matrix_mult_native()
{
    csr_matrix A( {4, 4}, {0,0}, 4);
    A[3][3] = 1.;
    A[2][1] = 3.;
    A[0][1] = 9.;

    // matrix-matrix
    {
	vector<double> b = {
		1.,2.,1., 5.,
		2.,5.,8., 7.,
		1.,8.,9., 9.,
		4.,1.,2., 3.
	};
	multi_array_ref<double, 2> B(b.data(), extents[4][4]);
	REQUIRE(B.num_elements() == b.size());

	vector<double> c(16);
	multi_array_ref<double, 2> C(c.data(), extents[4][4]);
	REQUIRE(C.num_elements() == c.size());

	ma::product(A, B, C); // C = A*B

	vector<double> c2 = {
		18., 45., 72., 63.,
		0., 0., 0., 0., 
		6., 15., 24., 21., 
		4., 1., 2., 3.
	};
	multi_array_ref<double, 2> C2(c2.data(), extents[4][4]);
	REQUIRE(C2.num_elements() == c2.size());
	verify_approx(C, C2);

	vector<double> d(16);
	multi_array_ref<double, 2> D(d.data(), extents[4][4]);
	REQUIRE(D.num_elements() == d.size());

        using ma::T;
	ma::product(T(A), B, D); // D = T(A)*B
	vector<double> d2 = {
		0, 0, 0, 0, 
		12, 42, 36, 72, 
		0, 0, 0, 0, 
		4, 1, 2, 3 
	};
	multi_array_ref<double, 2> D2(d2.data(), extents[4][4]);
	REQUIRE(D2.num_elements() == d2.size());
	verify_approx(D2, D);

    }
    
    // matrix-vector
    {
        vector<double> b = {1., 2., 1., 4.};
        multi_array_ref<double, 1> B(b.data(), extents[4]);
        REQUIRE(B.num_elements() == b.size());

        vector<double> c(4);
        multi_array_ref<double, 1> C(c.data(), extents[4]);
        REQUIRE(C.num_elements() == c.size());

        ma::product(A, B, C); // C = A*B

        vector<double> c2 = { 18., 0., 6., 4.};
        multi_array_ref<double, 1> C2(c2.data(), extents[4]);
        REQUIRE(C2.num_elements() == c2.size());
        verify_approx(C, C2);

        vector<double> d(4);
        multi_array_ref<double, 1> D(d.data(), extents[4]);
        REQUIRE(D.num_elements() == d.size());

        using ma::T;
        ma::product(T(A), B, D); // D = T(A)*B
        vector<double> d2 = { 0., 12., 0., 4.};
        multi_array_ref<double, 1> D2(d2.data(), extents[4]);
        REQUIRE(D2.num_elements() == d2.size());
        verify_approx(D2, D);
   } 

}

TEST_CASE("sparse_ma_operations_native", "[matrix_operations]")
{
  test_sparse_matrix_mult_native();
}

}
