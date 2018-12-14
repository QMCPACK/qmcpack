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


#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Configuration.h"

// Always test the fallback code, regardless of MKL definition
//#undef HAVE_MKL
#define MKL_INT         int
#define MKL_Complex8    std::complex<float> 
#define MKL_Complex16   std::complex<double>

// Avoid the need to link with other libraries just to get APP_ABORT
#undef APP_ABORT
#define APP_ABORT(x) {std::cout << x; exit(0);}

#include "AFQMC/Matrix/tests/matrix_helpers.h"

// Include the templates directly so all the needed types get instantiated
//  and so the undef of HAVE_MKL has an effect
#include "AFQMC/Numerics/ma_operations.hpp"

#include <stdio.h>
#include <string>
#include <complex>
#include <vector>
#include<boost/multi_array.hpp>

using std::string;
using std::complex;
using std::cout;
using std::endl;
using std::vector;
using boost::extents;
using boost::indices;
using range_t = boost::multi_array_types::index_range;
using boost::multi_array;
using boost::multi_array_ref;

namespace qmcplusplus
{

// tests dispatching through ma_operations
void test_dense_matrix_mult()
{
	{
		vector<double> m = {
			9.,24.,30.,
			4.,10.,12.,
			14.,16.,36.//,
		//	9., 6., 1.
		};
		multi_array_ref<double, 2> M(m.data(), extents[3][3]);
		REQUIRE(M.num_elements() == m.size());
		vector<double> x = {1.,2.,3.};
		multi_array_ref<double, 1> X(x.data(), extents[x.size()]);
		vector<double> y(3);
		multi_array_ref<double, 1> Y(y.data(), extents[y.size()]);

		using ma::T;
		ma::product(M, X, Y); // Y := M X
		
		vector<double> mx = {147., 60.,154.};
		multi_array_ref<double, 1> MX(mx.data(), extents[mx.size()]);
                verify_approx(MX, Y);
	}
	{
		vector<double> m = {
			9.,24.,30., 2.,
			4.,10.,12., 1.,
			14.,16.,36., 20.
		};
		multi_array_ref<double, 2> M(m.data(), extents[3][4]);
		REQUIRE(M.num_elements() == m.size());
		vector<double> x = {1.,2.,3., 4.};
		multi_array_ref<double, 1> X(x.data(), extents[x.size()]);
		vector<double> y(3);
		multi_array_ref<double, 1> Y(y.data(), extents[y.size()]);

		using ma::T;
		ma::product(M, X, Y); // Y := M X

		vector<double> mx = {155., 64.,234.};
		multi_array_ref<double, 1> MX(mx.data(), extents[mx.size()]);
		verify_approx( MX, Y );
	}
	{
		vector<double> m = {
			9.,24.,30., 2.,
			4.,10.,12., 1.,
			14.,16.,36., 20.
		};
		multi_array_ref<double, 2> M(m.data(), extents[3][4]);
		REQUIRE(M.num_elements() == m.size());
		vector<double> x = {1.,2.,3.};
		multi_array_ref<double, 1> X(x.data(), extents[x.size()]);
		vector<double> y(4);
		multi_array_ref<double, 1> Y(y.data(), extents[y.size()]);

		using ma::T;
		ma::product(T(M), X, Y); // Y := T(M) X
		
		vector<double> mx = {59., 92., 162., 64.};
		multi_array_ref<double, 1> MX(mx.data(), extents[mx.size()]);
		verify_approx( MX, Y );
	}
	{
		vector<double> m = {
			9.,24.,30., 9.,
			4.,10.,12., 7.,
			14.,16.,36., 1.
		};
		multi_array_ref<double, 2> M(m.data(), extents[3][4]);
		vector<double> x = {1.,2.,3., 4.};
		multi_array_ref<double, 1> X(x.data(), extents[x.size()]);
		vector<double> y = {4.,5.,6.};
		multi_array_ref<double, 1> Y(y.data(), extents[y.size()]);
		ma::product(M, X, Y); // y := M x
		
		vector<double> y2 = {183., 88.,158.};
		multi_array_ref<double, 1> Y2(y2.data(), extents[y2.size()]);
		verify_approx( Y, Y2 );
	}

	{
	vector<double> m = {
		1.,2.,1.,
		2.,5.,8.,
		1.,8.,9.
	};
	multi_array_ref<double, 2> M(m.data(), extents[3][3]);
	REQUIRE( ma::is_hermitian(M) );
	}{
	vector<double> m = {
		1.,0.  , 2.,0. ,  1.,0.,
		2.,0.  , 5.,0. ,  8.,-1.,
		1.,0.  , 8.,1. ,  9.,0.,
	};
	multi_array_ref<complex<double>, 2> M(reinterpret_cast<complex<double>*>(m.data()), extents[3][3]);
	REQUIRE( ma::is_hermitian(M) );
	}{
	vector<double> m = {
		1.,2.,1.,
		2.,5.,8.,
		1.,8.,9.
	};
	multi_array_ref<double, 2> M(m.data(), extents[3][3]);
	REQUIRE( ma::is_hermitian(M) );
	}
	{
	vector<double> a = {
		1.,0.,1.,
		3.,5.,8., 
		4.,8.,9.
	};
	multi_array_ref<double, 2> A(a.data(), extents[3][3]);
	REQUIRE( A.num_elements() == a.size() );
	vector<double> b = {
		6.,2.,8.,
		9.,5.,5.,
		1.,7.,9.
	};
	multi_array_ref<double, 2> B(b.data(), extents[3][3]);
	REQUIRE( B.num_elements() == b.size() );

	vector<double> c(9);
	multi_array_ref<double, 2> C(c.data(), extents[3][3]);
	REQUIRE( C.num_elements() == c.size() );
	
	ma::product(A, B, C);

	vector<double> ab = {
		7., 9., 17.,
		71., 87., 121.,
		105., 111., 153.
	};
	multi_array_ref<double, 2> AB(ab.data(), extents[3][3]);
	REQUIRE( AB.num_elements() == ab.size() );

	verify_approx(C, AB);

	using ma::N;
	ma::product(N(A), N(B), C); // same as ma::product(A, B, C);
	verify_approx(C, AB);

	using ma::T;
	
	ma::product(T(A), B, C);
	vector<double> atb = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
	multi_array_ref<double, 2> AtB(atb.data(), extents[3][3]);
	verify_approx(C, AtB);
	
	ma::product(A, T(B), C);
	vector<double> abt = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
	multi_array_ref<double, 2> ABt(abt.data(), extents[3][3]);
	verify_approx(C, ABt);

	ma::product(T(A), T(B), C);
	vector<double> atbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
	multi_array_ref<double, 2> AtBt(atbt.data(), extents[3][3]);
	verify_approx(C, AtBt);

	using ma::H;
        ma::product(H(A), T(B), C);
        vector<double> ahbt = {44., 44., 58., 74., 65., 107., 94., 94., 138.};
        multi_array_ref<double, 2> AhBt(ahbt.data(), extents[3][3]);
        verify_approx(C, AhBt);

        ma::product(A, H(B), C);
        vector<double> abh = {14., 14., 10., 92., 92., 110., 112., 121., 141.};
        multi_array_ref<double, 2> ABh(abh.data(), extents[3][3]);
        verify_approx(C, ABh);
	
	}
	{
		vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
		multi_array_ref<double, 2> A(a.data(), extents[3][3]);
		REQUIRE(A.num_elements() == a.size());
		multi_array<double, 2> B = A;
		ma::invert(A);

		multi_array<double, 2> Id(extents[3][3]);
		ma::set_identity(Id);

		multi_array<double, 2> Id2(extents[3][3]);
		ma::product(A, B, Id2);
						
		verify_approx(Id, Id2);
	}
        {
                std::vector<double> WORK;
                boost::multi_array<double,1> TAU(extents[3]);

                vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
                multi_array_ref<double, 2> A(a.data(), extents[3][3]);
                REQUIRE(A.num_elements() == a.size());
                WORK.reserve(  ma::gelqf_optimal_workspace_size(A) );
                WORK.reserve(  ma::glq_optimal_workspace_size(A) );
                ma::gelqf(A,TAU,WORK);
                ma::glq(A,TAU,WORK);

                multi_array<double, 2> Id(extents[3][3]);
                ma::set_identity(Id);

                using ma::H;
                multi_array<double, 2> Id2(extents[3][3]);
                ma::product(H(A), A, Id2);

                verify_approx(Id, Id2);
        }
        {
                std::vector<double> WORK;
                boost::multi_array<double,1> TAU(extents[4]);

                vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.,10.,23.,35.};
                multi_array_ref<double, 2> A(a.data(), extents[4][3]);
                REQUIRE(A.num_elements() == a.size());
                WORK.reserve(  ma::gelqf_optimal_workspace_size(A) );
                WORK.reserve(  ma::glq_optimal_workspace_size(A) );
                ma::gelqf(A,TAU,WORK);
                ma::glq(A,TAU,WORK);

                multi_array<double, 2> Id(extents[3][3]);
                ma::set_identity(Id);

                using ma::H;
                multi_array<double, 2> Id2(extents[3][3]);
                ma::product(H(A), A, Id2);

                verify_approx(Id, Id2);
        }
        {
                std::vector<double> WORK;
                boost::multi_array<double,1> TAU(extents[3]);

                vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.};
                multi_array_ref<double, 2> A(a.data(), extents[3][3]);
                REQUIRE(A.num_elements() == a.size());
                WORK.reserve(  ma::geqrf_optimal_workspace_size(A) );
                WORK.reserve(  ma::gqr_optimal_workspace_size(A) );
                ma::geqrf(A,TAU,WORK);
                ma::gqr(A,TAU,WORK);

                multi_array<double, 2> Id(extents[3][3]);
                ma::set_identity(Id);

                using ma::H;
                multi_array<double, 2> Id2(extents[3][3]);
                ma::product(H(A), A, Id2);

                verify_approx(Id, Id2);
        }
        {
                std::vector<double> WORK;
                boost::multi_array<double,1> TAU(extents[4]);

                vector<double> a = {37., 45., 59., 53., 81., 97., 87., 105., 129.,10.,23.,35.};
                multi_array_ref<double, 2> A(a.data(), extents[3][4]);
                REQUIRE(A.num_elements() == a.size());
                WORK.reserve(  ma::geqrf_optimal_workspace_size(A) );
                WORK.reserve(  ma::gqr_optimal_workspace_size(A) );
                ma::geqrf(A,TAU,WORK);
                ma::gqr(A,TAU,WORK);

                multi_array<double, 2> Id(extents[3][3]);
                ma::set_identity(Id);

                using ma::H;
                multi_array<double, 2> Id2(extents[3][3]);
                ma::product(A, H(A), Id2);

                verify_approx(Id, Id2);
        }
        {
                vector<double> a = {
                        9.,24.,30., 45.,
                        4.,10.,12., 12.
                };
                multi_array_ref<double, 2> A(a.data(), extents[2][4]);
                vector<double> at = {
                        9.,4.,
                        24.,10.,
                        30.,12.,
                        45.,12.
                };
                multi_array_ref<double, 2> AT(at.data(), extents[4][2]);
                multi_array<double, 2> B(extents[4][2]);
                ma::transpose(A,B);
                verify_approx( AT, B );
        }
        {
                using namespace std::complex_literals;
                vector<std::complex<double>> m_a = {
                    1.90000,   1.40000 + 0.90000i,   0.40000 + 0.80000i,
                    1.40000 - 0.90000i,   0.20000,   2.20000 + 0.60000i,
                    0.40000 - 0.80000i,   2.20000 - 0.60000i,   0.60000 
                };
                vector<std::complex<double>> m_b = {              
                    25.9622476651464 +  0.0000000000000i,   17.7794485121929 + 13.1574958765530i,   
                    11.2649352514491 + 16.4823940873968i,
                    17.7794485121928 - 13.1574958765530i,   20.5657808536051 -  0.0000000000000i,   
                    17.9925255171787 +  6.0065935802308i,
                    11.2649352514491 - 16.4823940873968i,   17.9925255171787 -  6.0065935802308i,   
                    17.9429273455619 -  0.0000000000000i
                };

                multi_array<std::complex<double>,2> A(extents[3][3]);
                multi_array<std::complex<double>,2> B(extents[3][3]);

                for(int i=0, k=0; i<A.shape()[0]; i++)
                    for(int j=0; j<A.shape()[1]; j++,k++)
                        A[i][j] = m_a[k];
                for(int i=0, k=0; i<A.shape()[0]; i++)
                    for(int j=0; j<A.shape()[1]; j++,k++)
                        B[i][j] = m_b[k];

                multi_array<std::complex<double>,2> C = ma::exp(A);
                verify_approx( C, B );
        }
}

TEST_CASE("dense_ma_operations", "[matrix_operations]")
{
  test_dense_matrix_mult(); 
}

}

