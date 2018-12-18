//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Alfredo Correa, correaa@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef MA_BLAS_HPP
#define MA_BLAS_HPP

#include "Numerics/OhmmsBlas.h"
#include<utility> //std::enable_if
#include<cassert>
#include<iostream>

namespace ma{

template<class MultiArray1Dx, 
         class MultiArray1Dy,
         typename = typename std::enable_if<std::decay<MultiArray1Dx>::type::dimensionality == 1>::type,
         typename = typename std::enable_if<std::decay<MultiArray1Dy>::type::dimensionality == 1>::type 
>
// decltype(BLAS::dot(x.size(), x.origin(), x.strides()[0], y.origin(), y.strides()[0]))
auto 
dot(MultiArray1Dx&& x, MultiArray1Dy&& y){
        assert(x.size() == y.size());
        return BLAS::dot(x.size(), x.origin(), x.strides()[0], y.origin(), y.strides()[0]);
}

template<class T, class MultiArray1D, typename = typename std::enable_if<std::decay<MultiArray1D>::type::dimensionality == 1>::type >
MultiArray1D scal(T a, MultiArray1D&& x){
	BLAS::scal(x.size(), a, x.origin(), x.strides()[0]);
	return std::forward<MultiArray1D>(x);
}

template<class T, class MultiArray1D>
auto operator*=(MultiArray1D&& x, T a) -> decltype(scal(a, std::forward<MultiArray1D>(x))){
	return scal(a, std::forward<MultiArray1D>(x));
}

template<class T, class MultiArray1DA, class MultiArray1DB, 
	typename = typename std::enable_if<std::decay<MultiArray1DA>::type::dimensionality == 1 and std::decay<MultiArray1DB>::type::dimensionality == 1>::type
>
MultiArray1DB axpy(T x, MultiArray1DA const& a, MultiArray1DB&& b){
	assert( a.shape()[0] == b.shape()[0] );
	BLAS::axpy(a.shape()[0], x, a.origin(), a.strides()[0], b.origin(), b.strides()[0]);
	return std::forward<MultiArray1DB>(b);
}

template<char IN, class T, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY,
	typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and MultiArray1DX::dimensionality == 1 and std::decay<MultiArray1DY>::type::dimensionality == 1>::type
>
MultiArray1DY gemv(T alpha, MultiArray2DA const& A, MultiArray1DX const& x, T beta, MultiArray1DY&& y){
        assert( (IN == 'N') || (IN == 'T') || (IN == 'C')  );
	if(IN == 'T' or IN == 'C') assert( x.shape()[0] == A.shape()[1] and y.shape()[0] == A.shape()[0]);
	else if(IN == 'N') assert( x.shape()[0] == A.shape()[0] and y.shape()[0] == A.shape()[1]);
	assert( A.strides()[1] == 1 ); // gemv is not implemented for arrays with non-leading stride != 1
	int M = A.shape()[1];
	int N = A.shape()[0];
	BLAS::gemv(IN, M, N, alpha, A.origin(), A.strides()[0], x.origin(), x.strides()[0], beta, y.origin(), y.strides()[0]);
	return std::forward<MultiArray1DY>(y);
} //y := alpha*A*x + beta*y,

template<char IN, class MultiArray2DA, class MultiArray1DX, class MultiArray1DY>
MultiArray1DY gemv(MultiArray2DA const& A, MultiArray1DX const& x, MultiArray1DY&& y){
	return gemv<IN>(1., A, x, 0., std::forward<MultiArray1DY>(y));
} //y := alpha*A*x

//	gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
//	gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
//	gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)
//	gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC, 
	typename = typename std::enable_if< MultiArray2DA::dimensionality == 2 and MultiArray2DB::dimensionality == 2 and std::decay<MultiArray2DC>::type::dimensionality == 2>::type
>
MultiArray2DC gemm(T alpha, MultiArray2DA const& a, MultiArray2DB const& b, T beta, MultiArray2DC&& c){
	assert( a.strides()[1] == 1 );
	assert( b.strides()[1] == 1 );
	assert( c.strides()[1] == 1 );
	assert( (TA == 'N') || (TA == 'T') || (TA == 'C')  );
	assert( (TB == 'N') || (TB == 'T') || (TB == 'C')  );
	int M = -1;
	int N = -1;
	int K = -1;
	if(TA == 'N' and TB == 'N'){
		M = a.shape()[1];
		N = b.shape()[0];
		K = a.shape()[0];
		assert(a.shape()[0] == b.shape()[1] and c.shape()[0] == b.shape()[0] and c.shape()[1] == a.shape()[1]);
	}
	if((TA == 'T' or TA == 'C') and (TB == 'T' or TB == 'C')){
		M = a.shape()[0];
		N = b.shape()[1];
		K = a.shape()[1];
		assert(a.shape()[1] == b.shape()[0] and c.shape()[0] == b.shape()[1] and c.shape()[1] == a.shape()[0]);
	}
	if((TA == 'T' or TA == 'C') and TB == 'N'){
		M = a.shape()[0];
		N = b.shape()[0];
		K = a.shape()[1];
		assert(a.shape()[1] == b.shape()[1] and c.shape()[0] == b.shape()[0] and c.shape()[1] == a.shape()[0]);
	}
	if(TA == 'N' and (TB == 'T' or TB == 'C')){
		M = a.shape()[1];
		N = b.shape()[1];
		K = a.shape()[0];
		assert(a.shape()[0] == b.shape()[0] and c.shape()[0] == b.shape()[1] and c.shape()[1] == a.shape()[1]);
	}
	BLAS::gemm(
		TA, TB, 
		M, N, K, alpha, 
		a.origin(), a.strides()[0], 
		b.origin(), b.strides()[0],
		beta, 
		c.origin(), c.strides()[0]
	);
	return std::forward<MultiArray2DC>(c);
}

template<char TA, char TB, class T, class MultiArray2DA, class MultiArray2DB, class MultiArray2DC>
MultiArray2DC gemm(MultiArray2DA const& a, MultiArray2DB const& b, MultiArray2DC&& c){
	return gemm(1., a, b, 0., std::forward<MultiArray2DC>(c));
}

}

#ifdef _TEST_MA_BLAS

#include<boost/multi_array.hpp>
#include<iostream>

using std::cout;

int main(){

	std::vector<double> v = {1.,2.,3.};
	{
		boost::multi_array_ref<double, 1> V(v.data(), boost::extents[v.size()]);
		ma::scal(2., V);
		{
			std::vector<double> v2 = {2.,4.,6.};
			boost::multi_array_ref<double, 1> V2(v2.data(), boost::extents[v2.size()]);
			assert( V == V2 );
		}
	}
	
	std::vector<double> m = {
		1.,2.,3.,
		4.,5.,6.,
		7.,8.,9.
	};
	boost::multi_array_ref<double, 2> M(m.data(), boost::extents[3][3]);
	assert( M.num_elements() == m.size());
	ma::scal(2., M[2]);
	{
		std::vector<double> m2 = {
			1.,2.,3.,
			4.,5.,6.,
			14.,16.,18.
		};
		boost::multi_array_ref<double, 2> M2(m2.data(), boost::extents[3][3]);
		assert( M == M2 );
	}
	typedef boost::multi_array_types::index_range range_t;
	boost::multi_array_types::index_gen indices;

	ma::scal(2., M[ indices[range_t(0,3)][2] ]);
	{
		std::vector<double> m2 = {
			1.,2.,6.,
			4.,5.,12.,
			14.,16.,36.
		};
		boost::multi_array_ref<double, 2> M2(m2.data(), boost::extents[3][3]);
		assert( M == M2 );
	}
	ma::scal(2., M[ indices[range_t(0,2)][1] ]);
	{
		std::vector<double> m2 = {
			1.,4.,6.,
			4.,10.,12.,
			14.,16.,36.
		};
		boost::multi_array_ref<double, 2> M2(m2.data(), boost::extents[3][3]);
		assert( M == M2 );
	}
	ma::axpy(2., M[1], M[0]); // M[0] += a*M[1]
	{
		std::vector<double> m2 = {
			9.,24.,30.,
			4.,10.,12.,
			14.,16.,36.
		};
		boost::multi_array_ref<double, 2> M2(m2.data(), boost::extents[3][3]);
		assert( M == M2 );
	}
	{
		std::vector<double> m = {
			9.,24.,30., 9.,
			4.,10.,12., 7.,
			14.,16.,36., 1.
		};
		boost::multi_array_ref<double, 2> M(m.data(), boost::extents[3][4]);
		assert( M[2][0] == 14. );
		std::vector<double> x = {1.,2.,3., 4.};
		boost::multi_array_ref<double, 1> X(x.data(), boost::extents[x.size()]);
		std::vector<double> y = {4.,5.,6.};
		boost::multi_array_ref<double, 1> Y(y.data(), boost::extents[y.size()]);
		ma::gemv<'T'>(1., M, X, 0., Y); // y := M x

		std::vector<double> y2 = {183., 88.,158.};
		boost::multi_array_ref<double, 1> Y2(y2.data(), boost::extents[y2.size()]);
		assert( Y == Y2 );
	}
	{
		std::vector<double> m = {
			9.,24.,30., 9.,
			4.,10.,12., 7.,
			14.,16.,36., 1.
		};
		boost::multi_array_ref<double, 2> M(m.data(), boost::extents[3][4]);
		typedef boost::multi_array_types::index_range range_t;
		boost::multi_array_types::index_gen indices;

		std::vector<double> x = {1.,2.};
		boost::multi_array_ref<double, 1> X(x.data(), boost::extents[x.size()]);
		std::vector<double> y = {4.,5.};
		boost::multi_array_ref<double, 1> Y(y.data(), boost::extents[y.size()]);
		
		auto const& mm = M[ indices[range_t(0,2,1)][range_t(0,2,1)] ];//, X, 0., Y); // y := M x
		ma::gemv<'T'>(1., M[ indices[range_t(0,2,1)][range_t(0,2,1)] ], X, 0., Y); // y := M x

		std::vector<double> y2 = {57., 24.};
		boost::multi_array_ref<double, 1> Y2(y2.data(), boost::extents[y2.size()]);
		assert( Y == Y2 );
	}
	{
		std::vector<double> m = {
			9.,24.,30.,
			4.,10.,12.,
			14.,16.,36.,
			4.,9.,1.
		};
		boost::multi_array_ref<double, 2> M(m.data(), boost::extents[4][3]);
		assert( M[2][0] == 14. );
		std::vector<double> x = {1.,2.,3.};
		boost::multi_array_ref<double, 1> X(x.data(), boost::extents[x.size()]);
		std::vector<double> y = {4.,5.,6., 7.};
		boost::multi_array_ref<double, 1> Y(y.data(), boost::extents[y.size()]);
		ma::gemv<'T'>(1., M, X, 0., Y); // y := M x
		std::vector<double> y2 = {147., 60.,154.,25.};
		boost::multi_array_ref<double, 1> Y2(y2.data(), boost::extents[y2.size()]);
		assert( Y == Y2 );
	}
	{
		std::vector<double> m = {
			9.,24.,30., 9.,
			4.,10.,12., 7.,
			14.,16.,36., 1.
		};
		boost::multi_array_ref<double, 2> M(m.data(), boost::extents[3][4]);
		assert( M[2][0] == 14. );
		std::vector<double> x = {1.,2.,3.};
		boost::multi_array_ref<double, 1> X(x.data(), boost::extents[x.size()]);
		std::vector<double> y = {4.,5.,6.,7.};
		boost::multi_array_ref<double, 1> Y(y.data(), boost::extents[y.size()]);
		ma::gemv<'N'>(1., M, X, 0., Y); // y := M^T x
		std::vector<double> y2 = {59., 92., 162., 26.};
		boost::multi_array_ref<double, 1> Y2(y2.data(), boost::extents[y2.size()]);
		assert( Y == Y2 );
	}
	{
		std::vector<double> a = {
			9.,24.,30., 2.,
			4.,10.,12., 9.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[2][4]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 6., 8., 
			4.,10., 2., 5.,
			14.,16., 9., 0.
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[3][4]);
		assert( B.num_elements() == b.size());

		std::vector<double> c(6);
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[3][2]);
		assert( C.num_elements() == c.size());

		ma::gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)

		std::vector<double> tab = {
			853., 420.,
			346., 185.,
			780., 324.
		};
		boost::multi_array_ref<double, 2> TAB(tab.data(), boost::extents[3][2]);
		assert( TAB.num_elements() == tab.size());

		assert( C == TAB );
	}
	{
		std::vector<double> a = {
			9.,24.,30.,
			4.,10.,12.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[2][3]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 6., 8., 
			4.,10., 2., 5.,
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[2][4]);
		assert( B.num_elements() == b.size());

		std::vector<double> c(12);
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[4][3]);
		assert( C.num_elements() == c.size());


		ma::gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B

		std::vector<double> tab = {
			97., 256., 318.,
			256., 676., 840.,
			62., 164., 204.,
			92., 242., 300.
		};
		boost::multi_array_ref<double, 2> TAB(tab.data(), boost::extents[4][3]);
		assert( TAB.num_elements() == tab.size());

		cout << "A = \n";
		for(int i = 0; i != A.shape()[0]; ++i, cout << '\n')
			for(int j = 0; j != A.shape()[1]; ++j)
				cout << A[i][j] << ' ';
		cout << '\n';
		cout << "B = \n";
		for(int i = 0; i != B.shape()[0]; ++i, cout << '\n')
			for(int j = 0; j != B.shape()[1]; ++j)
				cout << B[i][j] << ' ';
		cout << '\n';		
		cout << "C = \n";
		for(int i = 0; i != C.shape()[0]; ++i, cout << '\n')
			for(int j = 0; j != C.shape()[1]; ++j)
				cout << C[i][j] << ' ';
		cout << '\n';

		assert( C == TAB );
	}
	{
		std::vector<double> a = {
			9.,24.,30.,
			4.,10.,12.,
			3.,11.,45.,
			1.,2., 6.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[4][3]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 6., 8., 
			4.,10., 2., 5.,
			14.,16., 9., 0.
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[3][4]);
		assert( B.num_elements() == b.size());

		std::vector<double> c(9);
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[3][3]);
		assert( C.num_elements() == c.size());

		ma::gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)

		std::vector<double> tab = {
			203., 538., 876.,
			87., 228., 360.,
			217., 595., 1017.
		};
		boost::multi_array_ref<double, 2> TAB(tab.data(), boost::extents[3][3]);
		assert( TAB.num_elements() == tab.size());
		assert( C == TAB );
	}
	{
		std::vector<double> a = {
			9.,24.,30.,
			4.,10.,12.,
			14.,16.,36.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[3][3]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 4.,
			4.,10., 1.,
			14.,16.,3.
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[3][3]);
		assert( B.num_elements() == b.size());
		std::vector<double> c(9);
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[3][3]);
		assert( C.num_elements() == c.size());

		
		ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
		ma::gemm<'T', 'N'>(1., A, B, 0., C); // C = T(A*T(B)) = B*T(A) or T(C) = A*T(B)
		ma::gemm<'N', 'T'>(1., A, B, 0., C); // C =  T(T(A)*B) = T(B)*A or T(C) = T(A)*B
		ma::gemm<'N', 'N'>(1., A, B, 0., C); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
	}
	{
		std::vector<double> a = {
			9.,24.,30.,
			4.,10.,12.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[2][3]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 6., 8., 
			4.,10., 2., 5.,
			14.,16., 9., 0.
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[3][4]);
		assert( B.num_elements() == b.size());

		std::vector<double> c(8);
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[4][2]);
		assert( C.num_elements() == c.size());

		ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B

		std::vector<double> tab = {
			597, 244,
			936, 388,
			372, 152,
			192, 82
		};
		boost::multi_array_ref<double, 2> TAB(tab.data(), boost::extents[4][2]);
		assert( TAB.num_elements() == tab.size());
		assert( C == TAB );
	}
	{
		std::vector<double> a = {
			9.,24.,30., 45.,
			4.,10.,12., 12.
		};
		boost::multi_array_ref<double, 2> A(a.data(), boost::extents[2][4]);
		assert( A.num_elements() == a.size() );
		std::vector<double> b = {
			9.,24., 56.,
			4.,10., 78.,
			14.,16., 90.,
			6., 9., 18.
		};
		boost::multi_array_ref<double, 2> B(b.data(), boost::extents[4][3]);
		assert( B.num_elements() == b.size());
		std::vector<double> c = {
			9.,24., 8.,
			4.,10., 9.
		};
		boost::multi_array_ref<double, 2> C(c.data(), boost::extents[3][2]);
		assert( C.num_elements() == c.size());
		ma::gemm<'T', 'T'>(1., A, B, 0., C); // C = T(A*B) = T(B)*T(A) or T(C) = A*B
	}
	cout << "end test" << std::endl;
}

#endif
#endif

