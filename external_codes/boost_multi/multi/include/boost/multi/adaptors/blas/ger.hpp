// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2021 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_GER_HPP
#define MULTI_ADAPTORS_BLAS_GER_HPP

#include "../blas/core.hpp"

namespace boost::multi::blas {

using core::ger;

template<class T, class It1, class Size1, class It2, class Size2, class Out>
auto ger_n(T alpha, It1 x_first, Size1 x_n, It2 y_first, Size2 y_n, Out A_first) -> Out {
	assert( A_first->size() == x_n );
	assert( A_first->stride() == 1 );
	ger(x_n, y_n, alpha, base(x_first), stride(x_first), base(y_first), stride(y_first), base(A_first), stride(A_first));
	return A_first + y_n;
}

template<class T, class It1, class It2, class Out>
auto ger(T alpha, It1 x_first, It1 x_last, It2 y_first, It2 y_last, Out A_first) -> Out {
	assert( stride(x_first) == stride(x_last) );
	assert( stride(y_first) == stride(y_last) );
	return ger_n(alpha, x_first, std::distance(x_first, x_last), y_first, std::distance(y_first, y_last), A_first);
}

template<class T, class X1D, class Y1D, class A2D>
auto ger(T alpha, X1D const& x, Y1D const& y, A2D&& A) -> A2D&& {
	if(stride(A) == 1) {
		auto e = ger(alpha, begin(y), end(y), begin(x), end(x), begin(rotated(A)));
		assert( end(rotated(A)) == e );
	} else {
		assert( size(A) == size(y) );
		auto e = ger(alpha, begin(x), end(x), begin(y), end(y), begin(A));
		assert( end(A) == e );
	}
	return std::forward<A2D>(A);
}

template<class T, class It1, class Size1, class It2, class Size2, class Out>
auto gerc_n(T alpha, It1 x_first, Size1 x_n, It2 y_first, Size2 y_n, Out A_first) -> Out {
	assert( A_first->size() == x_n );
	assert( A_first->stride() == 1 );
	gerc(x_n, y_n, alpha, base(x_first), stride(x_first), base(y_first), stride(y_first), base(A_first), stride(A_first));
	return A_first + y_n;
}

template<class T, class It1, class It2, class Out>
auto gerc(T alpha, It1 x_first, It1 x_last, It2 y_first, It2 y_last, Out A_first) -> Out {
	assert( stride(x_first) == stride(x_last) );
	assert( stride(y_first) == stride(y_last) );
	return gerc_n(alpha, x_first, std::distance(x_first, x_last), y_first, std::distance(y_first, y_last), A_first);
}

template<class T, class X1D, class Y1D, class A2D>
auto gerc(T alpha, X1D const& x, Y1D const& y, A2D&& A) -> A2D {
	if(stride(A) == 1) {
		auto e = gerc(alpha, begin(y), end(y), begin(x), end(x), begin(rotated(A)));
		assert( end(rotated(A)) == e );
	} else {
		assert( size(A) == size(y) );
		auto e = gerc(alpha, begin(x), end(x), begin(y), end(y), begin(A));
		assert( end(A) == e );
	}
	return A;
}

template<class T, class It1, class Size1, class It2, class Size2, class Out>
auto geru_n(T alpha, It1 x_first, Size1 x_n, It2 y_first, Size2 y_n, Out A_first) -> Out {
	assert( A_first->size() == x_n );
	assert( A_first->stride() == 1 );
	geru(x_n, y_n, alpha, base(x_first), stride(x_first), base(y_first), stride(y_first), base(A_first), stride(A_first));
	return A_first + y_n;
}

template<class T, class It1, class It2, class Out>
auto geru(T alpha, It1 x_first, It1 x_last, It2 y_first, It2 y_last, Out A_first) -> Out {
	assert( stride(x_first) == stride(x_last) );
	assert( stride(y_first) == stride(y_last) );
	return geru_n(alpha, x_first, std::distance(x_first, x_last), y_first, std::distance(y_first, y_last), A_first);
}

template<class T, class X1D, class Y1D, class A2D>
auto geru(T alpha, X1D const& x, Y1D const& y, A2D&& A) -> A2D {
	if(stride(A) == 1) {
		auto e = geru(alpha, begin(y), end(y), begin(x), end(x), begin(rotated(A)));
		assert( end(rotated(A)) == e );
	} else {
		assert( size(A) == size(y) );
		auto e = geru(alpha, begin(x), end(x), begin(y), end(y), begin(A));
		assert( end(A) == e );
	}
	return A;
}

} // end namespace boost::multi::blas

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi blas ger"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../../array.hpp"
//#include "../../utility.hpp"

//#include<complex>
//#include<cassert>
//#include<iostream>
//#include<numeric>
//#include<algorithm>

//using std::cout;
//namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(multi_blas_ger){
//	namespace blas = multi::blas;
//	{
//		multi::array<double, 2> A = {
//			{0., 0. ,0.},
//			{0., 0., 0.}
//		};
//		multi::array<double, 1> const x = { 0., 0., 1.};
//		multi::array<double, 1> const y = { 0., 1.};
//		blas::ger(1., x, y, A); // A = a*A + (y^T)(x)
//		for(int i = 0; i != size(A); ++i){
//			for(int j = 0; j != size(A[i]); ++j)
//				std::cout << A[i][j] << ' ';
//			std::cout << std::endl;
//		}
//		std::cout << std::endl;
//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//	//	assert( A[1][1] == 6. );
//	}
//	{
//		multi::array<double, 2> A = {
//			{0., 0.},
//			{0., 0.},
//			{0., 0.}
//		};
//		multi::array<double, 1> const x = {0., 0., 1.};
//		multi::array<double, 1> const y = {0., 1.};
//		blas::ger(1., x, y, rotated(A)); // A^T = a*A^T + (y^T)(x) and A = a*A + (x^T)y
//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//		for(int i = 0; i != size(A); ++i){
//			for(int j = 0; j != size(A[i]); ++j)
//				std::cout << A[i][j] << ' ';
//			std::cout << std::endl;
//		}
////		std::cout << A[1][2] << std::endl;
////		assert( A[1][2] == 1. );
//	}
//	{
////		multi::array<double, 2> A = {
////			{2., 3., 6., 8.},
////			{4., 1., 6., 8.},
////			{0., 1., 6., 8.}
////		};
////		assert( A[1][2] == 6. );
////		multi::array<double, 1> const x = { 0., 1., 0.};
////		multi::array<double, 1> const y = { 0., 0., 1., 0.};
//		
//	//	multi::blas::ger(0., x, y, rotated(A)); // 

//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//	//	assert( A[1][1] == 4. );
//	}
//	{
//		multi::array<double, 2> A = {
//			{2., 3., 6., 8.},
//			{4., 1., 6., 8.},
//			{0., 1., 6., 8.}
//		};
//		multi::array<double, 1> const x = { 1., 2., 5.};
//		multi::array<double, 1> const y = {-2., 1., 1., 1.};
//		blas::ger(1., x, y, A); // 
//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//	//	assert( A[1][1] == 4. );
//	}
//	{
//		multi::array<double, 2> a = {
//			{2., 1., 1.},
//			{3., 4., 0.}
//		};
//		multi::array<double, 1> const x = { 1., 2., 5.};
//		multi::array<double, 1> const y = {-2., 1.};
//		blas::ger(1., x, y, rotated(a));
//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//		assert( a[1][1] == 6. );
//	}
//#if 0
//	{
//		multi::array<std::complex<double>, 2> a = {
//			{2., 3.}, 
//			{1., 4.}, 
//			{1.,0.}
//		};
//		multi::array<std::complex<double>, 1> const x = { 1., 2., 5.};
//		multi::array<std::complex<double>, 1> const y = {-2., 1.};
//		multi::blas::gerc(1., x, y, a);
//	//	a = {{2., 3.}, {1., 4.}, {1., 0.}}; GER[1, {1., 2., 5.}, {-2., 1.}, a]; Print[a] : {{0., 4.}, {-3., 6.}, {-9., 5.}}
//		assert( a[1][1] == 6. );
//	}
//	{
//		multi::array<std::complex<double>, 2> a = {{2. + 1.*I, 3. + 4.*I}, {1.+3.*I, 4. + 2.*I}, {1. + 7.*I, 0.}};
//		multi::array<std::complex<double>, 1> const x = { 1. + 1.*I, 2. + I*9., 5. + 4.*I};
//		multi::array<std::complex<double>, 1> const y = {-2. + 8.*I, 1. + 1.*I};
//		multi::blas::geru(1. + 2.*I, x, y, a); // a = alpha*outer(x, y) + a
////		a = {{2. + 1.*I, 3. + 4.*I}, {1. + 3.*I, 4. + 2.*I}, {1. + 7.*I, 0.}}; GER[1 + 2.*I, {1. + 1.*I, 2. + I*9., 5. + 4.*I}, {-2. + 8.*I,  1. + 1.*I}, a]; Print[a]; 
////		{{-20.-13. I,-1.+6. I},{-71.-151. I,-25.-1. I},{-105.-45. I,-17.+11. I}}
//		std::cout << "a11 " << a[1][1] << std::endl;
//		assert( a[1][1] == -25. - 1.*I );
//	}
//	{
//		multi::array<std::complex<double>, 2> a = {
//			{2. + 1.*I, 1. + 3.*I, 1. + 7.*I},
//			{3. + 4.*I, 4. + 2.*I, 0. + 0.*I}
//		};
//		std::cout << "a = " << size(a) << std::endl;
//		multi::array<std::complex<double>, 1> const x = { 1. + 1.*I, 2. + I*9., 5. + 4.*I};
//		multi::array<std::complex<double>, 1> const y = {-2. + 8.*I, 1. + 1.*I};
//		multi::blas::geru(1. + 2.*I, x, y, rotated(a)); // a = alpha*outer(x, y) + a
////		a = {{2. + 1.*I, 3. + 4.*I}, {1. + 3.*I, 4. + 2.*I}, {1. + 7.*I, 0.}}; GER[1 + 2.*I, {1. + 1.*I, 2. + I*9., 5. + 4.*I}, {-2. + 8.*I,  1. + 1.*I}, a]; Print[a]; 
////		{{-20.-13. I,-1.+6. I},{-71.-151. I,-25.-1. I},{-105.-45. I,-17.+11. I}}
//		std::cout << "here a11 " << a[1][1] << std::endl;
//		assert( a[1][1] == -25. - 1.*I );
//	}
//#endif

//}

//#endif
#endif
