#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -Wall -Wextra -D_TEST_MULTI_ADAPTORS_BLAS $0.cpp -o$0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2018-2020
#ifndef MULTI_ADAPTORS_BLAS_HPP
#define MULTI_ADAPTORS_BLAS_HPP

#include "../adaptors/blas/iamax.hpp"
#include "../adaptors/blas/asum.hpp"
#include "../adaptors/blas/axpy.hpp"
#include "../adaptors/blas/copy.hpp"
#include "../adaptors/blas/dot.hpp"
#include "../adaptors/blas/gemm.hpp"
#include "../adaptors/blas/syrk.hpp"
#include "../adaptors/blas/herk.hpp"
#include "../adaptors/blas/gemv.hpp"
#include "../adaptors/blas/ger.hpp"
#include "../adaptors/blas/nrm2.hpp"
#include "../adaptors/blas/trsm.hpp"
#include "../adaptors/blas/scal.hpp"
#include "../adaptors/blas/swap.hpp"

#if _TEST_MULTI_ADAPTORS_BLAS

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../utility.hpp"

#include<iostream>
#include<complex>
#include<numeric> // iota
#include<algorithm> // transform

namespace multi = boost::multi;
using complex = std::complex<double>;
complex const I{0, 1};

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex){

	using multi::blas::herk;
	{
		multi::array<complex, 2> const A = {
			{1. + 3.*I, 9. + 1.*I}, 
			{3. - 2.*I, 7. - 8.*I}, 
			{4. + 1.*I, 1. - 3.*I}
		};
		multi::array<complex, 2> C({3, 3}, 9999.);
		herk(1., A, C); // herk(A, C); // C†=C=AA†=(A†A)†
		BOOST_REQUIRE( C[1][2] == complex(41., 2.) );
		BOOST_REQUIRE( C[2][1] == conj(C[1][2]) );
	}
}

using T = std::complex<double>;

BOOST_AUTO_TEST_CASE(multi_blas_asum_complex){
	multi::array<T, 1> arr(1000000000, 0.);
//	std::iota(begin(arr), end(arr), -700.);
//	std::transform(cbegin(arr), cend(arr), begin(arr), [](auto&& a){return sqrt(a);});
	{
		using multi::blas::asum;
		BOOST_REQUIRE( asum(arr) == 0 );
	//	std::cout << asum(arr) << std::endl;
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_nrm2_complex){
	multi::array<T, 1> arr(1000000000, 0.);
//	std::iota(begin(arr), end(arr), -700.);
//	std::transform(cbegin(arr), cend(arr), begin(arr), [](auto&& a){return sqrt(a);});
	{
		using multi::blas::nrm2;
		BOOST_REQUIRE( nrm2(arr) == 0. );
	}
}

#if 0
	multi::array<double, 2> const CA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	{
		double const a0 = 2./3.; 
		double const b0 = 4./5.;
		double a = a0, b = b0;
		double c, s;
		using multi::blas::rotg;
		rotg(a, b, c, s);
		using std::abs; using std::sqrt;
		assert( abs(c - 5./sqrt(61.)) < 1e-15 );
		assert( abs(s - 6./sqrt(61.)) < 1e-15 );
		assert( abs(a - (b0>0?1:-1)*sqrt(a0*a0 + b0*b0)) < 1e-15 );
		assert( abs(  c*c  + s*s  - 1 ) < 1e-15 );
		assert( abs(  c*a0 + s*b0 - a ) < 1e-15 );
		assert( abs( -s*a0 + c*b0     ) < 1e-15 );
	}
	{
		using multi::blas::rotmg;
		double const x0 = 2./3.; 
		double const y0 = 4./5.;
		double const D1 = 1.;
		double const D2 = 1.;
		{
			double d1 = D1;
			double d2 = D2;
			double x1 = x0;
			double const y1 = y0;
			auto m = rotmg(d1, d2, x1, y1);
			assert( std::abs(x1 -( m.h()[0][0]*x0*std::sqrt(D1) + m.h()[0][1]*y0*std::sqrt(D2) )) < 1e-15 );
			assert( std::abs(      m.h()[1][0]*x0*std::sqrt(D1) + m.h()[1][1]*y0*std::sqrt(D2)  ) < 1e-15 );
		}
		{
			double x1 = x0;
			double const y1 = y0;
			double d1 = D1;
			double d2 = D2;
			multi::array<double, 1> X0 = {x0*std::sqrt(D1)};
			multi::array<double, 1> Y0 = {y0*std::sqrt(D2)};
			multi::array<double, 1> X1 = X0;
			multi::array<double, 1> Y1 = Y0;
			rotm(X1, Y1, rotmg(d1, d2, x1, y1));
			assert( std::abs( X1[0] - x1 ) <1e-15 );
			assert( Y1[0] == 0. );
		}
	}
	{
		multi::array<double, 1> X = CA[0];
		multi::array<double, 1> Y = CA[2];
		using multi::blas::rot;
		using std::cos; using std::sin;
		rot(X, Y, cos( 1.2), sin( 1.2));
		assert(X[1] == CA[0][1]*cos(1.2) + CA[2][1]*sin(1.2));
		assert(Y[1] == CA[2][1]*cos(1.2) - CA[0][1]*sin(1.2));
	}
	{
		multi::array<double, 1> const a0 = {2./3.};
		multi::array<double, 1> const b0 = {4./5.};
		using multi::blas::rotg;
		{
			double a = a0[0], b = b0[0];
			auto cs = rotg(a, b);
			multi::array<double, 1> a1 = a0;
			multi::array<double, 1> b1 = b0;
			rot(a1, b1, cs);
			assert( std::abs(a1[0] - a) < 1e-15 );
			assert( std::abs(b1[0]    ) < 1e-15 );
		}
		{
			double a = a0[0], b = b0[0];
			multi::array<double, 1> a1 = a0;
			multi::array<double, 1> b1 = b0;
			rot(a1, b1, rotg(a, b));
			assert( std::abs(a1[0] - a) < 1e-15 );
			assert( std::abs(b1[0]    ) < 1e-15 );
		}
	}
	{
		using multi::blas::dot;
		auto d = dot(CA[1], CA[2]);
		assert(d == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.));
	}
	using dcomplex = std::complex<double>;
	{
		multi::array<dcomplex, 2> A = CA;
		A[1][1] += dcomplex{1.1, 2.1};
		A[2][1] -= dcomplex{1.1, 2.1};
		using multi::blas::dotu;
		using multi::blas::dotc;
		using multi::blas::nrm2;
		using multi::blas::asum;
		assert(dotu(A[1], A[2]) == std::inner_product(begin(A[1]), begin(A[2]), end(A[1]), dcomplex{}, std::plus<>{}, [](auto&& a, auto&& b){return a*b;}));
		assert(dotc(A[1], A[2]) == std::inner_product(begin(A[1]), begin(A[2]), end(A[1]), dcomplex{}, std::plus<>{}, [](auto&& a, auto&& b){return conj(a)*b;}));
		assert(nrm2(A[1]) == std::sqrt(dotc(A[1], A[1])));
		assert(dotu(A[1], A[2]) == std::inner_product(begin(A[1]), begin(A[2]), end(A[1]), dcomplex{}, std::plus<>{}, [](auto&& a, auto&& b){return a*b;}));
		assert(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a + std::abs(real(b)) + std::abs(imag(b));}));
	}
	{
		auto const& A = CA.rotated(1)[1]; (void)A;
		using multi::blas::iamax;
		assert(iamax(A) == std::distance(begin(A), std::max_element(begin(A), end(A), [](auto&& a, auto&& b){
			return std::abs(a) < std::abs(b);
		})));
	}
	
///////////////////////////////////////////////////////////////////////////////
	{
		multi::array<double, 2> const M = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		multi::array<double, 1> const X = {1.1,2.1,3.1, 4.1};
		multi::array<double, 1> Y = {4.,5.,6.};
		multi::array<double, 1> Y2 = Y;
		multi::array<double, 1> Y3 = {214.02, 106.43, 188.37};
		double a = 1.1, b = 1.2;
		multi::blas::gemv(a, M, X, b, Y);
		multi::blas::gemv<double, double>(a, M, X, b, Y2);
		assert( Y == Y2 );
		assert( std::abs(Y[1] - Y3[1]) < 1e-14 );
	}
	{
		multi::array<double, 2> const M = {
			{ 9., 24., 30., 9.},
			{ 4., 10., 12., 7.},
			{14., 16., 36., 1.}
		};
		assert( M[2][0] == 14. );
		multi::array<double, 1> const X = {1.1,2.1,3.1};
		multi::array<double, 1> Y = {4.,5.,6., 7.};
		multi::array<double, 1> Y2 = Y;
		multi::array<double, 1> Y3 = {72.67, 112.7, 193.98, 38.87};
		double a = 1.1, b = 1.2;
		multi::blas::gemv(a, M.rotated(1), X, b, Y);
		multi::blas::gemv<double, double>(a, M.rotated(1), X, b, Y2);
		assert( std::abs(Y[1] - Y2[1]) < 1e-13 );
		assert( std::abs(Y[1] - Y3[1]) < 1e-13 );
	}
	auto const I = dcomplex{0.,1.};
	{
		multi::array<dcomplex, 2> const M = {
			{ 9. + 1.*I, 24. + 2.*I, 30. + 3.*I, 9. + 1.*I}, 
			{ 4. + 1.*I, 10. + 1.*I, 12. - 2.*I, 7. + 2.*I}, 
			{14. + 3.*I, 16. - 4.*I, 36. + 1.*I, 1. - 2.*I}
		};
		multi::array<dcomplex, 1> const X = {1.1+I*2., 2.1+I*1.1, 3.1+I*8. , 4.1+I*1.2};
		multi::array<dcomplex, 1> Y = {4.+I*3.1,5.-I*9.,6.+I*1.};
		multi::array<dcomplex, 1> Y2 = Y;
		multi::array<dcomplex, 1> const Y3 = {-486.81+698.69*I, -125.08+359.44*I, -504.21+707.01*I};
		dcomplex a = 1.1+I*2.1, b = 1.2+I*3.;
		cout<<">>"<<__LINE__ <<std::endl;
		multi::blas::gemv(a, M, X, b, Y);
		cout<<">>"<<__LINE__ <<std::endl;
		multi::blas::gemv<dcomplex, dcomplex>(a, M, X, b, Y2);
		using std::abs;
		assert( abs(Y[0] - Y3[0]) < 1e-12 && abs(Y[1] - Y3[1]) < 1e-12 && abs(Y[2] - Y3[2]) < 1e-12 );
	}
	{
		multi::array<dcomplex, 2> const M = {
			{9. + 1.*I, 4. + 1.*I, 14. + 3.*I}, 
			{24. + 2.*I, 10. + 1.*I, 16. - 4.*I}, 
			{30. + 3.*I, 12. - 2.*I, 36. + 1.*I}, 
			{9. + 1.*I,   7. + 2.*I, 1. - 2.*I}
		};
		multi::array<dcomplex, 1> const X = {1.1+I*2., 2.1+I*1.1, 3.1+I*8. , 4.1+I*1.2};
		multi::array<dcomplex, 1> Y = {4.+I*3.1,5.-I*9.,6.+I*1.};
		multi::array<dcomplex, 1> Y2 = Y;
		multi::array<dcomplex, 1> const Y3 = {-486.81+698.69*I, -125.08+359.44*I, -504.21+707.01*I};
		std::complex<double> a = 1.1+I*2.1, b = 1.2+I*3.;
		cout<<">>"<<__LINE__ <<std::endl;
		multi::blas::gemv(a, M.rotated(), X, b, Y);
		cout<<">>"<<__LINE__ <<std::endl;
		multi::blas::gemv<dcomplex, dcomplex>(a, M.rotated(), X, b, Y2);
		assert( abs(Y[0] - Y3[0]) < 1e-12 && abs(Y[1] - Y3[1]) < 1e-12 && abs(Y[2] - Y3[2]) < 1e-12 );
		assert( abs(Y[0] - Y2[0]) < 1e-12 && abs(Y[1] - Y2[1]) < 1e-12 && abs(Y[2] - Y2[2]) < 1e-12 );
	}
#if 0
	{
		multi::array<dcomplex, 2> const M = {
			{9. + 1.*I, 4. + 1.*I, 14. + 3.*I}, 
			{24. + 2.*I, 10. + 1.*I, 16. - 4.*I}, 
			{30. + 3.*I, 12. - 2.*I, 36. + 1.*I}, 
			{9. + 1.*I,   7. + 2.*I, 1. - 2.*I}
		};
		multi::array<dcomplex, 1> const X = {1.1+I*2., 2.1+I*1.1, 3.1+I*8. , 4.1+I*1.2};
		multi::array<dcomplex, 1> Y = {4.+I*3.1,5.-I*9.,6.+I*1.};
	//	multi::array<dcomplex, 1> Y2 = Y;
	//	multi::array<dcomplex, 1> const Y3 = {-486.81+698.69*I, -125.08+359.44*I, -504.21+707.01*I};
		std::complex<double> a = 1.1+I*2.1, b = 1.2+I*3.;
		cout<<">>"<<__LINE__ <<std::endl;
		multi::blas::gemv(a, M.rotated(), X, b, Y, multi::blas::conj{});
	//	cout<<">>"<<__LINE__ <<std::endl;
	//	multi::blas::gemv<dcomplex, dcomplex>(a, M.rotated(), X, b, Y2);
		cout<< Y[0] <<' '<< Y[1] <<' '<< Y[2] <<std::endl;
	//	assert( abs(Y[0] - Y3[0]) < 1e-12 && abs(Y[1] - Y3[1]) < 1e-12 && abs(Y[2] - Y3[2]) < 1e-12 );
	//	assert( abs(Y[0] - Y2[0]) < 1e-12 && abs(Y[1] - Y2[1]) < 1e-12 && abs(Y[2] - Y2[2]) < 1e-12 );
	}
#endif
	return 0;
//	return 0;
	{
//		multi::array<dcomplex, 2> const M = {
//			{ 9.+I*1., 24.+I*2., 30.+I*3.},
//			{ 4.+I*1., 10.+I*1., 12.-I*2.},
//			{14.+I*3., 16.-I*4., 36.+I*1.},
//			{ 9.+I*1.,  7.+I*2.,  1.-I*2.}
//		};
		multi::array<dcomplex, 2> const M = {
			{ 9. + 1.*I,  4. + 1.*I, 14. + 3.*I, 9. + 1.*I}, 
			{24. + 2.*I, 10. + 1.*I, 16. - 4.*I, 7. + 2.*I}, 
			{30. + 3.*I, 12. - 2.*I, 36. + 1.*I, 1. - 2.*I}
		};
		multi::array<dcomplex, 1> const X = {1.1+I*2., 2.1+I*1.1, 3.1+I*8., 4.1+I*1.2};
		multi::array<dcomplex, 1> Y = {4.+I*3.1,5.-I*9.,6.+I*1.};
		multi::array<dcomplex, 1> Y2 = Y;
		multi::array<dcomplex, 1> const Y3 = {-134.97+423.67*I, -265.81+431.55*I, -567.81+809.37*I};
		dcomplex const a = 1.1+I*2.1, b = 1.2+I*3.;
		cout<< "708" <<std::endl;
	//	multi::blas::gemv(a, M.rotated(), X, b, Y, multi::blas::conj<>{});
		zgemv_('N', std::get<0>(M.shape()), std::get<1>(M.shape()), a, M.base(), M.stride(), X.base(), X.stride(), b, Y.base(), Y.stride());
	//	zgemv_('T', std::get<1>(M.shape()), std::get<0>(M.shape()), a, M.base(), 2*std::get<0>(M.strides()), X.base(), stride(X), b, Y.base(), stride(Y));
		multi::blas::gemv<std::complex<double>, std::complex<double>>(a, M, X, b, Y2);
	//	multi::blas::gemv<std::complex<double>, std::complex<double>>(a, M.rotated(1), X, b, Y2);

	//	multi::blas::gemv<dcomplex, dcomplex>(a, M.rotated(), X, b, Y2, multi::blas::conj<>{});
		cout << Y[0] <<' '<< Y[1] <<' '<< Y[2] <<std::endl;
		cout << Y2[0] <<' '<< Y2[1] <<' '<< Y2[2] <<std::endl;
		cout << "finished" << std::endl;
	//	assert( std::abs(Y[1] - Y3[1]) < 1e-12 );
	//	assert( Y[1] == Y2[1] );
	}
//	assert(0);
#endif

#if 0

namespace boost{
namespace multi{
namespace blas{

template<class T> struct cs{
	T c; T s;
	operator multi::array<T, 2>() const{return {{c, s}, {-s, c}};}
};
template<class T> struct ab{T a; T b; using value_type = T;};
template<class T> struct modified_rotation{
	T data_[5];
	int flag() const{return data_[0];}
	multi::array<T, 2> h() const{
		switch(flag()){
			case -1: return {{data_[1], data_[2]}, {data_[3], data_[4]}};
			case  0: return {{T{+1}   , data_[2]}, {data_[3], T{+1}   }};
			case  1: return {{data_[1], T{+1}   }, {T{-1}   , data_[4]}};
			case -2: return {{T{+1}   , T{ 0}   }, {T{ 0}   , T{+1}   }};
			default: assert(0); return {};
		}
	}
};

template<class T>
auto rotg(T& a, T& b){
	cs<T> ret;
//	using blas::rotg;
	rotg(a, b, ret.c, ret.s );
	return ret;
}

template<class T>
modified_rotation<T> rotmg(T& d1, T& d2, T& x1, T const& y1){
	modified_rotation<T> ret;
	rotmg(d1, d2, x1, y1, ret.data_);
	return ret;
}

//template<class T>
//auto rotmg(T& d1, T& d2, T& b1, T const& b2){
//	modified_rotation<T> ret;
//	rotmg(d1, d2, b1, b2, ret);
//	return ret;
//}

template<class X1D, class Y1D, class T>
auto rot(X1D&& x, Y1D&& y, T const& c, T const& s){
	assert( size(x) == size(y) );
	assert( offset(x) == 0 and offset(y) == 0 );
//	using blas::rot;
	rot(size(x), origin(x), stride(x), origin(y), stride(y), c, s);
	return std::tie(x, y);
}
template<class X1D, class Y1D, class CS>
auto rot(X1D&& x, Y1D&& y, CS const& cs){
	return rot(std::forward<X1D>(x), std::forward<Y1D>(y), cs.c, cs.s);
}
template<class X1D, class Y1D, class M>
auto rotm(X1D&& x, Y1D&& y, M const& param){
	using boost::multi::size;
	assert( size(x) == size(y) );
	assert( offset(x) == 0 and offset(y) == 0);
	rotm(size(x), origin(x), stride(x), origin(y), stride(y), param.data_);
}

template<class It>
auto iamax(It first, It last){
	assert( stride(first) == stride(last) );
//	using blas::iamax;
	return iamax(std::distance(first, last), base(first), stride(first))
	#ifndef CBLAS_H
		- 1
	#endif
	;
}

template<class X1D> 
auto iamax(X1D const& x){
	assert( not offset(x) );
	return iamax(begin(x), end(x));
}

}}

#endif
#endif
#endif

