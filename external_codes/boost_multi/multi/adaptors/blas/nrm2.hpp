#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_NRM2_HPP
#define MULTI_ADAPTORS_BLAS_NRM2_HPP

#include "../blas/core.hpp"

#include "../../array.hpp"

#include<complex> // std::norm

namespace boost{
namespace multi{
namespace blas{

using core::nrm2;

template<class X1D, class R0D>
auto nrm2(X1D const& x, R0D&& r)
->decltype(nrm2(size(x), base(x), stride(x), base(r)), std::forward<R0D>(r)){
	return nrm2(size(x), base(x), stride(x), base(r)), std::forward<R0D>(r);}

using std::norm; // for some reason nvcc needs using std::norm/norm (and works in clang, gcc, culang, icc)

template<class X1D, 
	typename T = decltype(norm(std::declval<typename X1D::value_type>())),
	typename Alloc = typename std::allocator_traits<decltype(get_allocator(std::declval<X1D>()))>::template rebind_alloc<T>
>
auto nrm2(X1D const& x){
	return nrm2(x, multi::static_array<T, 0, Alloc>{}); // TODO: this supports only default constructible (deduced) allocator
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_NRM2

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS nrm2"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_real){
	namespace blas = multi::blas;
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};

	double n;
	BOOST_REQUIRE( blas::nrm2(rotated(cA)[1], n) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
	BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

	multi::array<double, 1> R(4);
	blas::nrm2( rotated(cA)[1], R[2]);
	BOOST_REQUIRE( R[2] ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

}

using complex = std::complex<double>; complex const I{0,1};

BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case){
	multi::array<complex, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};

	using multi::blas::nrm2;
	double n; 
	BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

	BOOST_REQUIRE( nrm2(rotated(cA)[1]) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
}

BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex){
	multi::array<complex, 2> const cA = {
		{1.,  2. + 1.*I,  3.,  4.},
		{5.,  6. + 4.*I,  7.,  8.},
		{9., 10. - 3.*I, 11., 12.}
	};

	using multi::blas::nrm2;
	double n;
	BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
	BOOST_REQUIRE( nrm2(rotated(cA)[1]) == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
}

#endif
#endif

