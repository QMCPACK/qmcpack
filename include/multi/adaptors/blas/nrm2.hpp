// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_NRM2_HPP
#define MULTI_ADAPTORS_BLAS_NRM2_HPP

#include "../blas/core.hpp"

#include "../../array.hpp"

#include<complex> // std::norm

namespace boost::multi::blas {

using core::nrm2;

using multi::base; 
using std::norm;  // nvcc11 needs using std::FUNCTION and the FUNCTION (and it works in clang, gcc, culang, icc)

template<class A1D, class A0D>
auto nrm2(A1D const& x, A0D&& res)  // NOLINT(readability-identifier-length) conventional BLAS naming
->decltype(nrm2(size(x), x.base(), x.stride(), base(res)), std::forward<A0D>(res)) {
	return nrm2(size(x), x.base(), x.stride(), base(res)), std::forward<A0D>(res); }

#if 0
template<class A1D>
auto nrm2(A1D const& x, double& r)
->decltype(nrm2(x.size(), x.base(), x.stride(), &r), r){
	return nrm2(x.size(), x.base(), x.stride(), &r), r;}

template<class A1D>
auto nrm2(A1D const& x, float& r)
->decltype(nrm2(x.size(), x.base(), x.stride(), &r), r){
	return nrm2(x.size(), x.base(), x.stride(), &r), r;}
#endif

template<
	class A1D, typename T = double, //  decltype(norm(std::declval<typename A1D::value_type>())), 
	class Alloc = typename std::allocator_traits<typename A1D::default_allocator_type>::template rebind_alloc<T>
>
NODISCARD("")
auto nrm2(A1D const& array)
//->std::decay_t<decltype(nrm2(array, multi::static_array<T, 0, Alloc>({}, x.get_allocator()) ))>{
->std::decay_t<decltype(nrm2(array, multi::static_array<T, 0, Alloc>({})))> { // array.get_allocator() in decltype doesn't work for icc
	return nrm2(array, multi::static_array<T, 0, Alloc>({}, array.get_allocator()));}

template<class Alloc, class A1D, typename T = decltype(norm(std::declval<typename A1D::value_type>())), 
	class AllocR = typename std::allocator_traits<typename A1D::default_allocator_type>::template rebind_alloc<T>
>
NODISCARD("")
auto nrm2(A1D const& array, AllocR const& alloc)
->std::decay_t<decltype(blas::nrm2(array, multi::static_array<T, 0, AllocR>({}, alloc)))> {
	return              blas::nrm2(array, multi::static_array<T, 0, AllocR>({}, alloc)) ; }

namespace operators {
	using std::norm;
	template<class A1D, class Real = decltype(norm(std::declval<typename A1D::value_type>()))>//decltype(norm(std::declval<typename A1D::value_type>()))> 
	NODISCARD("") auto operator^(A1D const& array, int n)
	->decltype(std::pow(Real{blas::nrm2(array)}, n)) {
		return std::pow(Real{blas::nrm2(array)}, n); }
} // end namespace operators

} // end namespace boost::multi::blas

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS nrm2"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>


//#include "../../array.hpp"
//#include "../../complex.hpp"

////#include<thrust/complex.h>

//#include<boost/mpl/list.hpp>

//namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_real){
//	namespace blas = multi::blas;
//	multi::array<double, 2> const cA = {
//		{1.,  2.,  3.,  4.},
//		{5.,  6.,  7.,  8.},
//		{9., 10., 11., 12.}
//	};

//	double n;
//	BOOST_REQUIRE( blas::nrm2(rotated(cA)[1], n) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//	BOOST_REQUIRE( n == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//	BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//	double n2 = blas::nrm2(rotated(cA)[1]);
//	BOOST_REQUIRE( n == n2 );

//	multi::array<double, 1> R(4);
//	blas::nrm2( rotated(cA)[1], R[2]);
//	BOOST_REQUIRE( R[2] ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//	multi::array<double, 0> R0;
//	blas::nrm2( rotated(cA)[1], R0);
//	BOOST_REQUIRE( R0 ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//	BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_blas_nrm2_operators){
//	multi::array<double, 1> X = {1.1,2.1,3.1, 4.1};
//	double n; multi::blas::nrm2(X, n);
//	BOOST_REQUIRE( n == multi::blas::nrm2(X) );

//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case){
//	using complex = std::complex<double>;
//	multi::array<complex, 2> const cA = {
//		{1.,  2.,  3.,  4.},
//		{5.,  6.,  7.,  8.},
//		{9., 10., 11., 12.}
//	};

//	using multi::blas::nrm2;
//	double n; 
//	BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//	BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

//#if 0
//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_thrust){
//	using complex = thrust::complex<double>;
//	multi::array<complex, 2> const cA = {
//		{1.,  2.,  3.,  4.},
//		{5.,  6.,  7.,  8.},
//		{9., 10., 11., 12.}
//	};

//	using multi::blas::nrm2;
//	double n;
//	BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//	BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_types){
//	boost::mpl::for_each<boost::mpl::list<
//		std   ::complex<double>, 
//		thrust::complex<double>//,
//	//	boost::multi::complex<double> // TODO make this work
//	>>([](auto cplx){
//		multi::array<decltype(cplx), 2> const cA = {
//			{1.,  2.,  3.,  4.},
//			{5.,  6.,  7.,  8.},
//			{9., 10., 11., 12.}
//		};

//		using multi::blas::nrm2;
//		double n;
//		BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//		BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//	});
//}
//#endif

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex){
//	using complex = std::complex<double>; complex const I{0,1};
//	multi::array<complex, 2> const cA = {
//		{1.,  2. + 1.*I,  3.,  4.},
//		{5.,  6. + 4.*I,  7.,  8.},
//		{9., 10. - 3.*I, 11., 12.}
//	};

//	using multi::blas::nrm2;
//	double n;
//	BOOST_REQUIRE( nrm2(rotated(cA)[1], n)   == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
//	BOOST_REQUIRE( nrm2(rotated(cA)[1])      == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );

//	using namespace multi::blas::operators;
//	BOOST_TEST_REQUIRE( (rotated(cA)[1]^-1) == 1/std::sqrt(norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1])) , boost::test_tools::tolerance(1e-15) );
//	BOOST_TEST_REQUIRE( (rotated(cA)[1]^2) == norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) , boost::test_tools::tolerance(1e-15) );
//}

//#endif
#endif
