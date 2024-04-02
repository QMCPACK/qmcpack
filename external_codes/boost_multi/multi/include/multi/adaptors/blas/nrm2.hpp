// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#ifndef MULTI_ADAPTORS_BLAS_NRM2_HPP
#define MULTI_ADAPTORS_BLAS_NRM2_HPP

#include <multi/adaptors/blas/core.hpp>

#include <multi/array.hpp>

#include<complex>  // std::norm

namespace boost::multi::blas {

using core::nrm2;

using multi::base;
using std::norm;  // nvcc11 needs using std::FUNCTION and the FUNCTION (and it works in clang, gcc, culang, icc)

template<class It, class Size, class A0D>
auto nrm2_n(It const& x, Size n, A0D res)  // NOLINT(readability-identifier-length) conventional BLAS naming
//->decltype(blas::default_context_of(x.base())->nrm2(n, x.base(), x.stride(), res), std::next(res)) {  // NOLINT(fuchsia-default-arguments-calls)
{   return blas::default_context_of(x.base())->nrm2(n, x.base(), x.stride(), res), std::next(res); }  // NOLINT(fuchsia-default-arguments-calls)

template<class A1D, class A0D>
auto nrm2(A1D const& x, A0D&& res)  // NOLINT(readability-identifier-length) conventional BLAS naming
//->decltype(nrm2_n(x.begin(), x.size(), &res)) {
{   return nrm2_n(x.begin(), x.size(), &res); }

template<class ItX, class Size>
class nrm2_ptr {
	ItX  x_first_;
	Size count_;

 protected:
	nrm2_ptr(ItX x_first, Size count) : x_first_{x_first}, count_{count} {}

 public:
	explicit operator bool() const {return true;}

	template<class ItOut, class Size2>
	friend constexpr auto copy_n(nrm2_ptr first, Size2 count, ItOut d_first) {
//  ->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {
		assert(count == 1); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }

	template<class ItOut, class Size2>
	friend constexpr auto uninitialized_copy_n(nrm2_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {assert(count == 1);
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }

	template<class ItOut, class Size2>
	static constexpr auto uninitialized_copy_n(nrm2_ptr first, Size2 count, ItOut d_first)
	->decltype(blas::nrm2_n(std::declval<ItX>(), Size2{}     , d_first), d_first + count) {assert(count == 1);
		return blas::nrm2_n(first.x_first_     , first.count_, d_first), d_first + count; }
};

template<class X, class Ptr = nrm2_ptr<typename X::const_iterator, typename X::size_type>>
struct nrm2_ref : private Ptr {
	using decay_type = decltype(norm(std::declval<typename X::value_type>()));
	explicit nrm2_ref(X const& x) : Ptr{begin(x), size(x)} {}  // NOLINT(readability-identifier-length) BLAS naming

	constexpr auto operator&() const& -> Ptr const& {return *this;}  // NOLINT(google-runtime-operator) reference type

	auto decay() const -> decay_type {decay_type ret; copy_n(operator&(), 1, &ret); return ret;}  // NOLINT(fuchsia-default-arguments-calls) complex
	operator decay_type()       const {return decay();}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions,hicpp-explicit-conversion) to allow terse syntax
// #if not defined(__CUDACC__) or not defined(__INTEL_COMPILER)
//  friend auto operator*(decay_type const& lhs, dot_ref const& self) {return lhs*self.decay();}
// #endif
	auto operator+() const -> decay_type {return decay();}

	auto operator==(nrm2_ref const& other) const -> bool {return decay() == other.decay();}
	auto operator!=(nrm2_ref const& other) const -> bool {return decay() != other.decay();}

	// template<class Other>
	// auto operator==(Other const& other) const
	// ->decltype(decay()==other) {
	//  return decay()==other; }
	// template<class Other>
	// auto operator!=(Other const& other) const
	// ->decltype(decay()!=other) {
	//  return decay()!=other; }
};

template<class X>
[[nodiscard]]
auto nrm2(X const& x) {  // NOLINT(readability-identifier-length) BLAS naming
	return nrm2_ref<X>{x};
}

namespace operators {
	using std::norm;
	template<class A1D, class Real = decltype(norm(std::declval<typename A1D::value_type>()))>//decltype(norm(std::declval<typename A1D::value_type>()))>
	[[nodiscard]] auto operator^(A1D const& array, int n)
	->decltype(std::pow(Real{blas::nrm2(array)}, n)) {
		return std::pow(Real{blas::nrm2(array)}, n); }

	template<class A1D>
	[[nodiscard]] auto abs(A1D const& array) {
		return blas::nrm2(array);
	}

	template<class A1D>
	[[nodiscard]] auto norm(A1D const& array) {
		auto const sqrt = +blas::nrm2(array);
		return sqrt*sqrt;
	}

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
//  namespace blas = multi::blas;
//  multi::array<double, 2> const cA = {
//      {1.,  2.,  3.,  4.},
//      {5.,  6.,  7.,  8.},
//      {9., 10., 11., 12.}
//  };

//  double n;
//  BOOST_REQUIRE( blas::nrm2(rotated(cA)[1], n) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( n == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//  double n2 = blas::nrm2(rotated(cA)[1]);
//  BOOST_REQUIRE( n == n2 );

//  multi::array<double, 1> R(4);
//  blas::nrm2( rotated(cA)[1], R[2]);
//  BOOST_REQUIRE( R[2] ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//  multi::array<double, 0> R0;
//  blas::nrm2( rotated(cA)[1], R0);
//  BOOST_REQUIRE( R0 ==  std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//  BOOST_REQUIRE( blas::nrm2(rotated(cA)[1]) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );

//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_blas_nrm2_operators){
//  multi::array<double, 1> X = {1.1,2.1,3.1, 4.1};
//  double n; multi::blas::nrm2(X, n);
//  BOOST_REQUIRE( n == multi::blas::nrm2(X) );

//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case){
//  using complex = std::complex<double>;
//  multi::array<complex, 2> const cA = {
//      {1.,  2.,  3.,  4.},
//      {5.,  6.,  7.,  8.},
//      {9., 10., 11., 12.}
//  };

//  using multi::blas::nrm2;
//  double n; 
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

//#if 0
//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_thrust){
//  using complex = thrust::complex<double>;
//  multi::array<complex, 2> const cA = {
//      {1.,  2.,  3.,  4.},
//      {5.,  6.,  7.,  8.},
//      {9., 10., 11., 12.}
//  };

//  using multi::blas::nrm2;
//  double n;
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//}

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex_real_case_types){
//  boost::mpl::for_each<boost::mpl::list<
//      std   ::complex<double>, 
//      thrust::complex<double>//,
//  //  boost::multi::complex<double> // TODO make this work
//  >>([](auto cplx){
//      multi::array<decltype(cplx), 2> const cA = {
//          {1.,  2.,  3.,  4.},
//          {5.,  6.,  7.,  8.},
//          {9., 10., 11., 12.}
//      };

//      using multi::blas::nrm2;
//      double n;
//      BOOST_REQUIRE( nrm2(rotated(cA)[1], n) == std::sqrt( 2.*2. + 6.*6 + 10.*10.) );
//      BOOST_REQUIRE( nrm2(rotated(cA)[1])    == n );
//  });
//}
//#endif

//BOOST_AUTO_TEST_CASE(multi_adaptor_multi_nrm2_complex){
//  using complex = std::complex<double>; complex const I{0,1};
//  multi::array<complex, 2> const cA = {
//      {1.,  2. + 1.*I,  3.,  4.},
//      {5.,  6. + 4.*I,  7.,  8.},
//      {9., 10. - 3.*I, 11., 12.}
//  };

//  using multi::blas::nrm2;
//  double n;
//  BOOST_REQUIRE( nrm2(rotated(cA)[1], n)   == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );
//  BOOST_REQUIRE( nrm2(rotated(cA)[1])      == std::sqrt( norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) ) );

//  using namespace multi::blas::operators;
//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^-1) == 1/std::sqrt(norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1])) , boost::test_tools::tolerance(1e-15) );
//  BOOST_TEST_REQUIRE( (rotated(cA)[1]^2) == norm(cA[0][1]) + norm(cA[1][1]) + norm(cA[2][1]) , boost::test_tools::tolerance(1e-15) );
//}

//#endif
#endif
