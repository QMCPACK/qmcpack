// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element transformed"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>
#include<numeric>

namespace multi = boost::multi;

using complex = std::complex<double>;
constexpr complex I{0, 1};  // NOLINT(readability-identifier-length) I imaginary unit

BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_function_reference) {
	multi::array<complex, 1> arr = { 1. + 2.*I,  3. +  4.*I};

	constexpr auto conj = static_cast<complex (&)(complex const&)>(std::conj<double>);

	auto const& conjd_arr = arr.element_transformed(conj);
	BOOST_REQUIRE( conjd_arr[0] == conj(arr[0]) );
	BOOST_REQUIRE( conjd_arr[1] == conj(arr[1]) );

//  Ac[0] = 5. + 4.*I;  // this doesn't compile, good!
	BOOST_REQUIRE( conjd_arr[0] == 1. - 2.*I );

	BOOST_REQUIRE( real(std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{0.})) == std::norm(arr[0]) + std::norm(arr[1]) );
	BOOST_REQUIRE( imag(std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{0.})) == 0.                                    );

	BOOST_REQUIRE( std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{0.}) == std::norm(arr[0]) + std::norm(arr[1]) );
}

BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_lambda) {
	multi::array<complex, 1> arr = { 1. + 2.*I,  3. +  4.*I};

	auto const& conjd_arr = arr.element_transformed([](auto const& cee) {return std::conj(cee);});
	BOOST_REQUIRE( conjd_arr[0] == std::conj(arr[0]) );
	BOOST_REQUIRE( conjd_arr[1] == std::conj(arr[1]) );

//	Ac[0] = 5. + 4.*I;  // this doesn't compile, good!
	BOOST_REQUIRE( conjd_arr[0] == 1. - 2.*I );
}

BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_lambda_with_const_return) {
	multi::array<complex, 1> arr = { 1. + 2.*I,  3. +  4.*I};

	auto const& conjd_arr = arr.element_transformed([](auto const& cee) -> auto const {return std::conj(cee);});  // NOLINT(readability-const-return-type) to disable assignment
	BOOST_REQUIRE( conjd_arr[0] == std::conj(arr[0]) );
	BOOST_REQUIRE( conjd_arr[1] == std::conj(arr[1]) );

//	Ac[0] = 5. + 4.*I;  // this doesn't compile, good!
	BOOST_REQUIRE( conjd_arr[0] == 1. - 2.*I );
}

template<typename ComplexRef> struct Conjd;

constexpr struct Conj_t {  // NOLINT(readability-identifier-naming) for testing
	template<class ComplexRef> constexpr auto operator()(ComplexRef&& zee) const {return Conjd<decltype(zee)>{zee};}
	template<class T> constexpr auto operator()(Conjd<T> const&) const = delete;
	template<class T> constexpr auto operator()(Conjd<T> &&) const = delete;
	template<class T> constexpr auto operator()(Conjd<T> &) const = delete;
} Conj;

template<typename ComplexRef>
struct Conjd {  // NOLINT(readability-identifier-naming) for testing
	using decay_type = decltype( + std::declval<ComplexRef>() );

	constexpr operator decay_type() const {return std::conj(c_);}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions)

	friend constexpr auto operator==(decay_type const& other, Conjd const& self) -> bool {return std::conj(self.c_) == other;}
	friend constexpr auto operator!=(decay_type const& other, Conjd const& self) -> bool {return std::conj(self.c_) != other;}

	friend constexpr auto operator==(Conjd const& self, decay_type const& other) -> bool {return other == std::conj(self.c_);}
	friend constexpr auto operator!=(Conjd const& self, decay_type const& other) -> bool {return other != std::conj(self.c_);}

	friend constexpr auto operator==(Conjd const& self, Conjd const& other) -> bool {return other.c_ == self.c_;}
	friend constexpr auto operator!=(Conjd const& self, Conjd const& other) -> bool {return other.c_ != self.c_;}

	constexpr auto operator=(decay_type const& other) && -> Conjd& {c_ = std::conj(other); return *this;}

 private:
	constexpr explicit Conjd(ComplexRef cee) : c_{cee} {}
	ComplexRef c_;
	friend decltype(Conj);
};

BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_proxy) {
	multi::array<complex, 1> const arr = { 1. + 2.*I,  3. +  4.*I};

	auto const& conj_arr = arr.element_transformed(Conj);
	BOOST_REQUIRE( std::conj(arr[0]) == conj_arr[0] );
	BOOST_REQUIRE( std::conj(arr[1]) == conj_arr[1] );

//  Ac[0] = 5. + 4.*I;  // not allowed, compile error, Ac is const
	BOOST_REQUIRE( conj_arr[0] == 1. - 2.*I );
}

BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_mutable_proxy) {
	multi::array<complex, 1> arr = { 1. + 2.*I,  3. +  4.*I};

	auto&& conj_arr = arr.element_transformed(Conj);  // NOLINT(readability-const-return-type) to disable assignment

	BOOST_REQUIRE( std::conj(arr[0]) == conj_arr[0] );
	BOOST_REQUIRE( std::conj(arr[1]) == conj_arr[1] );

	conj_arr[0] = 5. + 4.*I;
	BOOST_REQUIRE( conj_arr[0] == 5. + 4.*I );
	BOOST_REQUIRE(  arr[0] == 5. - 4.*I );
}

BOOST_AUTO_TEST_CASE(transform_ptr_single_value) {
	complex cee = 1. + 2.*I;

	constexpr auto conj_ro = [](auto const& zee) {return std::conj(zee);};  // NOLINT(readability-const-return-type,clang-diagnostic-ignored-qualifiers) to prevent assignment

	multi::transform_ptr<complex, decltype(conj_ro), complex*> conjd_ceeP{&cee, conj_ro};
	BOOST_REQUIRE( *conjd_ceeP == std::conj(1. + 2.*I) );
}

BOOST_AUTO_TEST_CASE(transform_ptr_1D_array) {
	multi::array<complex, 1> arr = { 1. + 2.*I,  3. +  4.*I};

	constexpr auto conj_ro = [](auto const& zee) {return std::conj(zee);};  // NOLINT(readability-const-return-type,clang-diagnostic-ignored-qualifiers) to prevent assignment

	auto const& conjd_arr = arr.element_transformed(conj_ro);
	BOOST_REQUIRE( conjd_arr[0] == conj_ro(arr[0]) );
	BOOST_REQUIRE( conjd_arr[1] == conj_ro(arr[1]) );

//  Ac[0] = 5. + 4.i;  // doesn't compile thanks to the `auto const` in the `conj` def
}

BOOST_AUTO_TEST_CASE(arthur_odwyer_array_transform_int) {
	struct S {  // NOLINT(readability-identifier-naming)
    	int a;
    	int b;
	};

	multi::array<S, 1> arr({2}, S{});
	auto&& ref = arr.element_transformed(&S::a);
	ref[0] = 99.;

	BOOST_REQUIRE( arr[0].a == 99. );

	auto const& cref = arr.element_transformed(&S::a);
	BOOST_REQUIRE( cref[0] == 99. );
//  cr[0] = 99.;  // compile error "assignment of read-only location"
}

BOOST_AUTO_TEST_CASE(arthur_odwyer_array_transform_int_array) {
	struct S {  // NOLINT(readability-identifier-naming)
    	int a[10];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing
    	int b;
	};

	multi::array<S, 1> vec({2}, S{});
	auto&& ref = vec.element_transformed(&S::a);
	ref[0][1] = 99.;

	BOOST_REQUIRE( ref[0][1] == 99. );
	BOOST_REQUIRE( vec[0].a[1] == 99. );

	auto const& cref = vec.element_transformed(&S::a);
	BOOST_REQUIRE( cref[0][1] == 99. );
//  cref[0][1] = 99.;  // compile error "assignment of read-only location"
}

BOOST_AUTO_TEST_CASE(indirect_transformed) {
	std::vector<double> vec = {0., 1.1, 2.2, 3.3, 4.4, 5.5};

	using index_t = std::vector<double>::size_type;

	multi::array<index_t, 1> const arr = {4, 3, 2, 1, 0};

	auto&& indirect_v = arr.element_transformed([&vec](index_t idx) noexcept -> double& {return vec[idx];});

	BOOST_REQUIRE(  indirect_v[1] ==  vec[3] );
	BOOST_REQUIRE( &indirect_v[1] == &vec[3] );

	indirect_v[1] = 99.;
	BOOST_REQUIRE(  vec[3] ==  99. );

	for(auto&& elem : indirect_v) {elem = 88.;}
	BOOST_REQUIRE(  vec[3] ==  88. );

	auto const& const_indirect_v = indirect_v;  (void)const_indirect_v;
//  const_indirect_v[1] = 999.;  // does not compile, good!
	BOOST_REQUIRE(  const_indirect_v[3] ==  88. );
}

BOOST_AUTO_TEST_CASE(indirect_transformed_carray) {
	double carr[5][3] = { // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing legacy types
		{ 0.0,  1.0,  2.0},
		{10.0, 11.0, 12.0},
		{20.0, 21.0, 22.0},
		{30.0, 31.0, 32.0},
		{40.0, 41.0, 42.0}
	};

	using index_t = std::vector<double>::size_type;
	multi::array<index_t, 1> const arr = {4, 3, 2, 1, 0};

	auto&& indirect_v = arr.element_transformed([&carr](index_t idx) noexcept -> double(&)[3] {return carr[idx];});  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)

	BOOST_REQUIRE( &indirect_v[1][2] ==  &carr[3][2] );
	BOOST_REQUIRE(  indirect_v[1][2] ==  32.0 );

	indirect_v[1][2] = 11111.0;
	BOOST_TEST   (  indirect_v[1][2] ==  11111.0 );

	auto const& const_indirect_v = indirect_v;

	BOOST_TEST(  const_indirect_v[1][2] ==  11111.0 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) testing legacy type
//  const_indirect_v[1][2] = 999.;  // doesn't compile, good!
}
