// Copyright 2022-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifdef _MSC_VER
#pragma warning(disable : 4244)  // warning C4244: 'initializing': conversion from '_Ty' to '_Ty', possible loss of data
#endif

#include <boost/multi/array.hpp>  // for transform_ptr, array, subarray

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

// IWYU pragma: no_include <algorithm>                        // for copy  // for GNU stdlib
// IWYU pragma: no_include <type_traits>                      // for declval  // for GNU stdlib
#include <complex>   // IWYU pragma: keep  // for complex, operator*, operator+
#include <iterator>  // IWYU pragma: keep  // for weakly_incrementable
#include <utility>   // IWYU pragma: keep  // for declval, forward
#include <vector>    // IWYU pragma: keep  // for vector

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<concepts>)
#include <concepts>  // IWYU pragma: keep
#include <ranges>    // IWYU pragma: keep
#endif

template<typename ComplexRef> struct Conjd;  // NOLINT(readability-identifier-naming) for testing

struct Conj_t {  // NOLINT(readability-identifier-naming,misc-use-internal-linkage) for testing
	template<class ComplexRef> constexpr auto operator()(ComplexRef&& zee) const noexcept { return Conjd<decltype(zee)>{std::forward<ComplexRef>(zee)}; }
	template<class T> constexpr auto          operator()(Conjd<T> const&) const = delete;
	template<class T> constexpr auto          operator()(Conjd<T>&&) const      = delete;
	template<class T> constexpr auto          operator()(Conjd<T>&) const       = delete;
};
inline constexpr Conj_t Conj;

template<typename ComplexRef>
struct Conjd {  // NOLINT(readability-identifier-naming) for testing
	using decay_type = decltype(+std::declval<ComplexRef>());

	// NOLINTNEXTLINE(google-explicit-constructor,hicpp-explicit-conversions)
	// explicit constexpr operator decay_type() const {return std::conj(c_); }  // NOSONAR(cpp:S1709)

	friend constexpr auto operator==(decay_type const& other, Conjd const& self) -> bool { return std::conj(self.c_) == other; }
	friend constexpr auto operator!=(decay_type const& other, Conjd const& self) -> bool { return std::conj(self.c_) != other; }

	friend constexpr auto operator==(Conjd const& self, decay_type const& other) -> bool { return other == std::conj(self.c_); }
	friend constexpr auto operator!=(Conjd const& self, decay_type const& other) -> bool { return other != std::conj(self.c_); }

	friend constexpr auto operator==(Conjd const& self, Conjd const& other) -> bool { return other.c_ == self.c_; }
	friend constexpr auto operator!=(Conjd const& self, Conjd const& other) -> bool { return other.c_ != self.c_; }

	constexpr auto operator=(decay_type const& other) && -> Conjd& {
		c_ = std::conj(other);
		return *this;
	}

 private:
	constexpr explicit Conjd(ComplexRef& cee) : c_{cee} {}
	ComplexRef& c_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) can be a reference
	friend decltype(Conj);
};

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_function_reference)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		// conste xpr
		auto const conj = static_cast<complex (&)(complex const&)>(std::conj<double>);

		auto const& conjd_arr = arr.element_transformed(conj);
		BOOST_TEST( conjd_arr[0] == conj(arr[0]) );
		BOOST_TEST( conjd_arr[1] == conj(arr[1]) );

		//  Ac[0] = 5. + 4.*I;  // this doesn't compile, good!
		BOOST_TEST( conjd_arr[0] == 1.0 - 2.0*I );

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<concepts>)
		auto conjd_arr_beg = conjd_arr.begin();
		conjd_arr_beg      = conjd_arr.end();

		static_assert(std::movable<decltype(conjd_arr_beg)>);
		static_assert(std::weakly_incrementable<decltype(conjd_arr_beg)>);  // NOLINT(misc-include-cleaner)
		BOOST_TEST( conjd_arr.begin() == std::ranges::begin(conjd_arr) );
#endif

		// BOOST_REQUIRE_CLOSE(real(std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{ 0.0, 0.0 })), std::norm(arr[0]) + std::norm(arr[1]), 1E-6);
		// BOOST_REQUIRE_CLOSE(imag(std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{ 0.0, 0.0 })), 0.0, 1E-6);

		// BOOST_TEST_REQUIRE( std::inner_product(arr.begin(), arr.end(), conjd_arr.begin(), complex{0.0, 0.0}) == std::norm(arr[0]) + std::norm(arr[1]) );
	}

	// BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_lambda)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		// g++ -std=20 needs the transformation (lambda) to be noexcept
		auto const& conjd_arr = arr.element_transformed([](auto const& cee) noexcept { return std::conj(cee); });
		BOOST_TEST( conjd_arr[0] == std::conj(arr[0]) );
		BOOST_TEST( conjd_arr[1] == std::conj(arr[1]) );

		//  Ac[0] = 5. + 4.*I;  // this doesn't compile, good!
		BOOST_TEST( conjd_arr[0] == 1.0 - 2.0*I );
	}

	// BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_lambda_with_const_return)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		// g++ -std=20 needs the transformation (lambda) to be noexcept
		auto&& conjd_arr = arr.element_transformed([](auto const& cee) noexcept { return std::conj(cee); });

		BOOST_TEST( conjd_arr[0] == std::conj(arr[0]) );
		BOOST_TEST( conjd_arr[1] == std::conj(arr[1]) );

		// conjd_arr[0] = 5.0 + 4.0*I;  // this compiles, but classes the implement operator= naively can be misleading here
		BOOST_TEST( conjd_arr[0] == 1.0 - 2.0*I );
	}

	// BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_proxy)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> const arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		auto const& conj_arr = arr.element_transformed(Conj);
		BOOST_TEST( std::conj(arr[0]) == conj_arr[0] );
		BOOST_TEST( std::conj(arr[1]) == conj_arr[1] );

		//  Ac[0] = 5. + 4.*I;  // not allowed, compile error, Ac is const
		BOOST_TEST( conj_arr[0] == 1.0 - 2.0*I );
	}

	// BOOST_AUTO_TEST_CASE(element_transformed_1D_conj_using_mutable_proxy)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		auto&& conj_arr = arr.element_transformed(Conj);  // NOLINT(readability-const-return-type) to disable assignment

		BOOST_TEST( std::conj(arr[0]) == conj_arr[0] );
		BOOST_TEST( std::conj(arr[1]) == conj_arr[1] );

		conj_arr[0] = 5.0 + 4.0 * I;
		BOOST_TEST( conj_arr[0] == 5.0 + 4.0*I );
		BOOST_TEST(  arr[0] == 5.0 - 4.0*I );
	}

	// BOOST_AUTO_TEST_CASE(transform_ptr_single_value)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		complex cee = 1.0 + 2.0 * I;

		// NOLINTNEXTLINE(readability-const-return-type,clang-diagnostic-ignored-qualifiers) to prevent assignment
		constexpr auto conj_ro = [](auto const& zee) noexcept {
			return std::conj(zee);
		};  // g++ -std=20 needs the transformation (lambda) to be noexcept

		multi::transform_ptr<complex, decltype(conj_ro), complex*> const conjd_ceeP{&cee, conj_ro};
		BOOST_TEST( *conjd_ceeP == std::conj(1.0 + 2.0*I) );
	}

	// BOOST_AUTO_TEST_CASE(transform_ptr_1D_array)
	{
		using complex = std::complex<double>;
		auto const I  = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) I imaginary unit

		multi::array<complex, 1> arr = {1.0 + 2.0 * I, 3.0 + 4.0 * I};

		// NOLINT(readability-const-return-type,clang-diagnostic-ignored-qualifiers) to prevent assignment
		constexpr auto conj_ro = [](auto const& zee) noexcept {
			return std::conj(zee);
		};  // g++ -std=20 needs the transformation (lambda) to be noexcept

		auto const& conjd_arr = arr.element_transformed(conj_ro);
		BOOST_TEST( conjd_arr[0] == conj_ro(arr[0]) );
		BOOST_TEST( conjd_arr[1] == conj_ro(arr[1]) );

		//  Ac[0] = 5. + 4.i;  // doesn't compile thanks to the `auto const` in the `conj` def
	}

	// BOOST_AUTO_TEST_CASE(arthur_odwyer_array_transform_int)
	{
		struct S {  // NOLINT(readability-identifier-naming)
			int a;
			int b;
		};

		multi::array<S, 1> arr({2}, S{});

		auto&& ref = arr.element_transformed(&S::a);

		ref[0] = 990;

		BOOST_TEST( arr[0].a == 990 );

		auto const& cref = arr.element_transformed(&S::a);
		BOOST_TEST( cref[0] == 990 );
		//  cr[0] = 99.;  // compile error "assignment of read-only location"
	}

	// BOOST_AUTO_TEST_CASE(arthur_odwyer_array_transform_int_array)
	{
		struct S {      // NOLINT(readability-identifier-naming)
			int a[10];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing
			int b;
		};

		multi::array<S, 1> vec({2}, S{});

		auto&& ref = vec.element_transformed(&S::a);

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		ref[0][1] = 990;

		BOOST_TEST( ref[0][1] == 990 );

		BOOST_TEST( vec[0].a[1] == 990 );

		auto const& cref = vec.element_transformed(&S::a);
		BOOST_TEST( cref[0][1] == 990 );
		//  cref[0][1] = 990;  // compile error "assignment of read-only location"

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(indirect_transformed)
	{
		std::vector<int> vec = {00, 11, 22, 33, 44, 55};  // std::vector NOLINT(fuchsia-default-arguments-calls)

		using index_t = std::vector<int>::size_type;

		multi::array<index_t, 1> const arr = {4, 3, 2, 1, 0};

		auto&& indirect_v = arr.element_transformed([&vec](index_t idx) noexcept -> int& { return vec[idx]; });

		BOOST_TEST(  indirect_v[1] ==  vec[3] );
		BOOST_TEST( &indirect_v[1] == &vec[3] );

		indirect_v[1] = 990;
		BOOST_TEST(  vec[3] ==  990 );

		//  for(auto&& elem : indirect_v) {elem = 88.;}
		//  std::fill(indirect_v.begin(), indirect_v.end(), 88.0);

#ifndef _MSC_VER
		indirect_v.fill(880);
		BOOST_TEST(  vec[3] ==  880 );

		auto const& const_indirect_v = indirect_v;
		(void)const_indirect_v;
		// const_indirect_v[1] = 9990;  // does not compile, good!
		BOOST_TEST(const_indirect_v[3] ==  880);
#endif
	}

	// BOOST_AUTO_TEST_CASE(indirect_transformed_carray)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing legacy types
		int carr[5][3] = {
			{ 00,  10,  20},
			{100, 110, 120},
			{200, 210, 220},
			{300, 310, 320},
			{400, 410, 420},
		};

		using index_t = std::vector<int>::size_type;

		multi::array<index_t, 1> const arr = {4, 3, 2, 1, 0};

		// cppcheck-suppress CastAddressToIntegerAtReturn ;  // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		auto&& indirect_v = arr.element_transformed([&carr](index_t idx) noexcept -> int (&)[3] { return carr[idx]; });
		// ^^^  TODO(correaa) investigate, portability: Returning an address value in a function with integer return type is not portable.

		BOOST_TEST( &indirect_v[1][2] ==  &carr[3][2] );
		BOOST_TEST(  indirect_v[1][2] ==  320 );

		indirect_v[1][2] = 111110;
		BOOST_TEST   (  indirect_v[1][2] ==  111110 );

		auto const& const_indirect_v = indirect_v;

		BOOST_TEST(  const_indirect_v[1][2] ==  111110 );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic) testing legacy type

		// const_indirect_v[1][2] = 999.0;  // doesn't compile, good!

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	{
		multi::array<int, 1> const arr = {4, 3, 2, 1, 0};

		auto const& arr_plus_1 = arr.element_transformed([](auto val) noexcept { return val + 1; });

		multi::array<int, 1> const arr2 = arr_plus_1;

		BOOST_TEST(( arr2 == multi::array<int, 1>{5, 4, 3, 2, 1} ));
	}
	{
		multi::array<int, 1> const arr = {4, 3, 2, 1, 0};

		auto const& arr_plus_1 = arr.element_transformed(multi::val([](auto val) noexcept { return val + 1; }));

		multi::array<int, 1> const arr2 = arr_plus_1;

		BOOST_TEST(( arr2 == multi::array<int, 1>{5, 4, 3, 2, 1} ));
	}

	return boost::report_errors();
}
