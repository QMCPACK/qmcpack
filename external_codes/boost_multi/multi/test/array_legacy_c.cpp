// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, rotated, subarray, dimens...

#include <boost/core/lightweight_test.hpp>

#include <array>    // for array
#include <complex>  // for complex
#include <tuple>    // for apply  // IWYU pragma: keep
// IWYU pragma: no_include <type_traits>  // for remove_reference<>::type
#include <utility>  // for move

namespace multi = boost::multi;

namespace fake {

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing a legacy interface
using fftw_complex = double[2];

namespace {

void fftw_plan_dft(int rank, int const* n, fftw_complex* in, fftw_complex* out, int sign, unsigned flags);

void fftw_plan_dft(int rank, int const* n, fftw_complex* in, fftw_complex* out, int sign, unsigned flags) {
	(void)rank, (void)n, (void)in, (void)out, (void)sign, (void)flags;
}

}  // end unnamed namespace

}  // end namespace fake

namespace {
constexpr auto f2(multi::array_ref<double, 1>&& array) -> double& { return std::move(array)[2]; }
}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(array_legacy_c)
	{
		using complex = std::complex<double>;

		multi::array<complex, 2> const in = {
			{{150.0, 0.0}, {16.0, 0.0}, {17.0, 0.0}, {18.0, 0.0}, {19.0, 0.0}},
			{  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0}},
			{{100.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}, {13.0, 0.0}, {14.0, 0.0}},
			{ {50.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0},  {9.0, 0.0}},
		};

		multi::array<std::complex<double>, 2> out(in.extents());

		BOOST_TEST( dimensionality(out) == dimensionality(in) );
		BOOST_TEST( sizes(out) == sizes(in) );

		static_assert(sizeof(complex) == sizeof(fake::fftw_complex), "!");
		fake::fftw_plan_dft(
			decltype(in)::dimensionality,
			std::apply([](auto... sizes) { return std::array<int, 2>{{static_cast<int>(sizes)...}}; }, in.sizes()).data(),
			// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) testing legacy code
			reinterpret_cast<fake::fftw_complex*>(const_cast<complex*>(in.data_elements())),  // NOSONAR
			// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): testing legacy code
			reinterpret_cast<fake::fftw_complex*>(out.data_elements()),
			1, 0
		);
	}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif
	// BOOST_AUTO_TEST_CASE(array_legacy_c_2)
	{
		double arr[5] = {150.0, 16.0, 17.0, 18.0, 19.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		BOOST_TEST( &f2(arr) == &arr[2] );
	}
#ifdef __clang__
#pragma clang diagnostic pop
#endif

	return boost::report_errors();
}
