// Copyright 2020-2024 Alfredo A. Correa

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/array.hpp>

#include <complex>
#include <numeric>  // for std::transform_reduce

namespace multi = boost::multi;

template<class M> auto power(M const& array) {
	return std::accumulate(array.elements().begin(), array.elements().end(), 0.0, [](auto e1, auto e2) { return std::move(e1) + std::norm(e2); });
	//  return std::transform_reduce(array.elements().begin(), array.elements().end(), 0.0, std::plus<>{}, [](auto zee) { return std::norm(zee); });
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	multi::fftw::environment const env;

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_move) {
		using complex                 = std::complex<double>;
		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in = {
			{100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{ 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		BOOST_TEST( in[1][1] == 7.0 - 4.0*I );

		auto const  in_copy = in;
		auto* const in_base = in.base();
		BOOST_TEST( in_base == in.base() );

		// in = multi::fftw::ref(in);

		// BOOST_TEST( in == in_copy );
		// BOOST_TEST( in_base == in.base() );  // prove no allocation
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed) {
		using complex                 = std::complex<double>;
		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in = {
			{100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{ 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		BOOST_TEST( in[1][1] == 7.0 - 4.0*I );

		auto const  in_copy = in;
		auto* const in_base = in.base();
		BOOST_TEST( in_base == in.base() );
		BOOST_TEST( in.size() == 5 );

		//  in = multi::fftw::ref(in).transposed();

		//  BOOST_TEST( in.size() == 3 );
		//  BOOST_TEST( in == in_copy.transposed() );  // prove correctness
		//  BOOST_TEST( in_base == in.base() );        // prove no allocation
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_naive) {
		using complex                 = std::complex<double>;
		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in = {
			{100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{ 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		BOOST_TEST( in[1][1] == 7.0 - 4.0*I );

		auto const  in_copy = in;
		auto* const in_base = in.base();
		BOOST_TEST( in_base == in.base() );
		BOOST_TEST( in.size() == 5 );

		in = in.transposed();  // this is UB

		BOOST_TEST( in.size() == 3 );
		//  BOOST_TEST( in != in_copy.transposed() );  // prove it is incorrect
		BOOST_TEST( in_base == in.base() );  // prove no allocation
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_naive_copy) {
		using complex                 = std::complex<double>;
		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in = {
			{100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{ 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		BOOST_TEST( in[1][1] == 7.0 - 4.0*I );

		auto const  in_copy = in;
		auto* const in_base = in.base();
		BOOST_TEST( in_base == in.base() );
		BOOST_TEST( in.size() == 5 );

		in = +in.transposed();

		BOOST_TEST( in.size() == 3 );
		BOOST_TEST( in == in_copy.transposed() );  // prove correctness
		BOOST_TEST( in_base != in.base() );        // prove no allocation
	}

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_fft_copy) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  multi::array<complex, 2> in2 = multi::fftw::fft(in);

	//  BOOST_TEST( power(in2)/num_elements(in2) - power(in_copy) < 1e-8 );
	//  BOOST_TEST( in2.base() != in_base );
	//  BOOST_TEST( not in.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// }

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_copyconstruct) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  multi::array<complex, 2> in2 = multi::fftw::ref(in).transposed();

	//  BOOST_TEST( in2 == in_copy.transposed() );
	//  BOOST_TEST( in2.base() != in_base );
	//  BOOST_TEST( in .base() == in_base );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// }

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_moveconstruct) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  multi::array<complex, 2> in2 = multi::fftw::ref(std::move(in)).transposed();

	//  BOOST_TEST( in2 == in_copy.transposed() );
	//  BOOST_TEST( in2.base() == in_base );
	//  BOOST_TEST( in.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// }

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_moveconstruct_implicit) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	// #if not defined(__INTEL_COMPILER)  // TODO(correaa) problem with icpc 2022.3.0.8751
	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  auto in2 = +multi::fftw::ref(std::move(in)).transposed();

	//  BOOST_TEST( in2 == in_copy.transposed() );
	// #if not defined(__NVCOMPILER)  // these tests fail with nvc++ 22.9, 23.1
	//  BOOST_TEST( in2.base() == in_base );
	//  BOOST_TEST( in.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// #endif
	// #endif
	// }

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_moveassign_from_temp) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	// #if not defined(__INTEL_COMPILER)  // TODO(correaa) problem with icpc 2022.3.0.8751
	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  multi::array<complex, 2> in2;
	//  in2 = static_cast<multi::array<complex, 2>>(multi::fftw::ref(std::move(in)).transposed());

	//  BOOST_TEST( in2 == in_copy.transposed() );
	// #if not defined(__NVCOMPILER)  // these tests fail with nvc++ 22.9, 23.1
	//  BOOST_TEST( in2.base() == in_base );
	//  BOOST_TEST( in.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// #endif
	// #endif
	// }

	// BOOST_AUTO_TEST_CASE(fftw_2D_const_range_transposed_moveassign) {
	//  using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	// #if not defined(__INTEL_COMPILER)  // TODO(correaa) problem with icpc 2022.3.0.8751
	//  multi::array<complex, 2> in = {
	//      {100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
	//      {  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
	//      {  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
	//      { 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	//  };

	//  auto const in_copy  = in;
	//  auto* const in_base = in.base();

	//  multi::array<complex, 2> in2;
	//  in2 = multi::fftw::ref(std::move(in)).transposed();

	//  BOOST_TEST( in2 == in_copy.transposed() );
	// #if not defined(__NVCOMPILER)  // these tests fail with nvc++ 22.9, 23.1
	//  BOOST_TEST( in2.base() == in_base );
	//  BOOST_TEST( in.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
	// #endif
	// #endif
	// }

	return boost::report_errors();
}
