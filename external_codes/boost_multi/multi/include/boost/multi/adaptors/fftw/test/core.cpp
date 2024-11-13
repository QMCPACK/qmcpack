// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/array.hpp>

#include <algorithm>    // for for_each, tra...
#include <complex>      // for operator*
#include <iterator>     // for begin, end
#include <numeric>      // for accumulate, iota
#include <random>       // std::mt19937_64
#include <string>       // for operator""s
#include <type_traits>  // for decay_t, enab...
#include <utility>      // for forward

namespace {

namespace multi = boost::multi;
namespace fftw  = multi::fftw;

template<class M> auto power(M const& elem) -> decltype(std::norm(elem)) { return std::norm(elem); }

template<class M, class = std::enable_if_t<(M::rank::value >= 1)>>  // DELETE((M::rank::value < 1))>
auto power(M const& array) {
	return accumulate(begin(array), end(array), 0.0, [](auto const& alpha, auto const& omega) { return alpha + power(omega); });
}

struct sum_power {
	template<class A, class B> auto operator()(A const& alpha, B const& omega) const { return alpha + power(omega); }
};

}  // end anonymous namespace

template<class T> class randomizer {
	std::mt19937_64 gen_;  // NOSONAR rng good enough for the test

 public:
	explicit randomizer(unsigned int seed) : gen_(seed) {}

	template<class M, class R = typename std::decay_t<M>::reference> void operator()(M&& arr) {
		std::for_each(std::begin(std::forward<M>(arr)), std::end(std::forward<M>(arr)), [self = this](R elem) { self->operator()(elem); });
	}
	void operator()(T& elem) {  // NOLINT(runtime/references) passing by reference
		std::normal_distribution<T> gauss;
		elem = gauss(gen_);
	}
};

template<class T> class randomizer<std::complex<T>> {
	std::mt19937_64 gen_;  // NOSONAR rng good enough for the test

 public:
	explicit randomizer(unsigned int seed) : gen_(seed) {}

	template<class M, class R = typename std::decay_t<M>::reference> void operator()(M&& arr) {
		std::for_each(std::begin(std::forward<M>(arr)), std::end(std::forward<M>(arr)), [self = this](R elem) { self->operator()(elem); });
	}
	void operator()(std::complex<T>& zee) {  // NOLINT(runtime/references) : passing by reference
		std::normal_distribution<T> gauss;
		zee = std::complex<T>(gauss(gen_), gauss(gen_));
	}
};

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE)  /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	fftw::environment const env;

	BOOST_AUTO_TEST_CASE(fftw_2D_identity_2) {  //, *boost::unit_test::tolerance(0.0001)) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const in = {
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		multi::array<complex, 2> out(extensions(in));

		multi::fftw::dft_forward({false, false}, in, out);  // out = in;

		BOOST_TEST( in[2][3].real() == out[2][3].real() );
		BOOST_TEST( in[2][3].imag() == out[2][3].imag() );

		BOOST_TEST( out == in );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_many) {  // , *boost::unit_test::tolerance(0.0001)) {
		using complex = std::complex<double>;

		auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const in = {
			{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
		multi::array<complex, 2> out(extensions(in));

		using multi::fftw::dft_forward;

		multi::fftw::dft_forward({false, false}, in.rotated(), out.rotated());
		BOOST_TEST( in == out );
	}

	BOOST_AUTO_TEST_CASE(fftw_many1_from_2) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		std::random_device       dev;
		multi::array<complex, 2> in({3, 10});
		randomizer<complex>{dev()}(in);
		multi::array<complex, 2> out({3, 10});
		fftw::dft_forward({false, true}, in, out);

		multi::array<complex, 2> out2({3, 10});
		std::transform(in.begin(), in.end(), out2.begin(), out2.begin(), [](auto const& in_elem, auto&& out2_elem) {
			fftw::dft_forward({true}, in_elem, out2_elem);
			return std::forward<decltype(out2_elem)>(out2_elem);
		});

		BOOST_TEST(out2 == out);
	}

	BOOST_AUTO_TEST_CASE(fftw_many2_from_3) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		std::random_device       dev;
		multi::array<complex, 3> in({3, 5, 6});
		randomizer<complex>{dev()}(in);
		multi::array<complex, 3> out({3, 5, 6});
		fftw::dft_forward({false, true, true}, in, out);

		multi::array<complex, 3> out2({3, 5, 6});
		std::transform(in.begin(), in.end(), out2.begin(), out2.begin(), [](auto const& in_elem, auto&& out2_elem) {
			fftw::dft_forward({true, true}, in_elem, out2_elem);
			return std::forward<decltype(out2_elem)>(out2_elem);
		});

		BOOST_TEST(out2 == out);
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_power_plan) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in({16, 16});
		std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		multi::array<complex, 2> out(extensions(in));
		auto const               pln = multi::fftw::plan::forward({true, true}, in.base(), in.layout(), out.base(), out.layout());
		pln.execute(in.base(), out.base());
		BOOST_TEST( power(in) - power(out)/num_elements(out) < 1e-7 );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in({16, 16});
		std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		multi::array<complex, 2> out(extensions(in));
		auto const               pln = multi::fftw::plan::forward({true, true}, in.base(), in.layout(), out.base(), out.layout());
		pln.execute(in.base(), out.base());
		BOOST_TEST( power(in) - power(out)/num_elements(out) < 1e-8 );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern_measure) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in({16, 16});
		std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		multi::array<complex, 2> out(extensions(in));
		auto const               pln = multi::fftw::plan::forward({true, true}, in.base(), in.layout(), out.base(), out.layout());
		pln.execute(in.base(), out.base());
		BOOST_TEST( power(in) - power(out)/num_elements(out) < 1e-8 );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_power_dft) {
		using complex                 = std::complex<double>;
		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in({16, 16});
		std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		multi::array<complex, 2> out(extensions(in));
		multi::fftw::dft_forward({true, true}, in, out);
		BOOST_TEST( power(in) - power(out)/num_elements(out) < 1e-8 );
	}

	BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place_over_ref_inplace) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 3> io({4, 4, 4});
		std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		auto const powerin = power(io);

		//  fftw::dft_inplace(multi::array_ref<complex, 3>(io.data(), io.extensions()), fftw::forward);

		fftw::dft_forward(
			{true, true, true},
			multi::array_ref<complex, 3>(data_elements(io), extensions(io)),
			multi::array_ref<complex, 3>(data_elements(io), extensions(io))
		);
		BOOST_TEST( powerin - power(io)/num_elements(io) < 1e-10 );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> const in = {
			{100.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
			{  3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
			{  4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
			{  3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
			{ 31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
		};
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_naive_square) {
		using complex = std::complex<double>;

		[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

		multi::array<complex, 2> in = {
			{100.0 + 2.0 * I, 9.0 - 1.0 * I, 2.0 + 4.0 * I},
			{  3.0 + 3.0 * I, 7.0 - 4.0 * I, 1.0 + 9.0 * I},
			{  4.0 + 1.0 * I, 5.0 + 3.0 * I, 2.0 + 4.0 * I},
		};
		multi::array<complex, 2> const in_transpose = in.transposed();
		in                                          = in.transposed();
		BOOST_TEST( in != in_transpose );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_nonpod) {
		using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s
		multi::array<std::string, 2> in = {
			{"100.0 + 2.0*I"s,  "9.0 - 1.0*I"s, "2.0 +  4.0*I"s},  // NOLINT(misc-include-cleaner) bug in clang-tidy 18
			{  "3.0 + 3.0*I"s,  "7.0 - 4.0*I"s, "1.0 +  9.0*I"s},
			{  "4.0 + 1.0*I"s,  "5.0 + 3.0*I"s, "2.0 +  4.0*I"s},
			{  "3.0 - 1.0*I"s,  "8.0 + 7.0*I"s, "2.0 +  1.0*I"s},
			{ "31.0 - 1.0*I"s, "18.0 + 7.0*I"s, "2.0 + 10.0*I"s},
		};
		multi::array<std::string, 2> const in_transpose = in.transposed();
		in                                              = in.transposed();
		BOOST_TEST( in != in_transpose );
	}

	BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_nonpod_square) {
		using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

		multi::array<std::string, 2> in = {
			{"100.0 + 2.0*I"s, "9.0 - 1.0*I"s, "2.0 +  4.0*I"s},
			{  "3.0 + 3.0*I"s, "7.0 - 4.0*I"s, "1.0 +  9.0*I"s},
			{  "4.0 + 1.0*I"s, "5.0 + 3.0*I"s, "2.0 +  4.0*I"s},
		};
		multi::array<std::string, 2> const in_transpose = in.transposed();
		in                                              = in.transposed();
		BOOST_TEST( in != in_transpose );
	}

	return boost::report_errors();
}
