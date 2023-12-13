// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW memory"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/fftw.hpp>
#include <multi/array.hpp>

#include<chrono>
#include<iostream>
#include<random>

#include<fftw3.h>

class watch : private std::chrono::high_resolution_clock{
	std::string label;
	time_point start = now();

 public:
	explicit watch(std::string label) : label{std::move(label)} {}
	~watch(){std::cerr<< label<<": "<< std::chrono::duration<double>(now() - start).count() <<" sec"<<std::endl;}
};

template<class T> struct randomizer {
	template<class M> void operator()(M&& arr) const {
		std::for_each(arr.begin(), arr.end(), [&self=*this](auto&& elem) {self.operator()(elem);});
	}
	void operator()(T& elem) const {  // NOLINT(runtime/references) passing by reference
		static std::random_device dev; static std::mt19937 gen{dev()}; static std::normal_distribution<T> gauss;
		elem = gauss(gen);
	}
};

template<class T> struct randomizer<std::complex<T>> {
	template<class M> void operator()(M&& arr) const {
		std::for_each(arr.begin(), arr.end(), [&self=*this](auto&& elem) {self.operator()(elem);});
	}
	void operator()(std::complex<T>& zee) const {  // NOLINT(runtime/references) : passing by reference
		static std::random_device dev; static std::mt19937 gen{dev()}; static std::normal_distribution<T> gauss;
		zee = std::complex<T>(gauss(gen), gauss(gen));
	}
};

using fftw_fixture = fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

BOOST_AUTO_TEST_CASE(fftw_3D) {
	using complex = std::complex<double>;  // TODO(correaa) make it work with thrust
	multi::array<complex, 3> in({10, 10, 10});
	in[2][3][4] = 99.0;
	multi::fftw::dft_forward(in);
	BOOST_REQUIRE(in[2][3][4] == 99.0);
}

BOOST_AUTO_TEST_CASE(fftw_1D_const) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const in = {1.0 + 2.0*I, 2.0 + 3.0 *I, 4.0 + 5.0*I, 5.0 + 6.0*I};

	auto fwd = multi::fftw::dft(in, fftw::forward);  // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2.0 - 2.0*I  );
	BOOST_REQUIRE( in[1]  == +2.0 + 3.0*I  );

	auto bwd = multi::fftw::dft(in, fftw::forward);  // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE( bwd[2] == -2.0 - 2.0*I  );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity_2, *boost::unit_test::tolerance(0.0001)) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft({false, false}, in, out, fftw::forward);  // out = in;
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity, *boost::unit_test::tolerance(0.0001)) {
	using complex = std::complex<double>;

	[[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{ 1.0 + 2.0 * I,  9.0 - 1.0 * I,  2.0 + 4.0 * I},
		{ 3.0 + 3.0 * I,  7.0 - 4.0 * I,  1.0 + 9.0 * I},
		{ 4.0 + 1.0 * I,  5.0 + 3.0 * I,  2.0 + 4.0 * I},
		{ 3.0 - 1.0 * I,  8.0 + 7.0 * I,  2.0 + 1.0 * I},
		{31.0 - 1.0 * I, 18.0 + 7.0 * I, 2.0 + 10.0 * I},
	};
	auto fwd = multi::fftw::dft({}, in, fftw::forward);
	BOOST_REQUIRE( fwd == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D, *boost::unit_test::tolerance(0.0001)) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};

	namespace fftw = multi::fftw;
	auto fwd = fftw::dft_forward(in);
	BOOST_TEST_REQUIRE( fwd[3][1].real() == -19.0455  );  // Fourier[in, FourierParameters -> {1, -1}][[4]][[2]]
	BOOST_TEST_REQUIRE( fwd[3][1].imag() == - 2.22717 );

	multi::array<complex, 1> const in0 = {1. + 2.*I, 9. - 1.*I, 2. + 4.*I};

	BOOST_REQUIRE( fftw::dft_forward(in[0]) == fftw::dft_forward(in0) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_rotated, *boost::unit_test::tolerance(0.0001)) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	using multi::array;
	array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	using multi::fftw::dft_forward;
	auto fwd = dft_forward(in);
	BOOST_REQUIRE(
		dft_forward(rotated(in)[0])
			== dft_forward(array<complex, 1>{1.0 + 2.0*I, 3.0 + 3.0*I, 4.0 + 1.0*I,  3.0 - 1.0*I, 31.0 - 1.0*I})
	);
	BOOST_REQUIRE( dft_forward(rotated(in)) == rotated(fwd) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_many, *boost::unit_test::tolerance(0.0001)) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	multi::array<complex, 2> out(extensions(in));

	using multi::fftw::dft_forward;

	multi::fftw::dft({fftw::none, fftw::forward}, in, out);
	BOOST_REQUIRE( dft_forward(in[0]) == out[0] );

	multi::fftw::dft({false, true}, rotated(in), rotated(out), fftw::forward);
	BOOST_REQUIRE( dft_forward(rotated(in)[0]) == rotated(out)[0] );

	multi::fftw::dft_forward({false, false}, rotated(in), rotated(out));
	BOOST_REQUIRE( in == out );

	multi::fftw::many_dft(in.begin(), in.end(), out.begin(), fftw::forward);
	BOOST_REQUIRE( dft_forward(in[0]) == out[0] );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_forward) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const in = {1.0 + 2.0*I, 2.0 + 3.0 *I, 4.0 + 5.0*I, 5.0 + 6.0*I};

	auto fwd = multi::fftw::dft_forward(in);  // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2.0 - 2.0*I  );
	BOOST_REQUIRE( in[1]  == +2.0 + 3.0*I  );

	auto bwd = multi::fftw::dft_forward(in);  // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE( bwd[2] == -2.0 - 2.0*I  );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_sign) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const in = {1.0 + 2.0*I, 2.0 + 3.0*I, 4.0 + 5.0*I, 5.0 + 6.0*I};

	auto const fwd = multi::fftw::dft(in, static_cast<multi::fftw::sign>(+1));  // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2. - 2.*I  );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_copy_by_false) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const in = {1.0 + 2.0*I, 2.0 + 3.0 *I, 4.0 + 5.0*I, 5.0 + 6.0*I};

	auto const out = multi::fftw::dft({false}, in, static_cast<multi::fftw::sign>(+1));
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_copy_by_false_forward) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto const out = multi::fftw::dft_forward({false}, in);
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_many1_from_2) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({3, 10}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({3, 10});
	fftw::dft({false, true}, in, out, fftw::forward);

	multi::array<complex, 2> out2({3, 10});
	std::transform(in.begin(), in.end(), out2.begin(), out2.begin(), [](auto const& in_elem, auto&& out2_elem) {
		fftw::dft_forward(in_elem, out2_elem);
		return std::forward<decltype(out2_elem)>(out2_elem);
	});

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_3) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 3> in ({3, 5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 3> out({3, 5, 6});
	fftw::dft_forward({false, true, true}, in, out);

	multi::array<complex, 3> out2({3, 5, 6});
	std::transform(in.begin(), in.end(), out2.begin(), out2.begin(), [](auto const& in_elem, auto&& out2_elem) {
		fftw::dft_forward(in_elem, out2_elem);
		return std::forward<decltype(out2_elem)>(out2_elem);
	});

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_2) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in ({5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({5, 6});
	fftw::dft({true, true}, in, out, static_cast<fftw::sign>(FFTW_FORWARD));

	multi::array<complex, 2> out2({5, 6});
	fftw::dft(in, out2, FFTW_FORWARD);
	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_4D) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 4> const in = [] {
		multi::array<complex, 4> in({6, 6, 6, 6}); in[2][3][4][5] = 99.0; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, true}, in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4][5] == 99.0);
}

BOOST_AUTO_TEST_CASE(fftw_4D_many) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	auto const in = [] {
		multi::array<complex, 4> in({7, 8, 9, 10}, {0.0, 0.0});
		in[2][3][4][5] = 99.0; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, false}, in, fftw::forward);
	BOOST_REQUIRE( in[2][3][4][5] == 99.0 );

	multi::array<complex, 4> out(extensions(in));
	multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
	BOOST_REQUIRE( out == fwd );
}

BOOST_AUTO_TEST_CASE(cufft_many_2D) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	auto const in = [] {
		multi::array<complex, 3> ret({10, 10, 10});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine{std::random_device{}()}, uniform_01 = std::uniform_real_distribution<>{}]() mutable{
				return complex{uniform_01(eng), uniform_01(eng)};
			}
		);
		return ret;
	}();
	multi::array<complex, 3> out(extensions(in));
	multi::fftw::many_dft((in.rotated()).begin(), (in.rotated()).end(), (out.rotated()).begin(), multi::fftw::forward);

	multi::array<complex, 3> out2(extensions(in));
	multi::fftw::dft_forward({true, false, true}, in, out2);
	BOOST_REQUIRE( out == out2 );
}

BOOST_AUTO_TEST_CASE(fftw_5D) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 5> in({4, 5, 6, 7, 8}, {0.0, 0.0});
	BOOST_REQUIRE( size(in) == 4 );

	in[2][3][4][5][6] = 99.0;
	auto const out_fwd = multi::fftw::dft(in, fftw::forward);

	BOOST_REQUIRE(in[2][3][4][5][6] == 99.0);
	BOOST_TEST_REQUIRE( power(in) - power(out_fwd)/num_elements(out_fwd) < 1e-5 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const pln{in, out, fftw::forward, fftw::preserve_input};
	pln();
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-7 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const pln{in.layout(), out.layout(), fftw::forward, fftw::preserve_input};
	pln(in.base(), out.base());
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern_measure) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const pln{in.layout(), out.layout(), fftw::forward, fftw::preserve_input};
	pln(in.base(), out.base());
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft_forward(in, out);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out_default) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 3> in({4, 4, 4}); std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);
	multi::array<complex, 3> const out = fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( std::abs(power(in) - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 3> io({4, 4, 4}); std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);
	auto powerin = power(io);
	fftw::dft_inplace(io, fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place_over_ref_inplace) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 3> io({4, 4, 4});
	std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto const powerin = power(io);
//  fftw::dft_inplace(multi::array_ref<complex, 3>(io.data(), io.extensions()), fftw::forward);
	fftw::dft_inplace(multi::array_ref<complex, 3>(data_elements(io), extensions(io)), fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_ref) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 3> in({4, 4, 4});
	std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);  //  NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 3> out({4, 4, 4});
	multi::array_ref<complex, 3>(data_elements(out), extensions(out)) = fftw::dft(multi::array_cref<complex, 3>(data_elements(in), extensions(in)), fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_temporary) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	double powerin = std::numeric_limits<double>::quiet_NaN();
	auto fun = [&]() {
		multi::array<complex, 3> in({4, 4, 4});
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		powerin = power(in);
		return in;
	};
	auto out = fftw::dft(fun(), fftw::forward);
	BOOST_REQUIRE( std::abs(powerin - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_transposition_square_inplace) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{ {11.0, 0.0}, {12.0, 0.0} },
		{ {21.0, 0.0}, {22.0, 0.0} }
	};
	BOOST_REQUIRE( in[1][0] == 21. );

	multi::fftw::copy(in, rotated(in));
	BOOST_TEST( in[0][1].real() == 21.0 );
	BOOST_TEST( in[0][1].imag() ==  0.0 );
}

BOOST_AUTO_TEST_CASE(fftw_4D_inq_poisson) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 4> const in = [] {
		multi::array<complex, 4> in({5, 10, 17, 1});
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		return in;
	}();

	multi::array<complex, 4> out(extensions(in));
	using multi::fftw::sign;
	multi::fftw::dft(
		{static_cast<sign>(0), static_cast<sign>(+1), static_cast<sign>(+1), static_cast<sign>(0)},
		in, out
	);

	using boost::multi::detail::get;
	BOOST_TEST( power(in) == power(out)/get<1>(sizes(out))/get<2>(sizes(out)) , boost::test_tools::tolerance(1e-10) );
}


BOOST_AUTO_TEST_CASE(fftw_1D_power_c_interface) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 1> in(16, {0.0, 0.0});
	BOOST_REQUIRE( size(in) == 16 );

	std::iota(begin(in), end(in), 1.0);
	BOOST_TEST_REQUIRE( power(in) == 1496.0 );

	multi::array<complex, 1> out(extensions(in));

	auto* pln = multi::fftw_plan_dft(in, out, fftw::forward, fftw::preserve_input);
	fftw_execute(pln);
	fftw_destroy_plan(pln);
	BOOST_TEST( power(in) == power(out)/num_elements(out), boost::test_tools::tolerance(1e-15) );
}

#if 0
BOOST_AUTO_TEST_CASE(fftw_2D_const_range_part1) {
	multi::array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};

	{
		multi::static_array<complex, 2> fwd(in.extensions());

		auto* data = fwd.data_elements();

		fwd = multi::fftw::fft(in);

		BOOST_REQUIRE( size(fwd) == size(in) );
		BOOST_REQUIRE( data == fwd.data_elements() );

		BOOST_TEST_REQUIRE( power(  in    ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
	{
		multi::array<complex, 2> fwd(in.extensions());

		auto* data = fwd.data_elements();

		fwd = multi::fftw::fft(in);

		BOOST_REQUIRE( size(fwd) == size(in) );
		BOOST_REQUIRE( data == fwd.data_elements() );

		BOOST_TEST_REQUIRE( power(  in    ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
	{
		multi::array<complex, 2> fwd;

		auto* data = fwd.data_elements();

		fwd = multi::fftw::fft(in);

		BOOST_REQUIRE( size(fwd) == size(in) );
		BOOST_REQUIRE( data != fwd.data_elements() );

		BOOST_TEST_REQUIRE( power(  in    ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_part2) {
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};

	{
		multi::array<complex, 2> fwd(multi::fftw::fft(in));

		BOOST_REQUIRE( size(fwd) == size(in) );

		BOOST_TEST_REQUIRE( power(in      ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
	{
		multi::array<complex, 2> fwd = multi::fftw::fft(in);

		BOOST_REQUIRE( size(fwd) == size(in) );

		BOOST_TEST_REQUIRE( power(in      ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
	{
		auto fwd = multi::array<complex, 2>(multi::fftw::fft(in));

		BOOST_REQUIRE( size(fwd) == size(in) );

		BOOST_TEST_REQUIRE( power(in      ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
	{
		auto fwd = + multi::fftw::fft(in);

		BOOST_REQUIRE( fwd.extensions() == in.extensions() );

		BOOST_TEST_REQUIRE( power(in      ) - power(  fwd    )/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[0]) - power((~fwd)[0])/size(fwd) < 1e-8 );
		BOOST_TEST_REQUIRE( power((~in)[2]) - power((~fwd)[2])/size(fwd) < 1e-8 );
	}
}
#endif

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};

	{
		multi::array<complex, 2> fwd(in.extensions());

		auto* data = fwd.data_elements();

		fwd = multi::fftw::ref(in);

		BOOST_REQUIRE( data == fwd.data_elements() );
		BOOST_REQUIRE( fwd == in );
	}
	{
		multi::array<complex, 2> const fwd = multi::fftw::ref(in);
		BOOST_REQUIRE( fwd == in );
	}
	{
		multi::array<complex, 2> fwd = in.transposed();
		BOOST_REQUIRE(   fwd .size() == 3 );
		BOOST_REQUIRE( (~fwd).size() == 5 );

		BOOST_TEST_REQUIRE( fwd[0][0] == in[0][0] );
		BOOST_TEST_REQUIRE( fwd[1][0] == in[0][1] );
		BOOST_TEST_REQUIRE( fwd[2][0] == in[0][2] );
	}
	{
		multi::array<complex, 2> fwd({3, 5}, {0.0, 0.0});
		fwd() = multi::fftw::ref(in.transposed());
		BOOST_REQUIRE(   fwd .size() == 3 );
		BOOST_REQUIRE( (~fwd).size() == 5 );

		BOOST_TEST_REQUIRE( fwd[0][0] == in[0][0] );
		BOOST_TEST_REQUIRE( fwd[1][0] == in[0][1] );
		BOOST_TEST_REQUIRE( fwd[2][0] == in[0][2] );
	}
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_part2) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> const in = {
		{  100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{    3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{    4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{    3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{   31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	{
		auto in2 = + multi::fftw::ref(in);
		BOOST_REQUIRE( in2 == in );
	}
	{
		auto in2 = + multi::fftw::ref(in.transposed());
		BOOST_REQUIRE( in2 == in.transposed() );
	}
	{
		auto in2 = + multi::fftw::ref(in.rotated());
		BOOST_REQUIRE( in2 == in.rotated() );
	}
	{
		auto&& ref = multi::fftw::ref(in);
		auto in2 =+ ref;
		BOOST_REQUIRE( in2 == in );

		multi::array<complex, 2> tt({3, 5});
		multi::fftw::dft({}, in, tt.transposed(), multi::fftw::forward);

		multi::array<complex, 2> in2t({3, 5}, {99.0, 0.0});
		in2t() = multi::fftw::ref(in).transposed();

		{
			std::for_each(in.begin(), in.end(), [](auto const& row) {
				std::copy(row.begin(), row.end(), std::ostream_iterator<complex>(std::cout, "\t"));
				std::cout<< std::endl;
			});
			std::cout<< std::endl;
		}
		{
			std::for_each(in2t.begin(), in2t.end(), [](auto const& row) {
				std::copy(row.begin(), row.end(), std::ostream_iterator<complex>(std::cout, "\t"));
				std::cout<< std::endl;
			});
			std::cout<< std::endl;
		}
		{
			std::for_each(tt.begin(), tt.end(), [](auto const& row) {
				std::copy(row.begin(), row.end(), std::ostream_iterator<complex>(std::cout, "\t"));
				std::cout<< std::endl;
			});
			std::cout<< std::endl;
		}

		BOOST_REQUIRE( tt   == in.transposed() );
		BOOST_REQUIRE( in2t == in.transposed() );
	}
	{
		auto in2 = + multi::fftw::ref(in).transposed();
		BOOST_REQUIRE( in2 == in.transposed() );
	}
	{
		multi::array<complex, 2> fwd({3, 5}, {0.0, 0.0});
		fwd() = multi::fftw::ref(in).transposed();
		BOOST_REQUIRE(   fwd .size() == 3 );
		BOOST_REQUIRE( (~fwd).size() == 5 );

		BOOST_TEST_REQUIRE( fwd[0][0] == in[0][0] );
		BOOST_TEST_REQUIRE( fwd[1][0] == in[0][1] );
		BOOST_TEST_REQUIRE( fwd[2][0] == in[0][2] );
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_many_new_interface) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	auto const in = [] {
		multi::array<complex, 4> in({17, 15, 10, 8}, {0.0, 0.0});
		in[2][3][4][5] = 99.0;
		return in;
	}();
	{
		auto fwd = + multi::fftw::ref(in)(fftw::forward, fftw::forward, fftw::forward, fftw::none);
		BOOST_REQUIRE( in[2][3][4][5] == 99.0 );

		multi::array<complex, 4> out(extensions(in));
		multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
		BOOST_REQUIRE( out == fwd );
	}
	{
		auto fwd = + multi::fftw::ref(in)(fftw::forward, fftw::forward, fftw::forward);
		BOOST_REQUIRE( in[2][3][4][5] == 99.0 );

		multi::array<complex, 4> out(extensions(in));
		multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
		BOOST_REQUIRE( out == fwd );
	}
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_naive) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{  100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{    3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{    4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{    3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{   31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	in = in.transposed();
//  BOOST_REQUIRE( in != in_transpose );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_naive_square) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	in = in.transposed();
	BOOST_REQUIRE( in != in_transpose );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	auto* in_base = in.base();
	in = multi::fftw::ref(in).transposed();

	BOOST_REQUIRE( in == in_transpose );
	BOOST_REQUIRE( in_base == in.base() );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_nested) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	auto* in_base = in.base();
	in = multi::fftw::ref(in.transposed());
	BOOST_REQUIRE( in == in_transpose );
	BOOST_REQUIRE( in_base == in.base() );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_square) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	auto* in_base = in.base();
	in = multi::fftw::ref(in).transposed();
	BOOST_REQUIRE( in == in_transpose );
	BOOST_REQUIRE( in_base == in.base() );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_square_nested) {
	using complex = std::complex<double>; [[maybe_unused]] auto const I = complex{0.0, 1.0};  // NOLINT(readability-identifier-length) imag unit

	multi::array<complex, 2> in = {
		{100.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I}
	};
	multi::array<complex, 2> const in_transpose = in.transposed();
	auto* in_base = in.base();
	in = multi::fftw::ref(in.transposed());
	BOOST_REQUIRE( in == in_transpose );
	BOOST_REQUIRE( in_base == in.base() );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_nonpod) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s
	multi::array<std::string, 2> in = {
		{  "100.0 + 2.0*I"s,  "9.0 - 1.0*I"s, "2.0 +  4.0*I"s},
		{    "3.0 + 3.0*I"s,  "7.0 - 4.0*I"s, "1.0 +  9.0*I"s},
		{    "4.0 + 1.0*I"s,  "5.0 + 3.0*I"s, "2.0 +  4.0*I"s},
		{    "3.0 - 1.0*I"s,  "8.0 + 7.0*I"s, "2.0 +  1.0*I"s},
		{   "31.0 - 1.0*I"s, "18.0 + 7.0*I"s, "2.0 + 10.0*I"s},
	};
	multi::array<std::string, 2> const in_transpose = in.transposed();
	in = in.transposed();
//  BOOST_REQUIRE( in != in_transpose );
}

BOOST_AUTO_TEST_CASE(fftw_2D_const_range_ref_transposed_nonpod_square) {
	multi::array<std::string, 2> in = {
		{  "100.0 + 2.0*I",  "9.0 - 1.0*I", "2.0 +  4.0*I"},  // std::string NOLINT(fuchsia-default-arguments-calls)
		{    "3.0 + 3.0*I",  "7.0 - 4.0*I", "1.0 +  9.0*I"},  // std::string NOLINT(fuchsia-default-arguments-calls)
		{    "4.0 + 1.0*I",  "5.0 + 3.0*I", "2.0 +  4.0*I"},  // std::string NOLINT(fuchsia-default-arguments-calls)
	};
	multi::array<std::string, 2> const in_transpose = in.transposed();
	in = in.transposed();
	BOOST_REQUIRE( in != in_transpose );
}
