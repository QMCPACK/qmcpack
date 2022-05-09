// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2020-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW transpose"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

//#include "../../../adaptors/../complex.hpp"
#include "../../../adaptors/fftw.hpp"
#include "../../../array.hpp"

#include<chrono>
#include<iostream>
#include<random>
//#include<thrust/complex.h> // TODO(correaa) make lib work with thrust complex

namespace {

namespace multi = boost::multi;
namespace fftw = multi::fftw;

using complex = std::complex<double>; MAYBE_UNUSED complex const I{0, 1};

template<class M> auto power(M const& m)->decltype(std::norm(m)) {return std::norm(m);}

template<class M, DELETE((M::rank_v < 1))> auto power(M const& m) {
	return accumulate(begin(m), end(m), 0., [](auto const& a, auto const& b){return a + power(b);});
}

struct sum_power{
	template<class A, class B> auto operator()(A const& a, B const& b) const{return a+power(b);}
};

MAYBE_UNUSED constexpr int N = 16;

} // end anonymous namespace

class watch : private std::chrono::high_resolution_clock{
	std::string label;
	time_point start = now();

 public:
	explicit watch(std::string label) : label{std::move(label)} {}
	watch(watch const&) = delete;
	watch(watch&&) = default;
	auto operator=(watch const&) = delete;
	auto operator=(watch&&) -> watch& = default; // NOLINT(fuchsia-trailing-return):
	~watch(){
		std::cerr<< label<<": "<< std::chrono::duration<double>(now() - start).count() <<" sec"<<std::endl;
	}
};

template<class T> struct randomizer{
	template<class M> void operator()(M&& m) const {for(auto&& e : m) {operator()(e);}}
	void operator()(T& e) const{ // NOLINT(runtime/references) : passing by reference
		static std::random_device r; static std::mt19937 g{r()}; static std::normal_distribution<T> d;
		e = d(g);
	}
};

template<class T> struct randomizer<std::complex<T>>{
	template<class M> void operator()(M&& m) const {for(auto&& e : m) {operator()(e);}}
	void operator()(std::complex<T>& e) const{ // NOLINT(runtime/references) : passing by reference
		static std::random_device r; static std::mt19937 g{r()}; static std::normal_distribution<T> d;
		e = std::complex<T>(d(g), d(g));
	}
};

struct fftw_fixture : fftw::environment{
//	void setup(){}
//	void teardown(){}//fftw_cleanup();}
};

BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

BOOST_AUTO_TEST_CASE(fftw_3D) {
	using complex = std::complex<double>;  // TODO(correaa) make it work with thrust
	multi::array<complex, 3> in({10, 10, 10});
	in[2][3][4] = 99.;
	multi::fftw::dft_forward(in);
	BOOST_REQUIRE(in[2][3][4] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_1D_const) {
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto fwd = multi::fftw::dft(in, fftw::forward); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2. - 2.*I  );
	BOOST_REQUIRE( in[1]  == +2. + 3.*I  );

	auto bwd = multi::fftw::dft(in, fftw::forward); // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE( bwd[2] == -2. - 2.*I  );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity_2, *boost::unit_test::tolerance(0.0001)) {
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft({false, false}, in, out, fftw::forward); // out = in;
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D_identity, *boost::unit_test::tolerance(0.0001)) {
	multi::array<complex, 2> const in = {
		{ 1. + 2.*I, 9. - 1.*I, 2. + 4.*I},
		{ 3. + 3.*I, 7. - 4.*I, 1. + 9.*I},
		{ 4. + 1.*I, 5. + 3.*I, 2. + 4.*I},
		{ 3. - 1.*I, 8. + 7.*I, 2. + 1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	auto fwd = multi::fftw::dft({}, in, fftw::forward);
	BOOST_REQUIRE( fwd == in );
}

BOOST_AUTO_TEST_CASE(fftw_2D, *boost::unit_test::tolerance(0.0001)) {
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};

	namespace fftw = multi::fftw;
	auto fwd = fftw::dft_forward(in);
	BOOST_TEST_REQUIRE( fwd[3][1].real() == -19.0455  ); // Fourier[in, FourierParameters -> {1, -1}][[4]][[2]]
	BOOST_TEST_REQUIRE( fwd[3][1].imag() == - 2.22717 );

	multi::array<complex, 1> const in0 = {1. + 2.*I, 9. - 1.*I, 2. + 4.*I};

	BOOST_REQUIRE( fftw::dft_forward(in[0]) == fftw::dft_forward(in0) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_rotated, *boost::unit_test::tolerance(0.0001)) {
	using multi::array;
	array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	using multi::fftw::dft_forward;
	auto fwd = dft_forward(in);
	BOOST_REQUIRE(
		dft_forward(rotated(in)[0])
			== dft_forward(array<complex, 1>{1.+2.*I, 3.+3.*I, 4. + 1.*I,  3. - 1.*I, 31. - 1.*I})
	);
	BOOST_REQUIRE( dft_forward(rotated(in)) == rotated(fwd) );
}

BOOST_AUTO_TEST_CASE(fftw_2D_many, *boost::unit_test::tolerance(0.0001)) {
	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
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
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto fwd = multi::fftw::dft_forward(in); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2. - 2.*I  );
	BOOST_REQUIRE( in[1]  == +2. + 3.*I  );

	auto bwd = multi::fftw::dft_forward(in); // InverseFourier[in, FourierParameters -> {-1, -1}]
	BOOST_REQUIRE( bwd[2] == -2. - 2.*I  );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_sign) {
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto const fwd = multi::fftw::dft(in, +1); // Fourier[in, FourierParameters -> {1, -1}]
	BOOST_REQUIRE( size(fwd) == size(in) );
	BOOST_REQUIRE( fwd[2] == -2. - 2.*I  );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_copy_by_false) {
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto const out = multi::fftw::dft({false}, in, +1);
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_1D_const_copy_by_false_forward) {
	multi::array<complex, 1> const in = {1. + 2.*I, 2. + 3. *I, 4. + 5.*I, 5. + 6.*I};

	auto const out = multi::fftw::dft_forward({false}, in);
	BOOST_REQUIRE( out == in );
}

BOOST_AUTO_TEST_CASE(fftw_many1_from_2) {
	multi::array<complex, 2> in({3, 10}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({3, 10});
	fftw::dft({false, true}, in, out, fftw::forward);

	multi::array<complex, 2> out2({3, 10});
	for(int i = 0; i!=size(in); ++i) {
		fftw::dft_forward(in[i], out2[i]);
	}

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_3) {
	multi::array<complex, 3> in({3, 5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 3> out({3, 5, 6});
	fftw::dft_forward({false, true, true}, in, out);

	multi::array<complex, 3> out2({3, 5, 6});
	for(int i = 0; i!=size(in); ++i) {
		fftw::dft_forward(in[i], out2[i]);
	}

	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_many2_from_2) {
	multi::array<complex, 2> in({5, 6}); randomizer<complex>{}(in);
	multi::array<complex, 2> out({5, 6});
	fftw::dft({true, true}, in, out, FFTW_FORWARD);

	multi::array<complex, 2> out2({5, 6});
	fftw::dft(in, out2, FFTW_FORWARD);
	BOOST_REQUIRE(out2 == out);
}

BOOST_AUTO_TEST_CASE(fftw_4D) {
	multi::array<complex, 4> const in = [] {
		multi::array<complex, 4> in({10, 10, 10, 10}); in[2][3][4][5] = 99.; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, true}, in, fftw::forward);
	BOOST_REQUIRE(in[2][3][4][5] == 99.);
}

BOOST_AUTO_TEST_CASE(fftw_4D_many) {
	auto const in = [] {
		multi::array<complex, 4> in({97, 95, 101, 10}, 0.);
		in[2][3][4][5] = 99.; return in;
	}();
	auto fwd = multi::fftw::dft({true, true, true, false}, in, fftw::forward);
	BOOST_REQUIRE( in[2][3][4][5] == 99. );

	multi::array<complex, 4> out(extensions(in));
	multi::fftw::many_dft(begin(unrotated(in)), end(unrotated(in)), begin(unrotated(out)), fftw::forward);
	BOOST_REQUIRE( out == fwd );
}

BOOST_AUTO_TEST_CASE(cufft_many_2D) {
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
	multi::array<complex, 5> in({4, 5, 6, 7, 8}, 0.);
	BOOST_REQUIRE( size(in) == 4 );

	in[2][3][4][5][6] = 99.;
	auto const out_fwd = multi::fftw::dft(in, fftw::forward);

	BOOST_REQUIRE(in[2][3][4][5][6] == 99.);
	BOOST_REQUIRE( power(in) - power(out_fwd)/num_elements(out_fwd) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const p{in, out, fftw::forward, fftw::preserve_input};
	p(); //execute(p); //p.execute();
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const p{in.layout(), out.layout(), fftw::forward, fftw::preserve_input};
	p(in.base(), out.base()); //execute(p); //p.execute();
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_plan_modern_measure) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::plan const p{in.layout(), out.layout(), fftw::forward, fftw::preserve_input};
	p(in.base(), out.base()); //execute(p); //p.execute();
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 2> out(extensions(in));
	multi::fftw::dft_forward(in, out);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_power_dft_out_default) {
	multi::array<complex, 2> in({16, 16});
	std::iota(data_elements(in), data_elements(in) + num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto out = multi::fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-8 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power) {
	multi::array<complex, 3> in({4, 4, 4}); std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);
	multi::array<complex, 3> out = fftw::dft(in, fftw::forward);
	BOOST_REQUIRE( std::abs(power(in) - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place) {
	multi::array<complex, 3> io({4, 4, 4}); std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2);
	auto powerin = power(io);
	fftw::dft_inplace(io, fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_in_place_over_ref_inplace) {
	multi::array<complex, 3> io({4, 4, 4});
	std::iota(io.data_elements(), io.data_elements() + io.num_elements(), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	auto const powerin = power(io);
//	fftw::dft_inplace(multi::array_ref<complex, 3>(io.data(), io.extensions()), fftw::forward);
	fftw::dft_inplace(multi::array_ref<complex, 3>(data_elements(io), extensions(io)), fftw::forward);
	BOOST_REQUIRE( powerin - power(io)/num_elements(io) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_ref) {
	multi::array<complex, 3> in({4, 4, 4});
	std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2); //  NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
	multi::array<complex, 3> out({4, 4, 4});
	multi::array_ref<complex, 3>(data_elements(out), extensions(out)) = fftw::dft(multi::array_cref<complex, 3>(data_elements(in), extensions(in)), fftw::forward);
	BOOST_REQUIRE( power(in) - power(out)/num_elements(out) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_3D_power_out_of_place_over_temporary) {
	double powerin = NAN;
	auto f = [&](){
		multi::array<complex, 3> in({4, 4, 4});
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		powerin = power(in);
		return in;
	};
	auto out = fftw::dft(f(), fftw::forward);
	BOOST_REQUIRE( std::abs(powerin - power(out)/num_elements(out)) < 1e-10 );
}

BOOST_AUTO_TEST_CASE(fftw_2D_transposition_square_inplace) {
	multi::array<complex, 2> in = {
		{11., 12.},
		{21., 22.}
	};
	BOOST_REQUIRE( in[1][0] == 21. );

	multi::fftw::copy(in, rotated(in));
	BOOST_TEST( in[0][1].real() == 21. );
	BOOST_TEST( in[0][1].imag() ==  0. );
}

BOOST_AUTO_TEST_CASE(fftw_4D_inq_poisson) {
	multi::array<complex, 4> const in = [] {
		multi::array<complex, 4> in({50, 100, 137, 1});
		std::iota(data_elements(in), data_elements(in)+num_elements(in), 1.2); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): test code
		return in;
	}();

	multi::array<complex, 4> out(extensions(in));
	multi::fftw::dft({0, 1, 1, 0}, in, out);

	using boost::multi::detail::get;
	BOOST_TEST( power(in) == power(out)/get<1>(sizes(out))/get<2>(sizes(out)) , boost::test_tools::tolerance(1e-10) );
}

BOOST_AUTO_TEST_CASE(fftw_1D_power) {
	multi::array<complex, 1> in(N, 0.);
	BOOST_REQUIRE( size(in) == N );

	std::iota(begin(in), end(in), 1.);
	BOOST_TEST_REQUIRE( power(in) == 1496. );

	multi::array<complex, 1> out(extensions(in));

	auto* p = multi::fftw_plan_dft(in, out, fftw::forward, fftw::preserve_input);
	fftw_execute(p);
	fftw_destroy_plan(p);
	BOOST_TEST( power(in) == power(out)/num_elements(out), boost::test_tools::tolerance(1e-15) );
}

