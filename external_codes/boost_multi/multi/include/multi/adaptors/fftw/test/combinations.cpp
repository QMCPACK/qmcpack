// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor"
#include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"

#include<chrono>
#include<complex>
#include<iostream>
#include<random>

namespace multi = boost::multi;

namespace utf = boost::unit_test::framework;

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

class watch : private std::chrono::high_resolution_clock {
	std::string label;
	time_point start = now();

 public:
	explicit watch(std::string label) : label{std::move(label)} {}
	watch(watch const&) = delete;
	watch(watch&&) = default;
	auto operator=(watch const&) = delete;
	auto operator=(watch&&) -> watch& = default;
	auto elapsed_sec() const {return std::chrono::duration<double>(now() - start).count();}
	~watch() {
		std::cerr
			<< label <<": "
			<< elapsed_sec() <<" sec"
			<<std::endl
		;
	}
};

BOOST_AUTO_TEST_CASE(fft_combinations, *boost::unit_test::tolerance(0.00001) ) {

	using complex = std::complex<double>;

	auto const in = [] {
		multi::array<complex, 4> ret({10 , 11 , 12 , 13 });
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine {std::random_device {}()},
				uniform_01 = std::uniform_real_distribution<>{}]() mutable{
				return complex{uniform_01(eng), uniform_01(eng)};
			});
		return ret;
	}();

	std::vector<std::array<bool, 4>> which_cases = {
		{false, true , true , true },
		{false, true , true , false},
		{true , false, false, false},
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	using std::cout;
	for(auto which : which_cases) {
		cout<<"case ";
		copy(begin(which), end(which), std::ostream_iterator<bool>{cout, ", "});
		cout<<"\n";

		multi::array<complex, 4> out = in;
		{
			watch unnamed{"cpu_oplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(which, in, out);
		}
		{
			multi::fftw::plan pln{which, in, out, multi::fftw::forward};
			watch unnamed{"cpu_oplac planned %ws wall, CPU (%p%)\n"};
			pln();
		}
		{
			auto in_rw = in;
			watch unnamed{"cpu_iplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(which, in_rw);
		}
		{
			auto in_rw = in;
			multi::fftw::plan pln{which, in_rw, in_rw, multi::fftw::forward};
			watch unnamed{"cpu_iplac planned %ws wall, CPU (%p%)\n"};
			pln();
		}
		{
			auto in_rw = in;
			multi::fftw::plan pln{which, in_rw, in_rw, multi::fftw::forward};
			watch unnamed{"cpu_iplac planned measured %ws wall, CPU (%p%)\n"};
			pln();
		}
		{
			watch unnamed{"cpu_alloc %ws wall, CPU (%p%)\n"};
			auto out_cpy = multi::fftw::dft_forward(which, in);
			BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		}
		{
			auto in_rw = in;
			watch unnamed{"cpu_move %ws wall, CPU (%p%)\n"};
			auto out_cpy = multi::fftw::dft_forward(which, std::move(in_rw));
			BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		}
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark, *boost::unit_test::enabled() ) {
	using complex = std::complex<double>;
	namespace fftw = multi::fftw;

	auto exts = multi::array<complex, 4>::extensions_type({6, 12, 12, 12});
	multi::array<complex, 4> in(exts);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_REQUIRE(in[0][0][0][0] == 1.2);
	std::array<bool, 4> which = {false, true, true, true};
	[&, unnamed = watch{utf::current_test_case().full_name()+" inplace FTTT"}] {
		fftw::dft(which, in, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" inplace FTTT"}] {
		fftw::dft(which, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_REQUIRE(in0000 != 1.2);

	multi::array<complex, 4> out(exts);
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}] {
		multi::array<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	[&, unnamed = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}] {
		multi::array<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	BOOST_REQUIRE(in0000 == in[0][0][0][0]);
}


BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark_syntax) {
	std::vector<std::array<bool, 4>> which_cases = {
		{false, true , true , true },
		{false, true , true , false},
		{true , false, false, false},
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};
	using complex = std::complex<double>;

	auto const in = [] {
		multi::array<complex, 4> ret({6, 12, 12, 12});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine {std::random_device {}()},
				uniform_01 = std::uniform_real_distribution<>{}]() mutable{
				return complex{uniform_01(eng), uniform_01(eng)};
			});
		return ret;
	}();

	auto io = in; (void)io;
	BOOST_REQUIRE( io.extensions() == in.extensions() );

	namespace fftw = multi::fftw;
	using clock = std::chrono::high_resolution_clock;
	{
		auto const tick = clock::now();
		multi::array<complex, 4> out({6, 12, 12, 12});
		out = multi::fftw::ref(in)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
		BOOST_REQUIRE( out.extensions() == in.extensions() );
		auto time = std::chrono::duration<double>(clock::now() - tick);
		std::cout<<"allocate and copy assign (out-of-place fft) : "<< time.count() <<std::endl;
	}
	{
		auto const tick = clock::now();
		auto const out = +multi::fftw::ref(in)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
		BOOST_REQUIRE( out.extensions() == in.extensions() );
		auto time = std::chrono::duration<double>(clock::now() - tick);
		std::cout<<"copy construct (out-of-place fft) : "<< time.count() <<std::endl;
	}
	{
		auto const tick = clock::now();
		io = multi::fftw::ref(io)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
		BOOST_REQUIRE( io.extensions() == in.extensions() );
		auto time = std::chrono::duration<double>(clock::now() - tick);
		std::cout<<"self copy assign (in-place fft) : "<< time.count() <<std::endl;
	}
	{
		auto const tick = clock::now();
		multi::array<complex, 4> out = multi::fftw::move(io)(fftw::none, fftw::forward, fftw::forward, fftw::forward);
		BOOST_REQUIRE( io.is_empty() );
		auto time = std::chrono::duration<double>(clock::now() - tick);
		std::cout<<"move construct (in-place fft) : "<< time.count() <<std::endl;
	}
}
