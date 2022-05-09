// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2020-2021
// Copyright 2020-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"

#include <boost/timer/timer.hpp>

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

	std::vector<std::array<bool, 4>> cases = {
		{false, true , true , true },
		{false, true , true , false},
		{true , false, false, false},
		{true , true , false, false},
		{false, false, true , false},
		{false, false, false, false},
	};

	using std::cout;
	for(auto c : cases) {
		cout<<"case ";
		copy(begin(c), end(c), std::ostream_iterator<bool>{cout, ", "});
		cout<<"\n";

		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"cpu_oplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in, out);
		}
		{
			multi::fftw::plan p{c, in, out, multi::fftw::forward};
			boost::timer::auto_cpu_timer t{"cpu_oplac planned %ws wall, CPU (%p%)\n"};
			p();
		}
		{
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_iplac %ws wall, CPU (%p%)\n"};
			multi::fftw::dft_forward(c, in_rw);
		}
		{
			auto in_rw = in;
			multi::fftw::plan p{c, in_rw, in_rw, multi::fftw::forward};
			boost::timer::auto_cpu_timer t{"cpu_iplac planned %ws wall, CPU (%p%)\n"};
			p();
		}
		{
			auto in_rw = in;
			multi::fftw::plan p{c, in_rw, in_rw, multi::fftw::forward};
			boost::timer::auto_cpu_timer t{
				"cpu_iplac planned measured %ws wall, CPU (%p%)\n"
			};
			p();
		}
		{
			boost::timer::auto_cpu_timer t{"cpu_alloc %ws wall, CPU (%p%)\n"};
			auto out_cpy = multi::fftw::dft_forward(c, in);
			BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		}
		 {
			auto in_rw = in;
			boost::timer::auto_cpu_timer t{"cpu_move %ws wall, CPU (%p%)\n"};
			auto out_cpy = multi::fftw::dft_forward(c, std::move(in_rw));
			BOOST_TEST(abs(out_cpy[5][4][3][1] - out[5][4][3][1]) == 0.);
		}
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark, *boost::unit_test::enabled() ) {
	using complex = std::complex<double>;
	namespace fftw = multi::fftw;

	auto x = multi::array<complex, 4>::extensions_type({6, 12, 12, 12});
	multi::array<complex, 4> in(x);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_REQUIRE(in[0][0][0][0] == 1.2);
	std::array<bool, 4> c = {false, true, true, true};
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}] {
		fftw::dft(c, in, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" inplace FTTT"}] {
		fftw::dft(c, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_REQUIRE(in0000 != 1.2);


	multi::array<complex, 4> out(x);
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace FTTT"}] {
		fftw::dft(c, in, out, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}] {
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	[&, _ = watch{utf::current_test_case().full_name()+" outofplace+alloc FTTT"}] {
		multi::array<complex, 4> out2(x);
		fftw::dft(c, in, out2, fftw::forward);
	}();
	BOOST_REQUIRE(in0000 == in[0][0][0][0]);
}
