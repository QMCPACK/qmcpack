// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW transpose"
#include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"

#include<chrono>
#include<complex>
#include<random>

namespace multi = boost::multi;

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
	~watch() {std::cerr<< label <<": "<< elapsed_sec() <<" sec"<<std::endl;}
};

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

BOOST_AUTO_TEST_CASE(fftw_transpose) {
	using complex = std::complex<double>;

	 {
		auto const in = [] {
		//	multi::array<complex, 2> ret({10137, 9973});
		//	multi::array<complex, 2> ret({1013, 997});
			multi::array<complex, 2> ret({101, 99});
			std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
				[eng = std::default_random_engine{std::random_device{}()}, uniform_01 = std::uniform_real_distribution<>{}]() mutable{
					return complex{uniform_01(eng), uniform_01(eng)};
				}
			);
		//	std::cout<<"memory size "<< ret.num_elements()*sizeof(complex)/1e6 <<" MB\n";
			return ret;
		}();

		{
			multi::array<complex, 2> out = in;
			 {
				watch unnamed{"transposition with aux   %ws wall, CPU (%p%)\n"};
				multi::array<complex, 2> aux = ~out;
				out = std::move(aux);
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
		 {
			multi::array<complex, 2> out = in;
			auto* out_data = out.data_elements();
			 {
				watch unnamed{"fftw transpose fun thread  %ws wall, CPU (%p%)\n"};
				multi::fftw::transpose( out );
				BOOST_REQUIRE( out.data_elements() == out_data );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
	}
}
