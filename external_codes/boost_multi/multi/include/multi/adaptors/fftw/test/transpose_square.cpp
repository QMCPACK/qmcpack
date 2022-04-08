#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS -Ofast $0 -o $0x -DHAVE_FFTW3_THREADS -lfftw3 -lfftw3_threads -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW transpose"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include<boost/timer/timer.hpp>

#include "../../fftw.hpp"

#include<complex>
#include<random>

namespace multi = boost::multi;

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

using complex = std::complex<double>;

BOOST_AUTO_TEST_CASE(fftw_transpose) {
	multi::fftw::initialize_threads();
	 {
		auto const in = [] {
		//	multi::array<complex, 2> ret({819, 819});
			multi::array<complex, 2> ret({81, 81});
			std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
				[eng = std::default_random_engine{std::random_device{}()}, uniform_01 = std::uniform_real_distribution<>{}]() mutable{
					return complex{uniform_01(eng), uniform_01(eng)};
				}
			);
		//	std::cout<<"memory size "<< ret.num_elements()*sizeof(complex)/1e6 <<" MB\n";
			return ret;
		}();
	//	multi::fftw::plan::with_nthreads(1);
		 {
			multi::array<complex, 2> out = in;
			auto* p = out.data_elements();
			 {
				boost::timer::auto_cpu_timer t{"fftw trans mve 1 thread  %ws wall, CPU (%p%)\n"};
				multi::fftw::transpose( out );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
//		{
//			multi::array<complex, 2> out = in;
//			auto p = out.data_elements();
//			{
//				boost::timer::auto_cpu_timer t{"fftw trans mve 1 thread  %ws wall, CPU (%p%)\n"};
//				out = multi::fftw::copy( transposed( move(out) ) );
//				BOOST_REQUIRE( out.data_elements() == p );
//				BOOST_REQUIRE( out[35][79] == in[79][35] );
//			}
//			BOOST_REQUIRE( out == ~in );
//		}
	//	multi::fftw::plan::with_nthreads(2);
//		{
//			multi::array<complex, 2> out = in;
//			auto p = out.data_elements();
//			{
//				boost::timer::auto_cpu_timer t{"fftw trans mve 2 thread  %ws wall, CPU (%p%)\n"};
//				out = multi::fftw::copy( ~move(out) );
//				BOOST_REQUIRE( out.data_elements() == p );
//				BOOST_REQUIRE( out[35][79] == in[79][35] );
//			}
//			BOOST_REQUIRE( out == ~in );
//		}
	//	multi::fftw::plan::with_nthreads(4);
//		{
//			multi::array<complex, 2> out = in;
//			auto p = out.data_elements();
//			{
//				boost::timer::auto_cpu_timer t{"fftw trans mve 4 thread  %ws wall, CPU (%p%)\n"};
//				out = multi::fftw::copy( ~move(out) );
//				BOOST_REQUIRE( out.data_elements() == p );
//				BOOST_REQUIRE( out[35][79] == in[79][35] );
//			}
//			BOOST_REQUIRE( out == ~in );
//		}
		{
			multi::array<complex, 2> out = in;
			multi::array<complex, 2> aux(extensions(out));
			 {
				boost::timer::auto_cpu_timer t{"auxiliary copy           %ws wall, CPU (%p%)\n"};
				aux = ~out;
				out = std::move(aux);
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			 {
				boost::timer::auto_cpu_timer t{"transposition with loop   %ws wall, CPU (%p%)\n"};
				for(auto i : extension(out)) {
					for(auto j = 0; j != i; ++j) {
						std::swap(out[i][j], out[j][i]);
					}
				}
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
		 {
			multi::array<complex, 2> out = in;
			 {
				boost::timer::auto_cpu_timer t{"transposition with loop 2 %ws wall, CPU (%p%)\n"};
				for(auto i = 0; i != out.size(); ++i) {
					for(auto j = i + 1; j != out.size(); ++j) {
						std::swap(out[i][j], out[j][i]);
					}
				}
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
	}
}

