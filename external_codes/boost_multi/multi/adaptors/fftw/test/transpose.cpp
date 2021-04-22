#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS -O3 $0 -o $0x -DHAVE_FFTW3_THREADS -lfftw3 -lfftw3_threads -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW transpose"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include<boost/timer/timer.hpp>

#include "../../fftw.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(fftw_transpose){

//	multi::fftw::initialize_threads();
	multi::fftw::plan::with_nthreads(1);

	using complex = std::complex<double>;


	{
		auto const in = []{
			multi::array<complex, 2> ret({10137, 9973});
			std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
				[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
			);
			std::cout<<"memory size "<< ret.num_elements()*sizeof(complex)/1e6 <<" MB\n";
			return ret;
		}();
		
		{
			multi::array<complex, 2> out = in;
			{
				boost::timer::auto_cpu_timer t{"transposition with aux   %ws wall, CPU (%p%)\n"};
				multi::array<complex, 2> aux = ~out;
				out = std::move(aux);
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw trans mve 1 thread  %ws wall, CPU (%p%)\n"};
				multi::array<complex, 2> out2 = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out2.data_elements() == p );
				BOOST_REQUIRE( out2[35][79] == in[79][35] );
			}
		}
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw transpose fun thread  %ws wall, CPU (%p%)\n"};
				multi::fftw::transpose( out );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw transpose 1 thread  %ws wall, CPU (%p%)\n"};
				out = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
		multi::fftw::plan::with_nthreads(2);
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw trans mve 2 thread  %ws wall, CPU (%p%)\n"};
				multi::array<complex, 2> out2 = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out2.data_elements() == p );
				BOOST_REQUIRE( out2[35][79] == in[79][35] );
			}
		}
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw transpose 2 threads %ws wall, CPU (%p%)\n"};
				out = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
		multi::fftw::plan::with_nthreads(3);
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw transpose 3 threads %ws wall, CPU (%p%)\n"};
				out = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
		multi::fftw::plan::with_nthreads(4);
		{
			multi::array<complex, 2> out = in;
			auto p = out.data_elements();
			{
				boost::timer::auto_cpu_timer t{"fftw transpose 4 threads %ws wall, CPU (%p%)\n"};
				out = multi::fftw::copy( ~move(out) );
				BOOST_REQUIRE( out.data_elements() == p );
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
		}
	}
}

