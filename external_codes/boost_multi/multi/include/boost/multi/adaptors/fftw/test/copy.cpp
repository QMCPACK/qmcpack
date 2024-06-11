#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXXX $CXXFLAGS -O3 $0 -o $0x -DHAVE_FFTW3_THREADS -lfftw3 -lfftw3_threads -lboost_unit_test_framework -lboost_timer&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW copy"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include<boost/timer/timer.hpp>

#include "../../fftw.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(fftw_copy){

	using complex = std::complex<double>;
	auto const in = []{
		multi::array<complex, 4> ret({96, 96, 96, 96});
		std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(), 
			[](){return complex{std::rand()*1./RAND_MAX, std::rand()*1./RAND_MAX};}
		);
		return ret;
	}();
	std::cout<<"memory size "<< in.num_elements()*sizeof(complex)/1e6 <<" MB\n";
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			boost::timer::auto_cpu_timer t{"fftw_copy in-inorder %ws wall, CPU (%p%)\n"};
			multi::fftw::copy(in, rotated(out));
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( rotated(out) == in );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			boost::timer::auto_cpu_timer t{"fftw_copy out-inorder %ws wall, CPU (%p%)\n"};
			multi::fftw::copy(unrotated(in), out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( rotated(out) == in   );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			boost::timer::auto_cpu_timer t{"assignment in-inorder %ws wall, CPU (%p%)\n"};
			rotated(out) = in;
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			boost::timer::auto_cpu_timer t{"assignment out-inorder %ws wall, CPU (%p%)\n"};
			out = unrotated(in);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"assignment inplace out-inorder %ws wall, CPU (%p%)\n"};
			out = unrotated(out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"assignment inplace in-inorder %ws wall, CPU (%p%)\n"};
			rotated(out) = out;
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	//	BOOST_REQUIRE( rotated(out) == in );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"assignment inplace with copy out-inorder %ws wall, CPU (%p%)\n"};
			out = unrotated(multi::array<complex, 4>{out});
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"assignment inplace with copy in-inorder %ws wall, CPU (%p%)\n"};
			rotated(out) = multi::array<complex, 4>{out};
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"fftw copy inplace in-inorder %ws wall, CPU (%p%)\n"};
			multi::fftw::copy(out, rotated(out));
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			boost::timer::auto_cpu_timer t{"fftw copy inplace out-inorder %ws wall, CPU (%p%)\n"};
			multi::fftw::copy(unrotated(out), out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		auto p = out.data_elements();
		{
			boost::timer::auto_cpu_timer t{"fftw move construct inplace in-inorder %ws wall, CPU (%p%)\n"};
			multi::array<complex, 4> out2 = multi::fftw::copy( out.move().unrotated() );
			BOOST_REQUIRE( out.empty() );
			BOOST_REQUIRE( p == out2.data_elements() );
			BOOST_TEST( out2[1][2][3][4].real() == in[2][3][4][1].real() );
		}
	}
	{
		multi::array<complex, 4> out = in;
		auto p = out.data_elements();
		multi::array<complex, 4> out2;
		{
			boost::timer::auto_cpu_timer t{"fftw move assign inplace in-inorder %ws wall, CPU (%p%)\n"};
			out2 = multi::fftw::copy( out.move().unrotated() );
			BOOST_REQUIRE( out.empty() );
			BOOST_REQUIRE( p == out2.data_elements() );
			BOOST_TEST( out2[1][2][3][4].real() == in[2][3][4][1].real() );
		}
	}
	{
		multi::array<complex, 4> out = in;
		auto p = out.data_elements();
		{
			boost::timer::auto_cpu_timer t{"fftw move self-assign inplace in-inorder %ws wall, CPU (%p%)\n"};
			out = multi::fftw::copy( out.move().unrotated() );
			BOOST_REQUIRE( p == out.data_elements() );
			BOOST_TEST( out[1][2][3][4].real() == in[2][3][4][1].real() );
		}
	}
	

}

