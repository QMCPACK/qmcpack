// Copyright 2020-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

sass

#include "../../fftw.hpp"

#include<chrono>
#include<complex>
#include<iostream>
#include<string>
#include<utility>  // move

namespace multi = boost::multi;

namespace {
// minimal RAII wall-clock timer (replaces boost::timer::auto_cpu_timer)
class auto_timer {
	std::string                           label_;
	std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();

 public:
	explicit auto_timer(std::string label = {}) : label_{std::move(label)} {}
	auto_timer(auto_timer const&)                    = delete;
	auto_timer(auto_timer&&)                         = delete;
	auto operator=(auto_timer const&) -> auto_timer& = delete;
	auto operator=(auto_timer&&) -> auto_timer&      = delete;
	~auto_timer() {
		std::cerr << label_ << std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count() << " s (wall)\n";
	}
};
}  // namespace

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
			auto_timer t{"fftw_copy in-inorder: "};
			multi::fftw::copy(in, rotated(out));
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( rotated(out) == in );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			auto_timer t{"fftw_copy out-inorder: "};
			multi::fftw::copy(unrotated(in), out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( rotated(out) == in   );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			auto_timer t{"assignment in-inorder: "};
			rotated(out) = in;
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	}
	{
		multi::array<complex, 4> out(extensions(in), 0.);
		{
			auto_timer t{"assignment out-inorder: "};
			out = unrotated(in);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"assignment inplace out-inorder: "};
			out = unrotated(out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"assignment inplace in-inorder: "};
			rotated(out) = out;
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
	//  BOOST_REQUIRE( rotated(out) == in );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"assignment inplace with copy out-inorder: "};
			out = unrotated(multi::array<complex, 4>{out});
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"assignment inplace with copy in-inorder: "};
			rotated(out) = multi::array<complex, 4>{out};
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"fftw copy inplace in-inorder: "};
			multi::fftw::copy(out, rotated(out));
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		{
			auto_timer t{"fftw copy inplace out-inorder: "};
			multi::fftw::copy(unrotated(out), out);
		}
		BOOST_REQUIRE( out[1][2][3][4] == in[2][3][4][1] );
		BOOST_REQUIRE( out == unrotated(in) );
	}
	{
		multi::array<complex, 4> out = in;
		auto p = out.data_elements();
		{
			auto_timer t{"fftw move construct inplace in-inorder: "};
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
			auto_timer t{"fftw move assign inplace in-inorder: "};
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
			auto_timer t{"fftw move self-assign inplace in-inorder: "};
			out = multi::fftw::copy( out.move().unrotated() );
			BOOST_REQUIRE( p == out.data_elements() );
			BOOST_TEST( out[1][2][3][4].real() == in[2][3][4][1].real() );
		}
	}
	

}

