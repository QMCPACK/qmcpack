// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW transpose"

#include<boost/test/unit_test.hpp>
#include<boost/timer/timer.hpp>

#include <multi/adaptors/fftw.hpp>

#include <chrono>  // NOLINT(build/c++11)
#include <complex>
#include <iostream>
#include <random>

namespace multi = boost::multi;

using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE( fftw_fixture );

using complex = std::complex<double>;

class watch : private std::chrono::high_resolution_clock {
	std::string label;
	time_point start = now();

 public:
	template<class String>
	explicit watch(String&& label) : label{std::forward<String>(label)} {}  // NOLINT(fuchsia-default-arguments-calls)
	watch(watch const&) = delete;
	watch(watch&&) = delete;

	auto operator=(watch const&) = delete;
	auto operator=(watch&&) = delete;

	auto elapsed_sec() const {return std::chrono::duration<double>(now() - start).count();}
	~watch() { std::cerr<< label <<": "<< elapsed_sec() <<" sec"<<std::endl; }
};

BOOST_AUTO_TEST_CASE(fftw_transpose) {
	multi::fftw::initialize_threads();
	{
		auto const in = [] {
		//  multi::array<complex, 2> ret({819, 819});
			multi::array<complex, 2> ret({81, 81});
			std::generate(ret.data_elements(), ret.data_elements() + ret.num_elements(),
				[eng = std::default_random_engine{std::random_device{}()}, uniform_01 = std::uniform_real_distribution<>{}]() mutable{
					return complex{uniform_01(eng), uniform_01(eng)};
				}
			);
		//  std::cout<<"memory size "<< ret.num_elements()*sizeof(complex)/1e6 <<" MB\n";
			return ret;
		}();
	//  multi::fftw::plan::with_nthreads(1);
		// {
		//  multi::array<complex, 2> out = in;
		//  auto* data = out.data_elements();
		//   {
		//      watch const unnamed{"fftw trans mve 1 thread  %ws wall, CPU (%p%)\n"s};
		//      multi::fftw::transpose( out );
		//      BOOST_REQUIRE( out.data_elements() == data );
		//      BOOST_REQUIRE( out[35][79] == in[79][35] );
		//  }
		//  BOOST_REQUIRE( out == ~in );
		// }
//    {
//      multi::array<complex, 2> out = in;
//      auto p = out.data_elements();
//      {
//        boost::timer::auto_cpu_timer t{"fftw trans mve 1 thread  %ws wall, CPU (%p%)\n"};
//        out = multi::fftw::copy( transposed( move(out) ) );
//        BOOST_REQUIRE( out.data_elements() == p );
//        BOOST_REQUIRE( out[35][79] == in[79][35] );
//      }
//      BOOST_REQUIRE( out == ~in );
//    }
//    multi::fftw::plan::with_nthreads(2);
//    {
//      multi::array<complex, 2> out = in;
//      auto p = out.data_elements();
//      {
//        boost::timer::auto_cpu_timer t{"fftw trans mve 2 thread  %ws wall, CPU (%p%)\n"};
//        out = multi::fftw::copy( ~move(out) );
//        BOOST_REQUIRE( out.data_elements() == p );
//        BOOST_REQUIRE( out[35][79] == in[79][35] );
//      }
//      BOOST_REQUIRE( out == ~in );
//    }
//    multi::fftw::plan::with_nthreads(4);
//    {
//      multi::array<complex, 2> out = in;
//      auto p = out.data_elements();
//      {
//        boost::timer::auto_cpu_timer t{"fftw trans mve 4 thread  %ws wall, CPU (%p%)\n"};
//        out = multi::fftw::copy( ~move(out) );
//        BOOST_REQUIRE( out.data_elements() == p );
//        BOOST_REQUIRE( out[35][79] == in[79][35] );
//      }
//      BOOST_REQUIRE( out == ~in );
//    }
		{
			multi::array<complex, 2> out = in;
			multi::array<complex, 2> aux(extensions(out));
			{
				watch const unnamed{"auxiliary copy           %ws wall, CPU (%p%)\n"s};
				aux = ~out;
				out = std::move(aux);
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			{
				watch const unnamed{"transposition with loop   %ws wall, CPU (%p%)\n"s};
				std::for_each(extension(out).begin(), extension(out).end(), [&out](auto idx) {
					auto ext = multi::extension_t(0L, idx);
					std::for_each(ext.begin(), ext.end(), [&out, idx](auto jdx) {
						std::swap(out[idx][jdx], out[jdx][idx]);
					});
				});
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			{
				watch const unnamed{"transposition with loop 2 %ws wall, CPU (%p%)\n"s};
				std::for_each(extension(out).begin(), extension(out).end(), [&out](auto idx) {
					auto ext = multi::extension_t(idx + 1, out.size());
					std::for_each(ext.begin(), ext.end(), [&out, idx](auto jdx) {
						std::swap(out[idx][jdx], out[jdx][idx]);
					});
				});
				BOOST_REQUIRE( out[35][79] == in[79][35] );
			}
			BOOST_REQUIRE( out == ~in );
		}
	}
}
