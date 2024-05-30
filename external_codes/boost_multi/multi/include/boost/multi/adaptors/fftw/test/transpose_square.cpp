// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <boost/multi/adaptors/fftw.hpp>

#include <chrono>  // NOLINT(build/c++11)
#include <complex>
#include <iostream>
#include <random>

namespace multi = boost::multi;

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE(fftw_fixture);

using complex = std::complex<double>;

class watch : private std::chrono::high_resolution_clock {  // NOSONAR(cpp:S4963) this class will report timing on destruction
	std::string label_;
	time_point  start_ = now();

 public:
	explicit watch(std::string label) : label_{std::move(label)} {}

	watch(watch const&) = delete;
	watch(watch&&)      = delete;

	auto operator=(watch const&) = delete;
	auto operator=(watch&&)      = delete;

	auto elapsed_sec() const { return std::chrono::duration<double>(now() - start_).count(); }
	~watch() { std::cerr << label_ << ": " << elapsed_sec() << " sec" << '\n'; }
};

BOOST_AUTO_TEST_CASE(fftw_transpose) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	multi::fftw::initialize_threads();
	{
		auto const in = [] {
			//  multi::array<complex, 2> ret({819, 819});
			multi::array<complex, 2> ret({81, 81});
			std::generate(
				ret.data_elements(), ret.data_elements() + ret.num_elements(),
				[eng        = std::default_random_engine{std::random_device{}()},
				 uniform_01 = std::uniform_real_distribution<>{}]() mutable {
					return complex{uniform_01(eng), uniform_01(eng)};
				}
			);
			//  std::cout<<"memory size "<< ret.num_elements()*sizeof(complex)/1e6 <<" MB\n";
			return ret;
		}();
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
