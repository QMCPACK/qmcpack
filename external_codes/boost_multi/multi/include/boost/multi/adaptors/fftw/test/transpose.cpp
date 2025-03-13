// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/test/unit_test.hpp>

#include <boost/multi/adaptors/fftw.hpp>

#include <chrono>  // NOLINT(build/c++11)
#include <complex>
#include <iostream>
#include <random>

namespace multi = boost::multi;

class watch  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions)
: private std::chrono::high_resolution_clock {
	std::string label_;
	time_point  start_ = now();

 public:
	explicit watch(std::string label) : label_{std::move(label)} {}  // NOLINT(fuchsia-default-arguments-calls)

	watch(watch const&) = delete;

	auto operator=(watch const&) -> watch& = delete;

	auto elapsed_sec() const { return std::chrono::duration<double>(now() - start_).count(); }
	~watch() { std::cerr << label_ << ": " << elapsed_sec() << " sec" << '\n'; }  // NOLINT(cpp:S4963)
};

using fftw_fixture = multi::fftw::environment;
BOOST_TEST_GLOBAL_FIXTURE(fftw_fixture);

BOOST_AUTO_TEST_CASE(fftw_transpose) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	using complex = std::complex<double>;

	auto const in = std::invoke([] {
		multi::array<complex, 2> ret({101, 99});  // ({1013, 997});  // ({10137, 9973});
		std::generate(
			ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng = std::default_random_engine{std::random_device{}()}, uniform_01 = std::uniform_real_distribution<>{}]() mutable {
				return complex{uniform_01(eng), uniform_01(eng)};
			}
		);
		return ret;
	});

	multi::array<complex, 2> out = in;

	watch const unnamed{"transposition with aux   %ws wall, CPU (%p%)\n"s};

	multi::array<complex, 2> aux = ~out;

	out = std::move(aux);
	BOOST_REQUIRE( out[35][79] == in[79][35] );
}
