// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/array.hpp>

#include <algorithm>   // for generate
#include <chrono>      // for operator-, duration, system...  // NOLINT(build/c++11)
#include <complex>     // for operator==, complex
#include <functional>  // for invoke  // IWYU pragma: keep
#include <iostream>    // for operator<<, basic_os...
#include <random>      // for linear_congruential_...
#include <string>      // for operator<<, operator""s
#include <utility>     // for move

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

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	multi::fftw::environment const env;

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

	watch const unnamed{"transposition with aux   %ws wall, CPU (%p%)\n"s};  //  NOLINT(misc-include-cleaner) bug in clang-tidy 18

	multi::array<complex, 2> aux{~out};

	out = std::move(aux);
	BOOST_TEST( out[35][79] == in[79][35] );

	return boost::report_errors();
}
