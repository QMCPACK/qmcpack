// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/fftw.hpp>  // for initialize_threads, environ...
#include <boost/multi/array.hpp>          // for array, subarray, layout_t

#include <algorithm>   // for for_each, generate
#include <chrono>      // for duration, operator-, high_r...  // NOLINT(build/c++11)
#include <complex>     // for operator==, complex
#include <functional>  // for invoke  // IWYU pragma: keep
#include <iostream>    // for basic_ostream, operator<<
#include <random>      // for uniform_real_distribution
#include <string>      // for char_traits, operator""s
#include <utility>     // for move, swap

namespace multi = boost::multi;

using complex = std::complex<double>;

class watch : private std::chrono::high_resolution_clock {  // NOSONAR(cpp:S4963) this class will report timing on destruction
	std::string label_;
	time_point  start_ = now();

 public:
	explicit watch(std::string label) : label_{ std::move(label) } {}

	watch(watch const&) = delete;
	watch(watch&&)      = delete;

	auto operator=(watch const&) = delete;
	auto operator=(watch&&)      = delete;

	auto elapsed_sec() const { return std::chrono::duration<double>(now() - start_).count(); }
	~watch() { std::cerr << label_ << ": " << elapsed_sec() << " sec" << '\n'; }
};

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
multi::fftw::environment const env;
BOOST_AUTO_TEST_CASE(fftw_transpose) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	multi::fftw::initialize_threads();
	{
		auto const in = std::invoke([] {
			//  multi::array<complex, 2> ret({819, 819});
			multi::array<complex, 2> ret({ 81, 81 });
			std::generate(
				ret.data_elements(), ret.data_elements() + ret.num_elements(),
				[eng        = std::default_random_engine{ std::random_device{}() },
				 uniform_01 = std::uniform_real_distribution<>{}]() mutable {
					return complex{ uniform_01(eng), uniform_01(eng) };
				}
			);
			return ret;
		});
		{
			multi::array<complex, 2> out = in;
			multi::array<complex, 2> aux(out.extensions());
			{
				watch const unnamed("auxiliary copy           %ws wall, CPU (%p%)\n"s);  // NOLINT(misc-include-cleaner) bug in clang-tidy 18
				aux = ~out;
				out = std::move(aux);
				BOOST_TEST( out[35][79] == in[79][35] );
			}
			BOOST_TEST( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			{
				watch const unnamed{ "transposition with loop   %ws wall, CPU (%p%)\n"s };
				std::for_each(extension(out).begin(), extension(out).end(), [&out](auto idx) {
					auto ext = multi::extension_t(0L, idx);
					std::for_each(ext.begin(), ext.end(), [&out, idx](auto jdx) {
						std::swap(out[idx][jdx], out[jdx][idx]);
					});
				});
				BOOST_TEST( out[35][79] == in[79][35] );
			}
			BOOST_TEST( out == ~in );
		}
		{
			multi::array<complex, 2> out = in;
			{
				watch const unnamed{ "transposition with loop 2 %ws wall, CPU (%p%)\n"s };
				std::for_each(extension(out).begin(), extension(out).end(), [&out](auto idx) {
					auto ext = multi::extension_t(idx + 1, out.size());
					std::for_each(ext.begin(), ext.end(), [&out, idx](auto jdx) {
						std::swap(out[idx][jdx], out[jdx][idx]);
					});
				});
				BOOST_TEST( out[35][79] == in[79][35] );
			}
			BOOST_TEST( out == ~in );
		}
	}
}
return boost::report_errors();}
