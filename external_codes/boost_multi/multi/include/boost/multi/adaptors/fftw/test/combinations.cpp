// Copyright 2020-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #include <boost/test/tools/fpc_tolerance.hpp>  // TODO(correaa) use lightweight
// #include <boost/test/included/unit_test.hpp>

#include <boost/multi/adaptors/fftw.hpp>

#include <boost/multi/array.hpp>

#include <algorithm>  // for generate, for_each
#include <array>      // for array
#include <chrono>     // for operator-, duration  // NOLINT(build/c++11)
#include <complex>    // for complex, operator==
#include <iostream>   // for operator<<, basic_os...
#include <numeric>    // for iota
#include <random>     // for linear_congruential_...
#include <string>     // for operator""s, operator+
#include <utility>    // for move
#include <vector>     // for vector

namespace multi = boost::multi;

template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;

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

template<class T, multi::dimensionality_type D> using marray = multi::array<T, D>;
constexpr auto exts                                          = multi::extensions_t<4>({6, 12, 24, 12});

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
multi::fftw::environment const env;

BOOST_AUTO_TEST_CASE(fft_combinations) {  // , *boost::unit_test::tolerance(0.00001)) {
	using complex = std::complex<double>;

	auto const in = [] {
		marray<complex, 4> ret(exts);
		std::generate(
			ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng        = std::default_random_engine{std::random_device{}()},
			 uniform_01 = std::uniform_real_distribution<>{}]() mutable {
				return complex{uniform_01(eng), uniform_01(eng)};
			}
		);
		return ret;
	}();

	// NOLINTNEXTLINE(fuchsia-default-arguments-calls)
	std::vector<std::array<bool, 4>> const which_cases = {
		{false,  true,  true,  true},
		{false,  true,  true, false},
		{ true, false, false, false},
		{ true,  true, false, false},
		{false, false,  true, false},
		{false, false, false, false},
	};

	using std::cout;
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	for(auto which : which_cases) {  // NOLINT(altera-unroll-loops)
		cout << "case ";
		std::for_each(which.begin(), which.end(), [](auto elem) { cout << elem << ", "; });
		// copy(which.begin(), which.end(), std::ostream_iterator<bool>{cout, ", "});
		cout << "\n";

		marray<complex, 4> out = in;
		{
			auto const  pln = multi::fftw::plan::forward(which, in.base(), in.layout(), out.base(), out.layout());
			watch const unnamed("cpu_oplac planned %ws wall, CPU (%p%)\n"s);  // NOLINT(misc-include-cleaner) bug in clang-tidy 18
			pln.execute(in.base(), out.base());
		}
		{
			auto in_rw = in;

			watch const unnamed{"cpu_iplac %ws wall, CPU (%p%)\n"s};
			multi::fftw::dft_forward(which, in_rw, in_rw);
		}
		{
			auto in_rw = in;

			auto const pln = multi::fftw::plan::forward(which, in_rw.base(), in_rw.layout(), in_rw.base(), in_rw.layout());

			watch const unnamed{"cpu_iplac planned %ws wall, CPU (%p%)\n"s};
			pln.execute(in_rw.base(), in_rw.base());
		}
		{
			auto        in_rw = in;
			auto const  pln   = multi::fftw::plan::forward(which, in_rw.base(), in_rw.layout(), in_rw.base(), in_rw.layout());
			watch const unnamed{"cpu_iplac planned measured %ws wall, CPU (%p%)\n"s};
			pln.execute(in_rw.base(), in_rw.base());
		}
	}
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark) {  //, *boost::unit_test::enabled()) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	using complex  = std::complex<double>;
	namespace fftw = multi::fftw;

	marray<complex, 4> in(exts);
	std::iota(in.data_elements(), in.data_elements() + in.num_elements(), 1.2);

	BOOST_TEST(in[0][0][0][0] == 1.2);
	std::array<bool, 4> which = {false, true, true, true};
	[&, unnamed = watch{"fftw_4D_power_benchmark inplace FTTT"s}] {
		fftw::dft(which, in, fftw::forward);
	}();
	[&, unnamed = watch{"fftw_4D_power_benchmark inplace FTTT"s}] {
		fftw::dft(which, in, fftw::forward);
	}();
	auto in0000 = in[0][0][0][0];
	BOOST_TEST(in0000 != 1.2);

	marray<complex, 4> out(exts);
	[&, unnamed = watch{"fftw_4D_power_benchmark outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{"fftw_4D_power_benchmark outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{"fftw_4D_power_benchmark outofplace FTTT"s}] {
		fftw::dft(which, in, out, fftw::forward);
	}();
	[&, unnamed = watch{"fftw_4D_power_benchmark outofplace+alloc FTTT"s}] {
		marray<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	[&, unnamed = watch{"fftw_4D_power_benchmark outofplace+alloc FTTT"s}] {
		marray<complex, 4> out2(exts);
		fftw::dft(which, in, out2, fftw::forward);
	}();
	BOOST_TEST(in0000 == in[0][0][0][0]);
}

BOOST_AUTO_TEST_CASE(fftw_4D_power_benchmark_syntax) {
	// NOLINTNEXTLINE(fuchsia-default-arguments-calls) use of std::vector
	std::vector<std::array<bool, 4>> const which_cases = {
		{false,  true,  true,  true},
		{false,  true,  true, false},
		{ true, false, false, false},
		{ true,  true, false, false},
		{false, false,  true, false},
		{false, false, false, false},
	};
	using complex = std::complex<double>;

	auto const in = [] {
		marray<complex, 4> ret(exts);
		std::generate(
			ret.data_elements(), ret.data_elements() + ret.num_elements(),
			[eng        = std::default_random_engine{std::random_device{}()},
			 uniform_01 = std::uniform_real_distribution<>{}]() mutable {
				return complex{uniform_01(eng), uniform_01(eng)};
			}
		);
		return ret;
	}();
}

return boost::report_errors();
}
