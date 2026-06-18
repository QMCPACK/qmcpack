// Copyright 2025-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/fftw.hpp>
#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <chrono>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>  // for is_const_v

namespace multi = boost::multi;

template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;

// NOLINTNEXTLINE(misc-use-internal-linkage)
class watch : private std::chrono::high_resolution_clock {  // NOSONAR(cpp:S4963) this class will report timing on destruction
	std::string label_;
	time_point  start_ = now();

 public:
	explicit watch(std::string_view label) : label_{label} {}

	watch(watch const&) = delete;
	watch(watch&&)      = default;

	auto operator=(watch const&) = delete;
	auto operator=(watch&&)      = delete;

	auto elapsed_sec() const { return std::chrono::duration<double>(now() - start_).count(); }
	~watch() { std::cerr << label_ << ": " << elapsed_sec() << " sec" << '\n'; }
};

namespace {
// clang-format off
template<class Tp>
inline
#ifdef _MSC_VER
__forceinline
#else
__attribute__((always_inline))
#endif
void DoNotOptimize(Tp& value) {  // NOLINT(readability-identifier-naming)
#ifdef _MSC_VER
	_ReadWriteBarrier(); (void)value;
#else
#if defined(__clang__) || defined(__circle_build__)
if constexpr(!std::is_const_v<Tp>) {  // NOLINT(bugprone-branch-clone)
	asm volatile("" : "+r,m"(value) : : "memory");  // NOLINT(hicpp-no-assembler)
} else {  // NOLINT(bugprone-branch-clone)
	asm volatile("" : "+m,r"(value) : : "memory");  // NOLINT(hicpp-no-assembler)
}
#else
	asm volatile("" : "+m,r"(value) : : "memory");  // NOLINT(hicpp-no-assembler)
#endif
#endif
}
// clang-format on
}  // end namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	multi::fftw::environment const env;  // destroyed last, runs fftw_cleanup() after the plan is gone

	using complex = std::complex<double>;

	auto in = [] {
#ifndef NDEBUG
		auto ret = multi::array<complex, 3>({100, 100, 100});
#else
		auto ret = multi::array<complex, 3>({10, 10, 10});
#endif
		std::generate(
			ret.elements().begin(), ret.elements().end(),
			[eng  = std::default_random_engine{std::random_device{}()},
			 dist = std::uniform_real_distribution<>{}]() mutable {
				return complex{dist(eng), dist(eng)};
			}
		);
		return ret;
	}();

	auto const pn = multi::fftw::plan::forward({{false, false, false}}, in.base(), in.layout(), in.base(), in.layout());

	DoNotOptimize(in);

	std::cout << pn.flops() << "FLOPS\n";

	[&, unnamed = watch{"3D *100x100x100"}] {
		for(int i = 0; i != 100; ++i) {  // NOLINT(altera-unroll-loops)
			pn.execute(in.base(), in.base());
		}
	}();

	DoNotOptimize(in);

	return boost::report_errors();
}
