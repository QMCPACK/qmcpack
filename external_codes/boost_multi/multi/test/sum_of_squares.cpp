// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <functional>   // IWYU pragma: keep   // for std::plus
#include <numeric>      // for std::transform_reduce, std::accumulate (fallback for libstdc++ < 9)
#include <type_traits>  // for std::is_same_v, std::decay_t
#include <utility>      // for std::forward

#if __has_include(<version>)
#include <version>  // IWYU pragma: keep   // for the feature-test macros used to guard the paths below
#endif

#ifdef __cpp_lib_ranges_fold
#include <algorithm>  // for std::ranges::fold_left  (C++23)
#include <ranges>     // for std::views::transform
#endif

// libstdc++'s Parallel STL (<execution>, in pstl/glue_numeric_impl.h) calls
//   __is_vectorization_preferred<_ExecutionPolicy, _It>(__exec)
// with EXPLICIT template args, turning `_ExecutionPolicy&&` from a forwarding
// reference into a fixed rvalue-ref, so the lvalue `__exec` no longer binds. This
// is a bug in the *library* headers, fixed in libstdc++ 14, so the discriminator is
// the libstdc++ version (_GLIBCXX_RELEASE), not the compiler: libstdc++ <= 13 breaks
// for every front end (gcc 13, clang 16..20, icpx, clang-cuda; confirmed in CI).
// nvcc is excluded too (cudafe++ mis-deduces the same helpers).  // && !defined(_MSC_VER)
#if defined(__cpp_lib_parallel_algorithm) && !defined(__NVCC__) && /* NOLINTNEXTLINE(misc-include-cleaner) */ \
	!(defined(__GLIBCXX__) && (_GLIBCXX_RELEASE < 14))             /* libstdc++ <= 13: broken pstl call site; fixed in libstdc++ 14 */
#define MULTI_HAS_PARALLEL_EXECUTION 1
#include <execution>  // for std::execution::par / parallel_policy
#endif

namespace multi = boost::multi;

namespace {

auto sos(int N) {  // NOLINT(readability-identifier-length)  // N is the number of integers to sum
	using multi::range;

#ifdef __cpp_lib_ranges_fold  // C++23 ranges form (sequential)
	return std::ranges::fold_left(
		range(0, N) | std::views::transform([](auto const& e) noexcept { return e * e; }),
		0,
		std::plus<>{}
	);
#elif !defined(__GLIBCXX__) || _GLIBCXX_RELEASE >= 9  // libstdc++ < 9 (gcc 7/8) has no std::transform_reduce
	return std::transform_reduce(
		range(0, N).begin(), range(0, N).end(),
		0,
		std::plus<>{},
		[](auto const& e) noexcept { return e * e; }
	);
#else                                                 // manual fallback: gcc 7/8
	return std::accumulate(
		range(0, N).begin(), range(0, N).end(),
		0,
		[](auto acc, auto const& e) noexcept { return acc + e * e; }
	);
#endif
}

struct no_policy_t {};  // `std::execution::parallel_policy` can't be copy-list-initialized from `{}` on MSVC

template<class ExecutionPolicy = no_policy_t>  // default policy is "none" -> the sequential branch below
auto sos(ExecutionPolicy&& ep, int N) {        // NOLINT(readability-identifier-length)  // N is the number of integers to sum
	using multi::range;

	if constexpr(std::is_same_v<ExecutionPolicy, no_policy_t>) {  // no <execution>: drop the policy, there is no (policy, ...) overload
		(void)std::forward<ExecutionPolicy>(ep);
		return sos(N);
	}
#ifdef MULTI_HAS_PARALLEL_EXECUTION  // avoid non-dependent std::transform_reduce lookup failing on gcc 7/8, where the branch is never instantiated but still parsed
	else {
		return std::transform_reduce(
			std::forward<ExecutionPolicy>(ep),
			range(0, N).begin(), range(0, N).end(), 0,
			std::plus<>{},
			[](auto const& e) { return e * e; }
		);
	}
#endif
}

}  // namespace

auto main() -> int {
	BOOST_TEST( sos(4) == 0 + 1 + 4 + 9 );
#ifdef MULTI_HAS_PARALLEL_EXECUTION
	BOOST_TEST( sos(std::execution::par, 4) == sos(4) );
#endif
	BOOST_TEST( sos({}, 4) == sos(4) );  // default policy via default template argument

	return boost::report_errors();
}
