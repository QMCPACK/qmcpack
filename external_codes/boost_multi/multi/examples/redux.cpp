// Copyright 2018-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wunused-macros"
#endif
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#pragma clang diagnostic ignored "-Wpadded"
#endif
#if defined(_MSC_VER)
#pragma warning(disable : 4244)  // warning C4244: 'initializing': conversion from '_Ty' to '_Ty', possible loss of data
#endif

#include <boost/multi/adaptors/blas.hpp>  // IWYU pragma: keep
#include <boost/multi/array.hpp>          // for array, implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

#include <algorithm>   // IWYU pragma: keep
#include <chrono>      // NOLINT(build/c++11)
#include <cmath>       // IWYU pragma: keep
#include <functional>  // IWYU pragma: keep
#include <iostream>
#include <numeric>  // IWYU pragma: keep
#include <random>   // IWYU pragma: keep
#include <string>
#include <string_view>
#include <type_traits>  // for std::decay_t
#include <utility>      // for move  // IWYU pragma: keep  // NOLINT(misc-include-cleaner) bug in clang-tidy 19

// IWYU pragma: no_include <pstl/glue_numeric_impl.h>         // for reduce, transform_reduce
// IWYU pragma: no_include <cstdlib>                          // for abs
// IWYU pragma: no_include <new>                              // for bad_alloc

#ifndef __NVCC__
#if defined(__has_include) && __has_include(<execution>) && (!defined(__INTEL_LLVM_COMPILER) || (__INTEL_LLVM_COMPILER > 20240000))
#if !(defined(__clang__) && defined(__CUDA__))
#include <execution>  // IWYU pragma: keep
#define HAS_STD_EXECUTION 1
#endif
#endif
#endif

namespace multi = boost::multi;

#if (__cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L))
#if __has_include(<generator>)
#include <generator>
#endif
#endif

#if defined(__cpp_lib_generator) && (__cpp_lib_generator >= 202207L) && !defined(_MSC_VER)
template<class Arr2D>
auto co_extensions_elements(Arr2D const& arr2d) -> std::generator<typename Arr2D::indexes> {
	auto const [is, js] = arr2d.extensions();
	for(auto const i : is) {
		for(auto const j : js) {
			co_yield typename Arr2D::indexes{i, j};
		}
	}
}

template<class Arr2D>
std::generator<typename Arr2D::element_cref>
co_celements(Arr2D const& arr2d) {
	auto const [is, js] = arr2d.extensions();
	for(auto const i : is) {
		for(auto const j : js) {
			co_yield arr2d[i][j];
		}
	}
}

#endif

class watch {
	std::chrono::time_point<std::chrono::high_resolution_clock> start_ = std::chrono::high_resolution_clock::now();

	std::string msg_;
	bool        running_ = true;

	template<class T>
#if defined(_MSC_VER)
	inline __forceinline static auto do_not_optimize_(T&& value) noexcept -> T&& {
		return std::forward<T>(value);
	}
#else
	inline __attribute__((always_inline))  // NOLINT(readability-redundant-inline-specifier)
	static auto
	do_not_optimize_(T&& value) noexcept -> T&& {
		if constexpr(std::is_pointer_v<T>) {
			asm volatile("" : "+m"(value)::"memory");  // NOLINT(hicpp-no-assembler)
		} else {
			asm volatile("" : "+r"(value)::);  // NOLINT(hicpp-no-assembler)
		}
		return std::forward<T>(value);
	}
#endif

 public:
	explicit watch(std::string_view msg) : msg_(msg) {}  // NOLINT(fuchsia-default-arguments-calls)
	template<class T>
	auto lap(T&& some) const -> T&& {
		do_not_optimize_(const_cast<std::decay_t<T>*>(&some));  // NOLINT(cppcoreguidelines-pro-type-const-cast)
		std::cerr << msg_ << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start_).count() << " ms\n";
		return std::forward<T>(some);
	}
	template<class T>
	void stop(T&& some) {
		if(running_) {
			running_ = false;
			lap(std::forward<T>(some));
		}
	}

	~watch() {
		if(running_) {
			stop(*this);
		}
	}
	watch(watch const&)          = delete;
	watch(watch&&)               = delete;
	auto operator=(watch const&) = delete;
	auto operator=(watch&&)      = delete;
	//  non-default destructor but does not define a copy constructor, a copy assignment operator, a move constructor or a move assignment operator
};

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// multi::array<double, 2>::size_type const maxsize = 39062;  // 390625;
	// multi::array<double, 2>::size_type const nmax    = 1000;   // 10000;

	// auto pp = [] /*__host__ __device__*/ (long ix, long iy) -> double { return double(ix) * double(iy); };

	auto const nx = 40000;  // nmax;     // for(long nx = 1; nx <= nmax; nx *= 10)
	auto const ny = 2000;   // maxsize;  // for(long ny = 1; ny <= maxsize; ny *= 5)

	// auto const nx = 4000;  // nmax;     // for(long nx = 1; nx <= nmax; nx *= 10)
	// auto const ny = 200;   // maxsize;  // for(long ny = 1; ny <= maxsize; ny *= 5)

	// auto total = nx*ny;

	// nx = 2;
	// ny = total / nx;

	multi::array<double, 2> K2D({nx, ny});

	for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {      // NOLINT(altera-id-dependent-backward-branch)
		for(multi::array<double, 2>::index iy = 0; iy != ny; ++iy) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops)
			K2D[ix][iy] = static_cast<double>(ix) * static_cast<double>(iy);
		}
	}

	{
		auto const accumulator = [&]() {
			multi::array<double, 1> ret({nx}, 0.0);
			for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {      // NOLINT(altera-id-dependent-backward-branch)
				for(multi::array<double, 2>::index iy = 0; iy != ny; ++iy) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops)
					ret[ix] += K2D[ix][iy];
				}
			}
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - (static_cast<double>(ix) * ny * (ny - 1.0) / 2.0) ) < 1.0e-8);
		}
	}

#if defined(NDEBUG) && !defined(RUNNING_ON_VALGRIND)

	{
		auto const accumulator = [&](watch = watch("raw loop")) {  // NOLINT(fuchsia-default-arguments-declarations)
			multi::array<double, 1> ret({nx}, 0.0);
			for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {      // NOLINT(altera-id-dependent-backward-branch)
				for(multi::array<double, 2>::index iy = 0; iy != ny; ++iy) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops)
					ret[ix] += K2D[ix][iy];
				}
			}
			return ret;
		}();  // NOLINT(fuchsia-default-arguments-calls)

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&] {
			watch const _("accumulate for");
			return std::accumulate(
				(~K2D).begin(), (~K2D).end(),
				multi::array<double, 1>(multi::array<double, 1>::extensions_type{K2D.extension()}, 0.0),
				[](auto const& acc, auto const& col) {
					multi::array<double, 1> res(acc.extensions());
					for(auto const i : col.extension()) {  // NOLINT(altera-unroll-loops)
						res[i] = acc[i] + col[i];
					}
					return res;
				}
			);
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&](auto init) {
			watch const _("accumulate move");
			return std::accumulate(
				(~K2D).begin(), (~K2D).end(), std::move(init), [](auto&& acc, auto const& col) {
					multi::array<double, 1> ret(std::forward<decltype(acc)>(acc));
					for(auto const i : col.extension()) {  // NOLINT(altera-unroll-loops)
						ret[i] += col[i];
					}
					return ret;
				}
			);
		}(multi::array<double, 1>(K2D.extension(), 0.0));

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&](auto init) {
			watch const _("accumulate forward");
			return std::accumulate(
				(~K2D).begin(), (~K2D).end(), std::move(init), [](auto&& acc, auto const& col) -> decltype(acc) {
					for(auto const i : col.extension()) {  // NOLINT(altera-unroll-loops)
						acc[i] += col[i];
					}
					return std::forward<decltype(acc)>(acc);
				}
			);
		}(multi::array<double, 1>(K2D.extension(), 0.0));

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&](auto init, watch = watch("accumulate transform forward")) {  // NOLINT(fuchsia-default-arguments-declarations)
			return std::accumulate(
				(~K2D).begin(), (~K2D).end(), std::move(init), [](auto&& acc, auto const& col) -> decltype(acc) {
					std::transform(col.begin(), col.end(), acc.begin(), acc.begin(), [](auto const& cole, auto&& acce) { return std::forward<decltype(acce)>(acce) + cole; });
					return std::forward<decltype(acc)>(acc);
				}
			);
		}(multi::array<double, 1>(K2D.extension(), 0.0));  // NOLINT(fuchsia-default-arguments-calls)

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

#if (!defined(__GLIBCXX__) || (__GLIBCXX__ >= 20190502))
	{
		auto const accumulator = [&] {
			watch const _("reduce transform forward");
			return std::reduce(
				(~K2D).begin(), (~K2D).end(), multi::array<double, 1>(K2D.extension(), 0.0), [](auto acc, auto const& col) {
					multi::array<double, 1> ret(std::move(acc));
					std::transform(col.begin(), col.end(), ret.begin(), ret.begin(), [](auto const& cole, auto&& acce) { return std::forward<decltype(acce)>(acce) + cole; });
					return ret;
				}
			);
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}
#endif

	{
		auto const accumulator = [&] {
			watch const _("transform accumulate element zero");

			multi::array<double, 1> ret(K2D.extension());
			std::transform(
				K2D.begin(), K2D.end(), ret.begin(), [](auto const& row) { return std::accumulate(row.begin(), row.end(), 0.0); }
			);
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

#if (!defined(__GLIBCXX__) || (__GLIBCXX__ >= 20190502))
	{
		auto const accumulator = [&] {
			watch const _("transform reduce element zero");

			multi::array<double, 1> ret(K2D.extension());
			std::transform(
				K2D.begin(), K2D.end(), ret.begin(), [](auto const& row) { return std::reduce(row.begin(), row.end(), 0.0); }
			);
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}
#endif

	{
		auto const accumulator = [&](auto&& init) {
			watch const _("transform accumulate");
			std::transform(
				K2D.begin(), K2D.end(), init.begin(), init.begin(), [](auto const& row, auto rete) { return std::accumulate(row.begin(), row.end(), std::move(rete)); }
			);
			return std::forward<decltype(init)>(init);
		}(multi::array<double, 1>(K2D.extension(), 0.0));

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

#if (!defined(__GLIBCXX__) || (__GLIBCXX__ >= 20200000))
	{
		auto const accumulator = [&](auto&& init) {
			watch const _("> transform reduce");
			std::transform(
				K2D.begin(), K2D.end(), init.begin(), init.begin(), [](auto const& row, auto rete) { return std::reduce(row.begin(), row.end(), std::move(rete)); }
			);
			return std::forward<decltype(init)>(init);
		}(multi::array<double, 1>(K2D.extension(), 0.0));

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}
#endif

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
#if (defined(__has_include) && __has_include(<execution>))
#if !defined(__NVCC__) && !defined(__NVCOMPILER) && !(defined(__clang__) && defined(__CUDA__))
#if (!defined(__clang_major__) || (__clang_major__ > 7))
#if (!defined(__GLIBCXX__) || (__GLIBCXX__ >= 20220000)) && !defined(_LIBCPP_VERSION)
#if !defined(__apple_build_version__) && (!defined(__INTEL_LLVM_COMPILER) || (__INTEL_LLVM_COMPILER > 20240000))
	{
		auto const accumulator = [&](watch = watch("transform reduce[unseq]")) {  // NOLINT(fuchsia-default-arguments-declarations)
			multi::array<double, 1> ret(K2D.extension(), 0.0);
			std::transform(
				K2D.begin(), K2D.end(),
				ret.begin(),
				ret.begin(),
				[](auto const& row, auto rete) { return std::reduce(std::execution::unseq, row.begin(), row.end(), std::move(rete)); }
			);
			return ret;
		}();  // NOLINT(fuchsia-default-arguments-calls)

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&] {
			watch const _("transform reduce[par]");

			multi::array<double, 1> ret(K2D.extension(), 0.0);
			std::transform(
				K2D.begin(), K2D.end(),
				ret.begin(),
				ret.begin(),
				[](auto const& row, auto rete) { return std::reduce(std::execution::par, row.begin(), row.end(), std::move(rete)); }
			);
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&] {
			watch const _("transform reduce[par_unseq]");

			multi::array<double, 1> ret(K2D.extension(), 0.0);
			std::transform(
				K2D.begin(), K2D.end(),
				ret.begin(),
				ret.begin(),
				[](auto const& row, auto rete) { return std::reduce(std::execution::par_unseq, row.begin(), row.end(), std::move(rete)); }
			);
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&]() {
			watch const _("transform[par] reduce");

			multi::array<double, 1> ret(K2D.extension(), 0.0);
			std::transform(
				std::execution::par,
				K2D.begin(), K2D.end(),
				ret.begin(),
				ret.begin(),
				[](auto const& row, auto rete) { return std::reduce(row.begin(), row.end(), std::move(rete)); }
			);
			return ret;
		}();

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&](auto ret) {
			watch const _("* transform[par] reduce[unseq]");
			std::transform(
				std::execution::par,
				K2D.begin(), K2D.end(),
				ret.begin(),
				ret.begin(),
				[](auto const& row, auto rete) { return std::reduce(std::execution::unseq, row.begin(), row.end(), std::move(rete)); }
			);
			return ret;
		}(multi::array<double, 1>(K2D.extension(), 0.0));

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		multi::array<double, 1> accumulator(K2D.extension(), 0.0);
		[&](auto acc_begin) {
			watch const _("transform[par] reduce[unseq] iterator");
			return std::transform(
				std::execution::par,
				K2D.begin(), K2D.end(),
				acc_begin, acc_begin,
				[](auto const& row, auto rete) { return std::reduce(std::execution::unseq, row.begin(), row.end(), std::move(rete)); }
			);
		}(accumulator.begin());

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
		}
	}

	{
		auto const accumulator = [&](auto zero_elem, watch = watch("transform[par] reduce[unseq] element zero")) {  // NOLINT(fuchsia-default-arguments-declarations)
			multi::array<double, 1> ret(K2D.extension());
			std::transform(
				std::execution::par,
				K2D.begin(), K2D.end(),
				ret.begin(),
				[zz = std::move(zero_elem)](auto const& row) { return std::reduce(std::execution::unseq, row.begin(), row.end(), std::move(zz)); }
			);
			return ret;
		}(0.0);  // NOLINT(fuchsia-default-arguments-calls)

		for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-5);
		}
	}
#endif
#endif
#endif
#endif
#endif
#endif  // __NVCC__
#endif

	// {
	//  auto const accumulator = [&](auto&& init) {
	//      watch const             _("blas gemv");
	//      multi::array<double, 1> ones({init.extension()}, 1.0);
	//      multi::blas::gemv_n(1.0, K2D.begin(), K2D.size(), ones.begin(), 0.0, init.begin());
	//      return std::forward<decltype(init)>(init);
	//  }(multi::array<double, 1>(K2D.extension(), 0.0));

	//  for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
	//      BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
	//  }
	// }

	// {
	//  auto const accumulator = [&](auto&& init) {
	//      watch const _("blas gemv smart");
	//      multi::blas::gemv_n(1.0, K2D.begin(), K2D.size(), init[0].begin(), 0.0, init[1].begin());
	//      return +init[1];
	//  }(multi::array<double, 2>({2, K2D.extension()}, 1.0));

	//  for(multi::array<double, 2>::index ix = 0; ix != nx; ++ix) {  // NOLINT(altera-unroll-loops)
	//      BOOST_TEST( std::abs( accumulator[ix] - static_cast<double>(ix) * ny * (ny - 1.0) / 2.0 ) < 1.0e-8);
	//  }
	// }

#if defined(NDEBUG)
	// Chris
	{
		multi::size_t const em  = 1200 / 2;
		multi::size_t const en  = 1000 / 2;
		multi::size_t const ell = 800 / 2;

		multi::array<double, 3> a3d({em, en, ell});
		multi::array<double, 2> b2d({en, ell});

		std::mt19937                     gen(42);  // NOLINT(cert-msc32-c,cert-msc51-cpp) use for unpredictable std::random_device{}
		std::uniform_real_distribution<> distrib;

		std::generate(a3d.elements().begin(), a3d.elements().end(), [&]() { return distrib(gen); });
		std::generate(b2d.elements().begin(), b2d.elements().end(), [&]() { return distrib(gen); });

		// std::iota(a3d.elements().begin(), a3d.elements().end(), 20.0);
		// std::iota(b2d.elements().begin(), b2d.elements().end(), 30.0);

		multi::array<double, 1> c_gold(em, 0.0);
		{
			for(multi::index k = 0; k != em; ++k) {
				for(multi::index j = 0; j != en; ++j) {       // NOLINT(altera-unroll-loops)
					for(multi::index i = 0; i != ell; ++i) {  // NOLINT(altera-unroll-loops)
						c_gold[k] += a3d[k][j][i] * b2d[j][i];
					}
				}
			}
			c_gold = multi::array<double, 1>(em, 0.0);

			for(multi::index const k : a3d.extension()) {
				auto const& a3dk = a3d[k];
				for(multi::index const j : b2d.extension()) {  // NOLINT(altera-unroll-loops)
					auto const& a3dkj = a3dk[j];
					auto const& b2dj  = b2d[j];
					for(multi::index const i : b2dj.extension()) {  // NOLINT(altera-unroll-loops)
						c_gold[k] += a3dkj[i] * b2dj[i];
					}
				}
			}
		}

		{
			watch _("chris raw 3-loop");

			multi::array<double, 1> c_flat(em, 0.0);

			for(multi::index const k : a3d.extension()) {
				auto const& a3dk = a3d[k];
				for(multi::index const j : b2d.extension()) {  // NOLINT(altera-unroll-loops)
					auto const& a3dkj = a3dk[j];
					auto const& b2dj  = b2d[j];
					for(multi::index const i : b2dj.extension()) {  // NOLINT(altera-unroll-loops)
						c_flat[k] += a3dkj[i] * b2dj[i];
					}
				}
			}
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris raw 2-loop flat");

			multi::array<double, 1> c_flat(em, 0.0);

			for(auto const k : c_flat.extension()) {
				auto const& a3d_rowes    = a3d[k].elements().base();
				auto const& b2d_elements = b2d.elements().base();
				for(auto const ji : b2d.elements().extension()) {   // NOLINT(altera-unroll-loops)
					c_flat[k] += a3d_rowes[ji] * b2d_elements[ji];  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				}
			}
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris raw 2-loop flat reversed");

			multi::array<double, 1> c_flat(em, 0.0);

			for(auto const ji : b2d.elements().extension()) {  // NOLINT(altera-unroll-loops)
				auto const& b2deji = b2d.elements().base()[ji];
				for(auto const k : c_flat.extension()) {  // NOLINT(altera-unroll-loops)
					c_flat[k] += a3d[k].elements().base()[ji] * b2deji;
				}
			}
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris raw 1-loop flat reversed transform");

			multi::array<double, 1> c_flat(em, 0.0);

			for(auto const ji : b2d.elements().extension()) {  // NOLINT(altera-unroll-loops)
				std::transform(c_flat.begin(), c_flat.end(), a3d.begin(), c_flat.begin(), [&](auto&& c_flat_elem, auto const& a3d_row) {
					return std::forward<decltype(c_flat_elem)>(c_flat_elem) + a3d_row.elements().base()[ji] * b2d.elements().base()[ji];
				});
			}
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris raw 1-loop transform_reduce");

			multi::array<double, 1> c_flat(em);

			for(auto const k : c_flat.extension()) {  // NOLINT(altera-unroll-loops)
				auto const& a3dkes = a3d[k].elements();
				c_flat[k]          = std::transform_reduce(
                    a3dkes.base(), a3dkes.base() + a3dkes.size(),
                    b2d.elements().base(), 0.0
                );
			}
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch                   _("chris transform transform_reduce pointer");
			multi::array<double, 1> c_flat(em);

			std::transform(
				a3d.begin(), a3d.end(), c_flat.begin(),
				[&](auto const& a3d_row) {
					return std::transform_reduce(
						a3d_row.base(), a3d_row.base() + a3d_row.elements().size(),
						b2d.base(), 0.0
					);
				}
			);
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}
#if defined(__cpp_lib_generator) && (__cpp_lib_generator >= 202207L) && !defined(_MSC_VER)
		{
			sasa
				watch               _("chris transform transform_reduce elements");
			multi::array<double, 1> c_flat(em);

			std::transform(
				a3d.begin(), a3d.end(), c_flat.begin(),
				[&](auto const& a3d_row) {
					return std::transform_reduce(
						co_celements(a3d_row).begin(), co_celements(a3d_row).end(),
						co_celements(b2d).begin(), 0.0
					);
				}
			);
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}
#endif

#if defined(HAS_STD_EXECUTION)
#if defined(__cpp_lib_execution)
		{
			watch                   _("chris transform(par) transform_reduce");
			multi::array<double, 1> c_flat(em);

			std::transform(
				std::execution::par,
				a3d.begin(), a3d.end(), c_flat.begin(),
				[&](auto const& a3d_row) {
					return std::transform_reduce(
						std::execution::par,
						a3d_row.base(), a3d_row.base() + a3d_row.elements().size(),
						b2d.base(), 0.0
					);
				}
			);
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}
#endif
#endif
		{
			watch _("chris accumulate");

			auto const c_flat = std::accumulate(
				b2d.elements().extension().begin(), b2d.elements().extension().end(),
				multi::array<double, 1>(em, 0.0),
				[&](auto&& acc, auto const ij) {
					for(auto const k : acc.extension()) {  // NOLINT(altera-unroll-loops)
						acc[k] += a3d[k].base()[ij] * b2d.base()[ij];
					}
					return std::forward<decltype(acc)>(acc);
				}
			);
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris accumulate move");

			auto const c_flat = std::accumulate(
				b2d.elements().extension().begin(), b2d.elements().extension().end(),
				multi::array<double, 1>(em, 0.0),
				[&](auto&& acc, auto const ij) {
					std::transform(
						acc.begin(), acc.end(), a3d.begin(), acc.begin(),
						[&](auto&& acce, auto const& a3de) { return std::forward<decltype(acce)>(acce) + a3de.base()[ij] * b2d.base()[ij]; }
					);
					return std::forward<decltype(acc)>(acc);
				}
			);
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}

		{
			watch _("chris transform reduce move transforms");

			auto const c_flat = [&] {
				return std::transform_reduce(
					b2d.elements().extension().begin(), b2d.elements().extension().end(),
					multi::array<double, 1>(em, 0.0),
					[&](auto&& acc, auto const& rhs) {
						std::transform(
							acc.begin(), acc.end(), rhs.begin(), acc.begin(),
							[&](auto&& acce, auto const& rhse) { return std::forward<decltype(acce)>(acce) + rhse; }
						);
						return std::forward<decltype(acc)>(acc);
					},
					[&](auto const ij) {
						multi::array<double, 1> ret(em);
						std::transform(
							a3d.begin(), a3d.end(), ret.begin(),
							[&](auto const& a3de) { return a3de.base()[ij] * b2d.base()[ij]; }
						);
						return ret;
					}
				);
			}();
			_.stop(c_flat);

			BOOST_TEST( std::transform_reduce(c_gold.begin(), c_gold.end(), c_flat.begin(), 0.0, std::plus<>{}, [](auto const& alpha, auto const& omega) { return std::abs(alpha - omega); }) < 1.0e-5 );
		}
	}
#endif
	return boost::report_errors();
}
