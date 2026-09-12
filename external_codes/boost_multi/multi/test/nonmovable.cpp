// Copyright 2022-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for subarray, array, range, operator!=

#include "boost/multi/restriction.hpp"  // for operator^

#include <boost/core/lightweight_test.hpp>

#include <array>
#include <new>  // IWYU pragma: keep
// IWYU pragma: no_include <utility>  // for declval, forward, move

namespace multi = boost::multi;

class ncnm {  // NOLINT(misc-use-internal-linkage)
	int val_;

 public:
	explicit ncnm(int val) : val_{val} {}

	ncnm(ncnm const&) = delete;
	ncnm(ncnm&&)      = delete;

	auto operator=(ncnm const&) = delete;
	auto operator=(ncnm&&)      = delete;

	auto val() const -> int { return val_; }

	~ncnm() = default;
};

class rando {  // NOLINT(misc-use-internal-linkage)
	int val_ = 0;

 public:
	rando() = default;
	explicit operator int() { return val_++; }
};

template<> inline constexpr bool boost::multi::force_element_trivial_default_construction<ncnm> = true;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		multi::array<ncnm, 2> arr2d({2, 2}, boost::multi::uninitialized_elements);
		new (&arr2d[1][1]) ncnm(99);
	}
	{
		std::array<ncnm, 4> arr1d{
			{ncnm(0), ncnm(1), ncnm(2), ncnm(3)}
		};
		multi::array_ref<ncnm, 2> arr2d_ref({2, 2}, arr1d.data());

		BOOST_TEST( arr2d_ref.size() == 2 );
		BOOST_TEST( arr2d_ref[1][1].val() == 3 );
	}
	{
		multi::array<int, 2> arr2d_res1 = [](auto ii, auto jj) { return static_cast<int>((ii * 2) + jj); } ^ multi::extents_t{2, 2};
		// multi::array<int, 2> arr2d_res2 = [rr = rando{}](auto, auto) mutable { return static_cast<int>(rr); } ^ multi::extents_t<2>({2, 2});

		// std::cout << arr2d_res1 << std::endl;
		// std::cout << arr2d_res2 << std::endl;

		BOOST_TEST( arr2d_res1[1][1] == 3 );
	}
	{
		multi::array arr2d_res1 = [](auto ii, auto jj) { return (ii * 2) + jj; } ^ multi::extents_t{2, 2};

		BOOST_TEST( arr2d_res1[1][1] == 3 );

#ifdef __cpp_multitdimensional_subscrip
		BOOST_TEST( arr2d_res1[1, 1] == 3 );
#endif
	}

	return boost::report_errors();
}
