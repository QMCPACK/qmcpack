// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, apply, operator==
#include <boost/multi/io.hpp>

// #include <boost/multi/detail/what.hpp>

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <iostream>

namespace ion {
template<class T>
struct segmented_iterator_traits;
}  // end namespace ion

// based on https://boostedcpp.net/2026/05/18/neoclassical-c-segmented-iterators-revisited-1/
namespace ion {
template<typename E, ::boost::multi::dimensionality_type D, typename EP, bool IsConst, bool IsMove, typename Stride>
struct segmented_iterator_traits<::boost::multi::detail::array_iterator<E, D, EP, IsConst, IsMove, Stride>> {
	using segment_iterator        = ::boost::multi::detail::array_iterator<E, 2, EP, IsConst, IsMove, typename Stride::stride2_type>;
	using segmented_iterator_type = ::boost::multi::detail::array_iterator<E, D, EP, IsConst, IsMove, Stride>;

	static auto segment(segmented_iterator_type self) { return self.outer(); }
	static auto local(segmented_iterator_type self) { return self.local(); }
};

namespace {
template<class SegIt, class T>
void fill_segmented(SegIt first, SegIt last, T x) {
	using traits = ion::segmented_iterator_traits<SegIt>;

	typename traits::segment_iterator       sf = traits::segment(first);
	typename traits::segment_iterator const sl = traits::segment(last);
	// typename traits::local_iterator   lf = traits::local(first);

	// auto sf = first.outer();
	// auto sl = last.outer();

	auto lf = first.local();

	while(true) {  // NOLINT(altera-unroll-loops)
		// typename traits::local_iterator le =
		//     (sf == sl) ? traits::local(last) : traits::end(sf);

		auto le = (sf == sl) ? last.local() : (*sf).end();

		std::fill(lf, le, x);
		if(sf == sl) {
			break;
		}
		// lf = traits::begin(++sf);
		lf = (*(++sf)).begin();
	}
}
}  // end namespace

}  // end namespace ion

namespace {

template<class SegIt, class T>
void fill_segmented(SegIt first, SegIt last, T x) {
	// typedef segmented_iterator_traits<SegIt> traits;
	// typename traits::segment_iterator sf = traits::segment(first);
	// typename traits::segment_iterator sl = traits::segment(last);
	// typename traits::local_iterator   lf = traits::local(first);

	auto sf = first.outer();
	auto sl = last.outer();

	auto lf = first.local();

	while(true) {  // NOLINT(altera-unroll-loops)
		// typename traits::local_iterator le =
		//     (sf == sl) ? traits::local(last) : traits::end(sf);

		auto le = (sf == sl) ? last.local() : (*sf).end();

		std::fill(lf, le, x);
		if(sf == sl) {
			break;
		}
		// lf = traits::begin(++sf);
		lf = (*(++sf)).begin();
	}
}

}  // namespace

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity)
	multi::array<int, 2> arr1({2, 3});
	std::fill(arr1().elements().begin(), arr1().elements().end(), 7);

	BOOST_TEST( arr1[1][1] == 7 );

	multi::array<int, 2> arr2 = {
		{1, 2, 3},
		{4, 5, 6}
	};
	BOOST_TEST( arr2.flattened().begin().segment().size() == 3 );

	BOOST_TEST( arr2.flattened().begin().segment()[0] == 1 );
	BOOST_TEST( arr2.flattened().begin().segment()[1] == 2 );
	BOOST_TEST( arr2.flattened().begin().segment()[2] == 3 );

	BOOST_TEST( arr2.flattened().begin().segment().num_elements() == 3 );

	// span to one-past-end of the flattened (segmented) layout: (nsegs-1) outer strides + last segment
	BOOST_TEST( arr2.flattened().layout().nelems() == 6 );  // contiguous 2x3: 6 - 3 + 3, NOT 6 + 3 + 3

	{
		multi::array<int, 2> big({2, 5}, 0);
		auto                 strided = big({0, 2}, {0, 3});  // 2x3 view with outer stride 5 (gap between rows)
		// must skip the inter-row gap past the last element: 5 (one outer stride) + 3 (last segment) == 8, NOT 5 + 5 + 3
		BOOST_TEST( strided.flattened().layout().nelems() == 8 );
	}

	BOOST_TEST( arr2[0].size() == 3 );

	BOOST_TEST( arr2[0][0] == 1 );
	BOOST_TEST( arr2[0][1] == 2 );
	BOOST_TEST( arr2[0][2] == 3 );

	BOOST_TEST( arr2[0] == arr2.flattened().begin().segment() );

	BOOST_TEST( (arr2.flattened().begin() + 0).segment() == arr2[0] );
	BOOST_TEST( (arr2.flattened().begin() + 1).segment() == arr2[0] );
	BOOST_TEST( (arr2.flattened().begin() + 2).segment() == arr2[0] );

	BOOST_TEST( (arr2.flattened().begin() + 3).segment() == arr2[1] );
	BOOST_TEST( (arr2.flattened().begin() + 4).segment() == arr2[1] );
	BOOST_TEST( (arr2.flattened().begin() + 5).segment() == arr2[1] );

	BOOST_TEST( (arr2.flattened().begin() + 0).local() == arr2[0].begin() );
	BOOST_TEST( (arr2.flattened().begin() + 1).local() == arr2[0].begin() + 1 );
	BOOST_TEST( (arr2.flattened().begin() + 2).local() == arr2[0].begin() + 2 );

	BOOST_TEST( (arr2.flattened().begin() + 3).local() == arr2[1].begin() );
	BOOST_TEST( (arr2.flattened().begin() + 4).local() == arr2[1].begin() + 1 );
	BOOST_TEST( (arr2.flattened().begin() + 5).local() == arr2[1].begin() + 2 );

	std::cout << "arr2 = " << arr2 << '\n';
	std::cout << "arr2.flattened().begin().segment() = " << arr2.flattened().begin().segment() << '\n';
	std::cout << "arr2.flattened().(begin() + 3).segment() = " << (arr2.flattened().begin() + 3).segment() << '\n';
	std::cout << "*(arr2.flattened().begin().global()) = " << *(arr2.flattened().begin().outer()) << '\n';
	std::cout << "*((arr2.flattened().begin() + 3).global()) = " << *((arr2.flattened().begin() + 3).outer()) << '\n';

	ion::fill_segmented(arr2.flattened().begin(), arr2.flattened().end(), 7);

	fill_segmented(arr2.flattened().begin(), arr2.flattened().end(), 7);

	// auto gf = arr2.flattened().begin().global();
	// auto gl = arr2.flattened().end().global();

	// // auto lf = arr2.flattened().begin().local();

	// std::cout << "*gf = " << *gf << std::endl;
	// std::cout << "*(gl - 1) = " << *(gl - 1) << std::endl;

	// BOOST_TEST( gf != gl );
	std::cout << arr2 << '\n';
	std::cout << arr1 << '\n';

	BOOST_TEST( arr1 == arr2 );

	return boost::report_errors();
}
