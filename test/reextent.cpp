// Copyright 2018-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, static_array, num_elements

#include <boost/core/lightweight_test.hpp>

#include <algorithm>
#include <initializer_list>  // for initializer_list
// #include <iterator>          // for size
#include <tuple>        // IWYU pragma: keep  // for get
#include <type_traits>  // for make_unsigned_t
#include <utility>      // for move
#include <vector>       // for vector

namespace multi = boost::multi;

namespace {
template<class T, class U>
constexpr auto comp_equal(T left, U right) noexcept -> bool {
	using UT = std::make_unsigned_t<T>;
	using UU = std::make_unsigned_t<U>;
	if constexpr(std::is_signed_v<T> == std::is_signed_v<U>) {
		return left == right;
	} else if constexpr(std::is_signed_v<T>) {
		return left < 0 ? false : static_cast<UT>(left) == right;
	} else {
		return right < 0 ? false : left == UU(right);
	}
#if !defined(__INTEL_COMPILER) && !defined(__NVCOMPILER) && !defined(_MSC_VER)
	__builtin_unreachable();
#endif
}
}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(array_reextent)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		multi::array<double, 2> arr3({2, 3});
		BOOST_TEST(size(arr3) == 2);
		BOOST_TEST(size(arr3[0]) == 3);

		arr.reextent({5, 4}, 990);
		BOOST_TEST( arr.num_elements()== 5L*4L );
		BOOST_TEST( arr[1][2] ==  60 );   // reextent preserves values when it can...
		BOOST_TEST( arr[4][3] == 990 );  // ...and gives selected value to the rest
	}
	{
		multi::array<int, 2> arr({2, 3});

		BOOST_TEST( arr.size() == 2 );

		arr = multi::array<int, 2>{};

		BOOST_TEST( arr.size() == 0 );
	}

#ifndef __circle_build__
	{
		multi::array<int, 2> const arr({2, 3});

		BOOST_TEST( arr.size() == 2 );

		multi::array<int, 2> const brr({});

		BOOST_TEST( brr.size() == 0 );
	}

	{
		multi::array<int, 2> const arr({2, 3});

		BOOST_TEST( arr.size() == 2 );

		multi::array<int, 2> const brr = {};

		BOOST_TEST( brr.size() == 0 );
	}

	{
		multi::array<int, 2> arr({2, 3});

		BOOST_TEST( arr.size() == 2 );

		arr = {};

		BOOST_TEST( arr.size() == 0 );
	}
	{
		// multi::array<int, 2> arr({2, 3});

		// BOOST_TEST( arr.size() == 2 );

		// arr = {{}};  // TODO(correaa) this syntax produces UB in clang

		// BOOST_TEST( arr.size() == 0 );
	}
#endif

	// BOOST_AUTO_TEST_CASE(array_reextent_noop)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		multi::array<double, 2> arr3({2, 3});
		BOOST_TEST(size(arr3) == 2);
		BOOST_TEST(size(arr3[0]) == 3);

		auto* const A_base = arr.base();
		arr.reextent({2, 3});
		BOOST_TEST( arr.num_elements()== 2L*3L );
		BOOST_TEST( arr[1][2] ==  60 );  // reextent preserves values when it can...

		BOOST_TEST( A_base == arr.base() );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_noop_with_init)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		multi::array<double, 2> arr3({2, 3});
		BOOST_TEST(size(arr3) == 2);
		BOOST_TEST(size(arr3[0]) == 3);

		auto* const A_base = arr.base();
		arr.reextent({2, 3}, 990);
		BOOST_TEST( arr.num_elements()== 2L*3L );
		BOOST_TEST( arr[1][2] ==  60 );  // reextent preserves values when it can...

		BOOST_TEST( A_base == arr.base() );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_moved)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		auto* const A_base = arr.base();

		arr = std::move(arr).reextent({2, 3});  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved

		BOOST_TEST( arr.size() == 2 );
		BOOST_TEST( arr.num_elements() == 2L*3L );
		BOOST_TEST( arr.num_elements()== 2L*3L );
		BOOST_TEST(arr[1][2] == 60);  // after move the original elments might not be the same

		BOOST_TEST( A_base == arr.base() );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		auto* const A_base = arr.base();

		arr = std::move(arr).reextent({2, 3});  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved

		BOOST_TEST( arr.num_elements()== 2L*3L );
		BOOST_TEST( arr[1][2] ==  60 );  // after move the original elments might not be the same

		BOOST_TEST( A_base == arr.base() );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial_change_extents)
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr[1][2] = 60;
		BOOST_TEST( arr[1][2] == 60 );  // cppcheck-suppress knownConditionTrueFalse ;

		// auto* const A_base = arr.base();  // to check bellow

		arr = std::move(arr).reextent({4, 5});

		BOOST_TEST( arr.num_elements()== 4L*5L );
		// BOOST_TEST( arr[1][2] !=  6.0 );  // after move the original elements might not be the same, but it is not 100% possible to check
		// BOOST_TEST( A_base != arr.base() );  // after move the base pointer might have changed but it is not 100% to check
	}
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr.reextent({3, 2});

		BOOST_TEST( arr.num_elements() == 6 );
	}
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr = std::move(arr).reextent({3, 2});

		BOOST_TEST( arr.num_elements() == 6 );
	}
	{
		multi::array<int, 2> arr({2, 3});
		BOOST_TEST( arr.num_elements() == 6 );

		arr = std::move(arr).reextent({30, 20});
		std::fill_n(arr.elements().begin(), arr.num_elements(), 99);

		BOOST_TEST( arr.num_elements() == 600 );
	}
	{
		multi::array<int, 2> arr({20, 30});
		BOOST_TEST( arr.num_elements() == 600 );

		arr = std::move(arr).reextent({3, 2});
		std::fill_n(arr.elements().begin(), arr.num_elements(), 99);

		BOOST_TEST( arr.num_elements() == 6 );
	}

	// BOOST_AUTO_TEST_CASE(array_move_clear)
	{
		multi::array<int, 2> arr({2, 3});

		arr = multi::array<int, 2>(extents(arr), 1230);
		BOOST_TEST( arr[1][2] == 1230 );

		arr.clear();
		BOOST_TEST( arr.num_elements() == 0 );
		BOOST_TEST( size(arr) == 0 );

		arr.reextent({5, 4}, 660);
		BOOST_TEST( arr[4][3] == 660 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_1d)
	{
		multi::array<int, 1> arr(multi::extents_t<1>{multi::iextension{10}}, 40);
		BOOST_TEST( size(arr) == 10 );
		BOOST_TEST( arr[9] == 40 );

		arr.reextent(multi::extents_t<1>{multi::iextension{20}});
		BOOST_TEST( size(arr) == 20 );
		BOOST_TEST( arr[9] == 40 );
		// BOOST_TEST( arr[19] == 0.0 );  // impossible to know since it is only sometimes 0.0

		arr.reextent(boost::multi::tuple<int>(22));
		BOOST_TEST( size(arr) == 22 );
		BOOST_TEST( arr[9] == 40 );

		arr.reextent({23});
		BOOST_TEST( size(arr) == 23 );
	}

	// BOOST_AUTO_TEST_CASE(tuple_decomposition)
	{
		boost::multi::tuple<int, int> const tup{1, 2};
		auto [t0, t1] = tup;
		BOOST_TEST( t0 == 1 );
		BOOST_TEST( t1 == 2 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_0D)
	{
		multi::array<int, 0> arr({}, 40);
		arr.reextent(arr.extents());
		BOOST_TEST( *arr.data_elements() == 40 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_1d_with_initialization)
	{
		multi::array<int, 1> arr(multi::extents_t<1>{multi::iextension{10}}, 40);
		BOOST_TEST( size(arr) == 10 );
		BOOST_TEST( arr[9] == 40 );

		arr.reextent(multi::extents_t<1>{multi::iextension{20}}, 80);
		BOOST_TEST( size(arr) == 20 );
		BOOST_TEST( arr[9] == 40 );
		BOOST_TEST( arr[19] == 80 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_2d)
	{
		multi::array<int, 2> arr({10, 20}, 40);
		BOOST_TEST( arr[1][2] == 40 );

		arr.clear();
		BOOST_TEST( arr.num_elements() == 0 );
		BOOST_TEST( size(arr) == 0 );

		arr.reextent({20, 30}, 90);
		BOOST_TEST( arr[1][2] == 90 );
		BOOST_TEST( arr[11][22] == 90 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_2d_with_move)
	{
		multi::array<int, 2> arr = {
			{1, 2, 3},
			{4, 5, 6},
		};
		BOOST_TEST( arr.size() == 2 );

		arr = std::move(arr).reextent({3, 2});

		BOOST_TEST( arr.size() == 3 );
	}

	// BOOST_AUTO_TEST_CASE(array_reextent_2d_array)
	{
		multi::array<int, 2> arr({10, 20}, 40);
		BOOST_TEST( arr[1][2] == 40 );

		arr.clear();
		BOOST_TEST( arr.num_elements() == 0 );
		BOOST_TEST( size(arr) == 0 );
	}

	// BOOST_AUTO_TEST_CASE(array_vector_size)
	{
		std::vector<double> const vec(100);  // std::vector NOLINT(fuchsia-default-arguments-calls)
		{
			multi::array<double, 1> const arr(static_cast<multi::ssize_t>(vec.size()));
			BOOST_TEST( comp_equal(arr.size(), vec.size()) );
		}
		{
			multi::array<double, 1> const arr(multi::iextensions<1>(static_cast<multi::ssize_t>(vec.size())));  // warning: sign-conversion
			BOOST_TEST( comp_equal(arr.size(), vec.size()) );
		}
	}

	// BOOST_AUTO_TEST_CASE(array_iota)
	{
		multi::array<double, 1> const Aarr(10);

		multi::array<multi::array<double, 2>::index, 1> const Barr(Aarr.extent().begin(), Aarr.extent().end());

		BOOST_TEST( Barr[0] == 0 );
		BOOST_TEST( Barr[1] == 1 );
		BOOST_TEST( Barr[9] == 9 );

		multi::array<multi::array<double, 1>::index, 1> const Carr(Aarr.extent());
		BOOST_TEST( Carr[0] == 0 );
		BOOST_TEST( Carr[1] == 1 );
		BOOST_TEST( Carr[9] == 9 );

		multi::array<multi::array<double, 1>::index, 1> const Darr(Aarr.extents());
		BOOST_TEST( Darr.extents() == Aarr.extents() );
	}

#ifndef __INTEL_COMPILER
	// BOOST_AUTO_TEST_CASE(extension_index_op)
	{
		multi::array<double, 2> const Aarr({11, 13});

		auto const Aext = Aarr.extents();

		using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		BOOST_TEST( get<0>(Aext[3][5]) == 3 );
		BOOST_TEST( get<1>(Aext[3][5]) == 5 );

		for(int i = 0; i != 3; ++i) {
			for(int j = 0; j != 5; ++j) {  // NOLINT(altera-unroll-loops)
				auto [ip, jp] = Aext[i][j];
				BOOST_TEST(ip == i);
				BOOST_TEST(jp == j);
			}
		}
	}
#endif
	{
		multi::array<int, 2> arr;

		std::vector<std::vector<int>> varr = {
			std::vector<int>{1, 2},
			std::vector<int>{3, 4}
		};

		arr.reextent({static_cast<multi::array<int, 2>::size_type>(varr.size()), static_cast<multi::array<int, 2>::size_type>(varr[0].size())});
	}
	return boost::report_errors();
}
