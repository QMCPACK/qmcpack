// Copyright 2019-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, apply, operator==, layout_t

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for fill
// #include <complex>    // for complex
#include <cstddef>  // for size_t
// #include <limits>   // for std::numeric_limits
// #include <iterator>  // for size
#include <memory>  // for std::allocator  // IWYU pragma: keep
// #include <stdexcept>  // for std::runtime_error
// IWYU pragma: no_include <type_traits>  // for decay_t
#include <utility>  // for move
#include <vector>   // for vector, allocator

namespace multi = boost::multi;

namespace {

struct element_test {  // NOLINT(cppcoreguidelines-special-member-functions,hicpp-special-member-functions) some of these are trully unused
	static std::size_t construct_count;
	static std::size_t destruct_count;

 private:
	int value_;

 public:
	explicit element_test(int val = 0) noexcept : value_(val) { ++construct_count; }
	~element_test() noexcept { ++destruct_count; }
	element_test(element_test const& other) noexcept : value_(other.value_) { ++construct_count; }  // NOLINT(bugprone-copy-constructor-init)
	// element_test(element_test&& other) noexcept : value_(other.value_) { ++construct_count; }
	auto operator=(element_test const& /*other*/) noexcept -> element_test& = default;
	// auto operator=(element_test&& /*other*/) noexcept -> element_test&      = default;
};
std::size_t element_test::construct_count = 0;
std::size_t element_test::destruct_count  = 0;

constexpr auto make_ref(int* ptr) {
	return multi::array_ref<int, 2>(ptr, {5, 7});
}

template<class T, class Allocator>
auto eye(multi::extents_t<2> exts, Allocator const& alloc) {
	multi::array<T, 2, Allocator> ret(exts, 0, alloc);
	std::fill(ret.diagonal().begin(), ret.diagonal().end(), T{1});
	return ret;
}

template<class T>
auto eye(multi::extents_t<2> exts) { return eye<T>(exts, std::allocator<T>{}); }

}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(equality_1D)
	{
		multi::array<int, 1> arr  = {10, 20, 30};
		multi::array<int, 1> arr2 = {10, 20, 30};

		BOOST_TEST(    arr == arr2  );
		BOOST_TEST( !(arr != arr2) );

		BOOST_TEST(    arr() == arr2()  );
		BOOST_TEST( !(arr() != arr2()) );
	}

	// BOOST_AUTO_TEST_CASE(equality_2D)
	{
		multi::array<int, 2> arr = {
			{10, 20, 30},
			{40, 50, 60},
		};
		multi::array<int, 2> arr2 = {
			{10, 20, 30},
			{40, 50, 60},
		};

		BOOST_TEST( arr == arr2 );
		BOOST_TEST( !(arr != arr2) );

		BOOST_TEST( arr() == arr2() );
		BOOST_TEST( !(arr() != arr2()) );

		BOOST_TEST( arr[0] == arr2[0] );
		BOOST_TEST( !(arr[0] != arr2[0]) );
	}

	// BOOST_AUTO_TEST_CASE(multi_copy_move)
	{
		multi::array<double, 2> arr({3, 3}, 0.0);
		multi::array<double, 2> arr2 = arr;
		BOOST_TEST( arr == arr2 );  // cppcheck-suppress knownConditionTrueFalse ;

		auto* arr_data = arr.data_elements();

		multi::array<double, 2> arr3 = std::move(arr);

		BOOST_TEST( arr3.data_elements() == arr_data );

		multi::array<double, 2> const arr4(std::move(arr2));
		BOOST_TEST( arr4.size() == 3 );
	}

	// BOOST_AUTO_TEST_CASE(range_assignment)
	{
		{
			auto const ext = multi::make_extension_t(10L);

			multi::array<multi::ssize_t, 1> vec(ext.begin(), ext.end());

			BOOST_TEST( ext.size() == vec.size() );
			BOOST_TEST( vec[1] == 1L );
		}
		// {
		// 	multi::array<multi::size_t, 1> vec(multi::extents_t<1>{multi::iextension{10}});

		// 	auto const ext = extension(vec);

		// 	vec.assign(ext.begin(), ext.end());
		// 	BOOST_TEST( vec[1] == 1 );
		// }
	}

#ifndef MULL
	// BOOST_AUTO_TEST_CASE(rearranged_assignment)
	{
		auto const ext5 = multi::extents_t<5>{2, 14, 14, 7, 2};

		[[maybe_unused]] auto const ext52 = ext5;

		[[maybe_unused]] multi::array<int, 5> const src_test(ext5);

		multi::array<int, 5> src({2, 14, 14, 7, 2});

		src[0][1][2][3][1] = 99;

		BOOST_TEST( src[0][1][2][3][1] == 99 );  // cppcheck-suppress knownConditionTrueFalse ; for test

		// BOOST_TEST( tmp.unrotated().partitioned(2).transposed().rotated().extensions() == src.extensions() );
		// BOOST_TEST( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
	}
#endif

	// BOOST_AUTO_TEST_CASE(rearranged_assignment_resize)
	{
		multi::array<double, 2> const arrA({4, 5});
		multi::array<double, 2>       arrB({2, 3});

		arrB = arrA;
		BOOST_TEST( arrB.size() == 4 );
	}

#ifndef _MSC_VER  // TODO(correaa) fix
				  // seems to produce a deterministic divide by zero
				  // Assertion failed: stride_ != 0, file D:\a\boost-multi\boost-root\boost/multi/detail/layout.hpp, line 767
				  // D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : error C2220: the following warning is treated as an error
				  // D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : warning C4723: potential divide by 0
				  // D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : warning C4723: potential divide by 0

	// BOOST_AUTO_TEST_CASE(rvalue_assignments) {
	//  using complex = std::complex<double>;

	//  std::vector<double> const vec1(200, 99.0);  // NOLINT(fuchsia-default-arguments-calls)
	//  std::vector<complex>      vec2(200);        // NOLINT(fuchsia-default-arguments-calls)

	//  auto linear1 = [&] {
	//      return multi::array_cptr<double, 1>(vec1.data(), 200);
	//  };
	//  auto linear2 = [&] {
	//      return multi::array_ptr<complex, 1>(vec2.data(), 200);
	//  };

	//  *linear2() = *linear1();
	// }
#endif

	// BOOST_AUTO_TEST_CASE(assignments)
	{
		{
			std::vector<int> vec(static_cast<std::size_t>(5 * 7), 99);  // NOLINT(fuchsia-default-arguments-calls)

			constexpr int val = 33;

			multi::array<int, 2> const arr({5, 7}, val);
			multi::array_ref<int, 2>(vec.data(), {5, 7}) = arr();  // arr() is a subarray

			BOOST_TEST( vec[9] == val );  // cppcheck-suppress containerOutOfBounds ;
			BOOST_TEST( !vec.empty() );
			BOOST_TEST( !arr.is_empty() );
		}
		{
			std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
			std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

			auto Bp = &multi::array_ref<int, 2>({5, 7}, wec.data());
			// multi::array_ptr<int, 2> const Bp(wec.data(), {5, 7});

			make_ref(vec.data()) = *Bp;  // cppcheck-suppress danglingTemporaryLifetime

			auto&& mref = make_ref(vec.data());
			// mref        = (*Bp).sliced(0, 5);

			mref = Bp->sliced(0, 5);

			// make_ref(vec.data()) = (*Bp).sliced(0, 5);
			make_ref(vec.data()) = Bp->sliced(0, 5);

			BOOST_TEST( vec[9] == 33 );
		}
		{
			std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
			std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

			make_ref(vec.data()) = make_ref(wec.data());

			BOOST_TEST( vec[9] == 33 );
		}
	}

	// BOOST_AUTO_TEST_CASE(assigment_temporary)
	{
		multi::array<int, 2> Id = eye<int>(multi::extents_t<2>({3, 3}));
		BOOST_TEST( Id == eye<double>({3, 3}) );
		BOOST_TEST( Id[1][1] == 1 );
		BOOST_TEST( Id[1][0] == 0 );
	}

	// BOOST_AUTO_TEST_CASE(self_assignment)
	{
		// NOLINTBEGIN(fuchsia-default-arguments-calls)
		multi::dynamic_array<std::vector<int>, 2> arr = {
			{std::vector<int>(10, 1), std::vector<int>(20, 2)},
			{std::vector<int>(30, 3), std::vector<int>(40, 4)},
		};
		// NOLINTEND(fuchsia-default-arguments-calls)

		BOOST_TEST( arr[1][1] == std::vector<int>(40, 4) );  // NOLINT(fuchsia-default-arguments-calls)
		auto* loc = &arr[1][1][5];

		auto const* arr_ptr = std::addressof(arr);

		arr = *arr_ptr;

		BOOST_TEST( arr[1][1] == std::vector<int>(40, 4) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( &arr[1][1][5] == loc );
	}

	// BOOST_AUTO_TEST_CASE(dynamic_array_move)
	{
		// NOLINTBEGIN(fuchsia-default-arguments-calls)
		multi::dynamic_array<std::vector<int>, 2> arr = {
			{std::vector<int>(10, 1), std::vector<int>(20, 2)},
			{std::vector<int>(30, 3), std::vector<int>(40, 4)},
		};
		// NOLINTEND(fuchsia-default-arguments-calls)

		BOOST_TEST( arr[1][1] == std::vector<int>(40, 4) );  // NOLINT(fuchsia-default-arguments-calls)

		multi::dynamic_array<std::vector<int>, 2> arr2(std::move(arr));
		BOOST_TEST( arr2[1][1] == std::vector<int>(40, 4) );  // NOLINT(fuchsia-default-arguments-calls)

		// cppcheck-suppress accessMoved ;
		BOOST_TEST( arr[0][0].empty() );  // NOLINT(clang-analyzer-cplusplus.Move,fuchsia-default-arguments-calls,bugprone-use-after-move,hicpp-invalid-access-moved)
		BOOST_TEST( arr[1][1].empty() );  // NOLINT(clang-analyzer-cplusplus.Move,fuchsia-default-arguments-calls,bugprone-use-after-move,hicpp-invalid-access-moved)
	}

	// BOOST_AUTO_TEST_CASE(move_assignment_no_leak)
	{
		element_test::construct_count = 0;
		element_test::destruct_count  = 0;
		{
			multi::array<element_test, 2> arr1({3, 4}, element_test(1));
			multi::array<element_test, 2> arr2({2, 2}, element_test(2));

			arr1 = std::move(arr2);  // move assignment: should destroy arr1's old buffer, steal arr2's
		}

		BOOST_TEST( element_test::construct_count == element_test::destruct_count );  // no leaks
	}

	// BOOST_AUTO_TEST_CASE(reextent_no_leak)
	{
		element_test::construct_count = 0;
		element_test::destruct_count  = 0;
		{
			multi::array<element_test, 2> arr({3, 4}, element_test(1));
			arr.reextent({5, 5});  // resize to larger dimensions: allocates new buffer, constructs new, copies intersection, destroys old
		}

		BOOST_TEST( element_test::construct_count == element_test::destruct_count );  // no leaks
	}

	return boost::report_errors();
}
