// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for array, implicit_cast, explicit_cast

#include <type_traits>  // for is_same_v, is_same

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		static_assert(std::is_trivially_default_constructible_v<multi::array<double, 0>::cursor>);
		static_assert(std::is_trivially_copy_constructible_v   <multi::array<double, 0>::cursor>);
		static_assert(std::is_trivially_copy_assignable_v      <multi::array<double, 0>::cursor>);
	}

	// BOOST_AUTO_TEST_CASE(iterator_1d) 
	{
		static_assert(std::is_trivially_default_constructible_v<multi::array<double, 1>::cursor>);
		static_assert(std::is_trivially_copy_constructible_v   <multi::array<double, 1>::cursor>);
		static_assert(std::is_trivially_copy_assignable_v      <multi::array<double, 1>::cursor>);
	}
	// BOOST_AUTO_TEST_CASE(iterator_2d) 
	{
		static_assert(std::is_trivially_default_constructible_v<multi::array<double, 2>::cursor>);
		static_assert(std::is_trivially_copy_constructible_v   <multi::array<double, 2>::cursor>);
		static_assert(std::is_trivially_copy_assignable_v      <multi::array<double, 2>::cursor>);

		// {
		//  multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, 99.0);
		//  BOOST_TEST( arr.size() == 100 );
		//  BOOST_TEST( arr.begin() < arr.end() );
		//  BOOST_TEST( arr.end() - arr.begin() == arr.size() );

		//  multi::array<double, 1>::const_iterator const cbarr = arr.cbegin();
		//  multi::array<double, 1>::iterator             barr  = arr.begin();

		//  [[maybe_unused]] multi::array<double, 1>::const_iterator const cbarr3{barr};

		//  BOOST_TEST(  barr == cbarr );  // problem in C++20
		//  BOOST_TEST( cbarr ==  barr );  // problem in C++20

		//  barr += 1;
		//  barr -= 1;
		//  BOOST_TEST( cbarr == barr );

		//  multi::array<double, 1>::const_iterator const cbarr2 = begin(arr);
		//  BOOST_TEST( cbarr2 == cbarr );
		// }
		// {
		//  multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, 99.0);
		//  BOOST_TEST( arr.size() == 100 );
		//  BOOST_TEST( arr.begin() < arr.end() );

		//  auto                                          arr2 = arr.begin();
		//  multi::array<double, 1>::const_iterator const cbb  = arr2;
		//  BOOST_TEST( cbb == arr2 );
		//  // BOOST_TEST( arr2 == cbb );  // TODO(correaa) problem in C++20
		// }
		// {
		//  multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, 99.0);
		//  BOOST_TEST( arr.size() == 100 );
		//  BOOST_TEST( arr.begin() < arr.end() );

		//  auto const arrend  = arr.end();
		//  auto const arrlast = arrend - 1;

		//  BOOST_TEST( arrlast + 1 == arrend );
		// }
	}

	// BOOST_AUTO_TEST_CASE(iterator_2d) {
	//  BOOST_TEST((std::is_trivially_copy_constructible_v   <multi::layout_t<2>>));
	//  BOOST_TEST((std::is_trivially_copy_assignable_v      <multi::layout_t<2>>));
	//  BOOST_TEST((std::is_trivially_default_constructible_v<multi::layout_t<2>>));

	//  BOOST_TEST((std::is_trivially_copy_constructible_v   <multi::subarray_ptr<double, 2>>));
	//  BOOST_TEST((std::is_trivially_copy_assignable_v      <multi::subarray_ptr<double, 2>>));
	//  BOOST_TEST((std::is_trivially_default_constructible_v<multi::subarray_ptr<double, 2>>));

	//  // BOOST_TEST((std::is_trivially_default_constructible_v<multi::array<double, 2>::iterator>));  // TODO(correaa)
	//  BOOST_TEST((std::is_trivially_copy_constructible_v   <multi::array<double, 2>::iterator>));
	//  BOOST_TEST((std::is_trivially_copy_assignable_v      <multi::array<double, 2>::iterator>));

	//  {
	//      multi::array<double, 2> const arr({120, 140}, 99.0);

	//      BOOST_TEST( arr.size() == 120 );
	//      BOOST_TEST( arr.cbegin() < arr.cend() );
	//      BOOST_TEST( arr.cend() - arr.cbegin() == arr.size() );

	//      using iter = multi::array<double, 2>::iterator;
	//      static_assert(std::is_same_v<iter::element, double>);
	//      static_assert(std::is_same_v<iter::value_type, multi::array<double, 1>>);
	//      static_assert(std::is_same_v<iter::reference, multi::subarray<double, 1>>);
	//      static_assert(std::is_same_v<iter::element_ptr, double*>);

	//      using citer = multi::array<double, 2>::const_iterator;
	//      static_assert(std::is_same_v<citer::element, double>);
	//      static_assert(std::is_same_v<citer::value_type, multi::array<double, 1>>);

	//      static_assert(std::is_same_v<citer::reference, multi::const_subarray<double, 1>>);
	//      static_assert(std::is_same_v<citer::element_ptr, double*>);

	//      auto const arrend  = arr.end();
	//      auto const arrlast = arrend - 1;

	//      BOOST_TEST( arrlast + 1 == arrend );
	//  }
	//  {
	//      std::vector<double>         vec(10000);  // std::vector NOLINT(fuchsia-default-arguments-calls)
	//      multi::array_ref<double, 2> arr(vec.data(), {100, 100});
	//      BOOST_TEST(arr.size() == 100);
	//      begin(arr)[4][3] = 2.0;
	//  }
	// }

	// BOOST_AUTO_TEST_CASE(iterator_interface) {
	//  multi::array<int, 3> arr = {
	//      { {12, 11},  {24, 10}},
	//      {{112, 30}, {344, 40}},
	//      { {12, 11},  {24, 10}}
	//  };

	//  BOOST_TEST( size(arr) == 3 );
	//  BOOST_TEST( size(arr[0]) == 2 );
	//  BOOST_TEST( size(arr[0][0]) == 2 );
	//  BOOST_TEST( arr[0][0][1] == 11 );

	//  BOOST_TEST( begin(arr) < end(arr) );
	//  BOOST_TEST( cbegin(arr) < cend(arr) );
	//  BOOST_TEST( begin(arr[0]) < end(arr[0]) );
	//  BOOST_TEST( begin(arr[0]) < end(arr[0]) );

	//  //  BOOST_TEST(( multi::array<double, 3>::reverse_iterator {A.begin()} == rend(A) ));

	//  //  BOOST_TEST( rbegin(A) < rend(A) );

	//  BOOST_TEST( end(arr) - begin(arr) == size(arr) );
	//  //  BOOST_TEST( rend(A) - rbegin(A) == size(A) );

	//  BOOST_TEST( size(*begin(arr)   ) == 2 );
	//  BOOST_TEST( size( begin(arr)[1]) == 2 );

	//  BOOST_TEST( &(arr[1][1].begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	//  BOOST_TEST( &arr[0][1][0] == &arr[0][1][0] );
	//  BOOST_TEST( &((*arr.begin())[1][0]) == &arr[0][1][0] );

	//  BOOST_TEST( &((*arr.begin()).operator[](1)[0]) == &arr[0][1][0] );
	//  BOOST_TEST( &(  arr.begin()->operator[](1)[0]) == &arr[0][1][0] );

	//  BOOST_TEST( &((*arr.begin()).operator[](1).begin()[0]) == &arr[0][1][0] );  // NOLINT(readability-container-data-pointer) test access
	//  BOOST_TEST( &(  arr.begin()->operator[](1).begin()[0]) == &arr[0][1][0] );    // NOLINT(readability-container-data-pointer) test access

	//  BOOST_TEST( &((*(arr.begin()+1)).operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	//  BOOST_TEST( &(  (arr.begin()+1)->operator[](1).begin()[0]) == &arr[1][1][0] );    // NOLINT(readability-container-data-pointer) test access

	//  BOOST_TEST( &((*(begin(arr)+1)).operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	//  BOOST_TEST( &((  begin(arr)+1)->operator[](1).begin()[0]) == &arr[1][1][0] );    // NOLINT(readability-container-data-pointer) test access

	//  BOOST_TEST( &((*(cbegin(arr)+1)).operator[](1).begin()[0]) == &arr[1][1][0] );  // NOLINT(readability-container-data-pointer) test access
	//  BOOST_TEST( &((  cbegin(arr)+1)->operator[](1).begin()[0]) == &arr[1][1][0] );    // NOLINT(readability-container-data-pointer) test access
	// }

	// BOOST_AUTO_TEST_CASE(iterator_semantics) {
	//  multi::array<double, 3> arr = {
	//      { {1.2, 1.1},  {2.4, 1.0}},
	//      {{11.2, 3.0}, {34.4, 4.0}},
	//      { {1.2, 1.1},  {2.4, 1.0}}
	//  };

	//  multi::array<double, 3>::iterator it;
	//  // BOOST_TEST(( multi::array<double, 3>::iterator{} == it ));  // `it` is uninitialized
	//  // BOOST_TEST(( it == multi::array<double, 3>::iterator{} ));

	//  it = begin(arr);
	//  BOOST_TEST( it == begin(arr) );

	//  it += 1;
	//  it -= 1;
	//  BOOST_TEST( it == begin(arr) );

	//  auto const& arrc = arr();
	//  BOOST_TEST( &arrc[0][0][0] == &arr[0][0][0] );

	//  auto const& arrc2 = arr();

	//  BOOST_TEST( arrc.addressof() == arrc2.addressof() );  // BOOST_TEST( &arrc == &arrc2 );

	//  multi::array<double, 3>::iterator const it2 = begin(arr);
	//  BOOST_TEST(it == it2);

	//  it = end(arr);
	//  BOOST_TEST(it != it2);
	//  BOOST_TEST(it > it2);

	//  multi::array<double, 3>::iterator const it3{it};
	//  BOOST_TEST( it3 == it );

	//  static_assert(std::is_same<multi::array<double, 3>::iterator::element_ptr, double*>{});

	//  // cit = it3;
	//  // BOOST_REQUIRE( cit == it3 );  // TODO(correaa)
	//  // BOOST_REQUIRE( it3 == cit );  // TODO(correaa)

	//  // cit = it3;
	//  // BOOST_TEST( cit == it3 );  // TODO(correaa)
	//  // BOOST_TEST( it3 == cit );  // TODO(correaa)

	//  BOOST_TEST( &arr[0][2][1] == &begin(arr)[0][2][1] );

	//  [[maybe_unused]] multi::array<double, 3>::const_iterator const cit2 = it3;

	//  static_assert(decltype(begin(arr))::rank_v == 3, "!");
	//  static_assert(decltype(begin(arr))::rank{} == 3, "!");

	//  // auto&& ref = multi::ref(begin(arr), end(arr));

	//  // BOOST_TEST( arr.base() == ref.base() );
	//  // BOOST_TEST(  arr[0][2][1] ==  ref[0][2][1] );
	//  // BOOST_TEST( &arr[0][2][1] == &ref[0][2][1] );
	//  // BOOST_TEST( arr.layout().stride() == ref.layout().stride());
	//  // BOOST_TEST( arr.layout().offset() == ref.layout().offset());
	//  // BOOST_TEST( arr.layout().nelems() == ref.layout().nelems());

	//  // BOOST_TEST( arr.num_elements() == ref.num_elements() );
	//  // BOOST_TEST( arr.stride() == ref.stride() );
	//  // BOOST_TEST( arr.layout() == ref.layout() );

	//  // BOOST_TEST( &multi::ref(begin(arr), end(arr)) == &arr );
	// }

	// BOOST_AUTO_TEST_CASE(iterator_arrow_operator) {
	//  // NOLINTBEGIN(fuchsia-default-arguments-calls) std::string has a default constructor
	//  multi::array<std::string, 2> arr = {
	//      {"00", "01"},
	//      {"10", "11"},
	//      {"20", "21"}
	//  };
	//  // NOLINTEND(fuchsia-default-arguments-calls)

	//  BOOST_TEST( arr[1][0] == "10" );

	//  BOOST_TEST( std::is_sorted(begin(arr), end(arr)) );                      // sorted by rows
	//  BOOST_TEST( std::is_sorted(begin(arr.rotated()), end(arr.rotated())) );  // sorted by cols

	//  BOOST_TEST( (*begin( arr           )).size() == arr[0].size() );
	//  BOOST_TEST(   begin( arr           )->size() == arr[0].size() );

	//  BOOST_TEST( (*begin( arr.rotated() )).size() == arr.size() );
	//  BOOST_TEST(   begin( arr.rotated() )->size() == arr.size() );

	//  BOOST_TEST( &((*begin( arr           )).operator[](1)) == &(arr[0][1]) );
	//  BOOST_TEST( &(  begin( arr           )->operator[](1)) == &(arr[0][1]) );

	//  BOOST_TEST( &((*begin( arr.rotated() )).operator[](1)) == &(arr[1][0]) );
	//  BOOST_TEST( &(  begin( arr.rotated() )->operator[](1)) == &(arr[1][0]) );
	// }

	// BOOST_AUTO_TEST_CASE(index_range_iteration) {
	//  multi::index_range irng(0, 5);  // semiopen interval
	//  std::ostringstream out;
	//  std::copy(irng.begin(), irng.end(), std::ostream_iterator<multi::index_range::value_type>{out, ","});
	//  BOOST_TEST_EQ(out.str(), std::string{"0,1,2,3,4,"});  // NOLINT(fuchsia-default-arguments-calls)

	//  BOOST_TEST( std::accumulate(begin(irng), end(irng), static_cast<multi::index_range::value_type>(0U)) == irng.size()*(irng.size()-1)/2 );

	//  auto const sum_of_cubes = [](auto&& acc, auto const& elem) {
	//      return std::forward<decltype(acc)>(acc) + elem * elem * elem;
	//  };
	//  BOOST_TEST( std::accumulate(begin(irng), end(irng), multi::index_range::value_type{}, sum_of_cubes) > 0 );
	// }

	// BOOST_AUTO_TEST_CASE(multi_reverse_iterator_1D) {
	//  multi::array<double, 1> arr(100, 66.0);
	//  BOOST_TEST( &arr[99] == &*std::make_reverse_iterator(arr.end()) );

	//  auto rbegin = std::make_reverse_iterator(arr.end());
	//  rbegin += 100;
	//  multi::array<double, 1>::iterator const begin{rbegin.base()};
	//  BOOST_TEST( begin  == arr.begin() );
	// }

	// BOOST_AUTO_TEST_CASE(multi_reverse_iterator_2D) {
	//  multi::array<int, 2> arr = {
	//      {  10,   20},
	//      { 100,  200},
	//      {1000, 2000}
	//  };
	//  BOOST_TEST( (*arr.begin())[1] == 20 );
	//  BOOST_TEST( arr.begin()->operator[](1) == 20 );

	//  auto rbegin = std::make_reverse_iterator(arr.end());

	//  BOOST_TEST( (*rbegin)[1] == 2000 );

	//  BOOST_TEST( arr.begin()   < arr.begin() + 1 );
	//  BOOST_TEST( arr.end() - 1 < arr.end()       );
	// }

	return boost::report_errors();
}
