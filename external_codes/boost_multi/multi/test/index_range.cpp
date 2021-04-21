#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -DMULTI_ACCESS_NDEBUG -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include <boost/hana/integral_constant.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include<cassert>
#include<iostream>
#include<vector>

namespace hana = boost::hana;

////https://stackoverflow.com/a/35110453/225186
//template<class T>constexpr std::remove_reference_t<T> const_aux(T&&t){return t;}
//template<bool b> struct logic_assert_aux;
//template<> struct logic_assert_aux<true>{
//	template<class T> static constexpr void _(T&& cond){static_assert(cond, "!");}
//};
//template<> struct logic_assert_aux<false>{
//	template<class T> static constexpr void _(T&& cond){assert(cond && "!");}
//};
//#if (not defined(__INTEL_COMPILER)) and (not defined(__NVCC__))
//#define logic_assert(ConD, MsG) logic_assert_aux<noexcept(const_aux(ConD))>::_(ConD);
//#else
//#define logic_assert(ConD, MsG) assert(ConD)
//#endif
#if 0
if constexpr(noexcept(const_aux(ConD))) static_assert(ConD, MsG);\
	else assert(ConD && MsG);
#endif

template<typename Integral, Integral const n>
struct integral_constant : private hana::integral_constant<Integral, n>{
//	using hana::integral_constant<Integral, n>::integral_constant;
	constexpr explicit operator Integral const&() const{
		return hana::integral_constant<Integral, n>::value;
	}
	integral_constant() = default;
	constexpr explicit integral_constant(Integral const& i){
		assert(i == n);
	}
	constexpr auto operator==(Integral const& o) const{return static_cast<Integral const&>(*this)==o;}
	constexpr auto operator==(integral_constant const&/*other*/){return std::true_type{};}
	template<Integral n2, typename = std::enable_if_t<n2!=n> >
	constexpr auto operator==(integral_constant<Integral, n2> const&/*other*/){return std::false_type{};}
	template<Integral n2>
	friend constexpr auto operator+(integral_constant const&/*a*/, integral_constant<Integral, n2> const&/*b*/){
		return integral_constant<Integral, hana::integral_constant<Integral, n>::value + n2>{};
	}
	template<Integral n2>
	friend constexpr auto operator-(integral_constant const&/*a*/, integral_constant<Integral, n2> const&/*b*/){
		return integral_constant<Integral, hana::integral_constant<Integral, n>::value - n2>{};
	}
//	constexpr auto operator=(Integral other) -> integral_constant&{assert(other == n); return *this;}
	friend constexpr auto operator>=(Integral const& a, integral_constant const&/*self*/){return a >= n;}
	friend constexpr auto operator<(Integral const& a, integral_constant const&/*self*/){return a < n;}
};

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_range){

#if defined(__cpp_deduction_guides) and __cpp_deduction_guides
	BOOST_REQUIRE(( multi::range{5, 5}.empty() ));
#else
	BOOST_REQUIRE(( multi::range<std::ptrdiff_t>{5, 5}.empty() ));
#endif
{
	auto r = multi::range<std::ptrdiff_t>{5, 10};
	std::vector<double> v(r.begin(), r.end());
	BOOST_REQUIRE( v[1] == 6 );
}
{
	auto r = multi::range<std::ptrdiff_t>{5, 10};
	auto f = [](auto x){return x+1;};
	std::vector<double> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);
	BOOST_REQUIRE( v[1] == 7 );
}
}

BOOST_AUTO_TEST_CASE(multi_range_with_hana_literals){
//	using namespace hana::literals; // contains the _c suffix
	
	static_assert(( integral_constant<int, 1234>{} == 1234 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == 1235 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == integral_constant<int, 1235>{} ), "!");
#if defined(__cpp_deduction_guides) and __cpp_deduction_guides
	static_assert(( multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}.size() == integral_constant<int, 5>{} ), "!");
	static_assert(( size(multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}) == integral_constant<int, 5>{} ), "!");
//	integral_constant<int, 5> five; five = 5;
#endif
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr){
	BOOST_REQUIRE( multi::extension_t<int>{5} == 5 ); // this is not a constexpr in cuda 10
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> rr{5, 12};
	BOOST_REQUIRE( rr.contains(6) );
	BOOST_REQUIRE( not rr.contains(12) );
}

