#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -DMULTI_ACCESS_NDEBUG -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
#include <boost/test/execution_monitor.hpp>  // for boost::execution_exception

#include "../array.hpp"

#include <boost/hana/integral_constant.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include<vector>
#include<cassert>
#include<iostream>

namespace multi = boost::multi;
namespace hana = boost::hana;

using std::cout;

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
	constexpr operator Integral const&() const{
		return hana::integral_constant<Integral, n>::value;
	}
	constexpr integral_constant(){}
	explicit constexpr integral_constant(Integral const& i){
		assert(i == n);
	}
	constexpr auto operator==(Integral const& o) const{return static_cast<Integral const&>(*this)==o;}
	constexpr std::true_type operator==(integral_constant const&){return {};}
	template<Integral n2, typename = std::enable_if_t<n2!=n> >
	constexpr std::false_type operator==(integral_constant<Integral, n2> const&){return {};}
	template<Integral n2>
	friend constexpr auto operator+(integral_constant const&, integral_constant<Integral, n2> const&){
		return integral_constant<Integral, hana::integral_constant<Integral, n>::value + n2>{};
	}
	template<Integral n2>
	friend constexpr auto operator-(integral_constant const&, integral_constant<Integral, n2> const&){
		return integral_constant<Integral, hana::integral_constant<Integral, n>::value - n2>{};
	}
	constexpr integral_constant& operator=(Integral other){assert(other == n); return *this;}
	friend constexpr auto operator>=(Integral const& a, integral_constant const&){return a >= n;}
	friend constexpr auto operator<(Integral const& a, integral_constant const&){return a < n;}
};

template<class T> T& evoke(T&& t){return t;}

template<class T> void what(T&&) = delete;

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_range){

#if defined(__cpp_deduction_guides) and __cpp_deduction_guides
	BOOST_REQUIRE(( multi::range{5, 5}.size() == 0 ));
#else
	BOOST_REQUIRE(( multi::range<std::ptrdiff_t>{5, 5}.size() == 0 ));
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
	using namespace hana::literals; // contains the _c suffix
	static_assert(( integral_constant<int, 1234>{} == 1234 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == 1235 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == integral_constant<int, 1235>{} ), "!");
#if defined(__cpp_deduction_guides) and __cpp_deduction_guides
	static_assert(( multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}.size() == integral_constant<int, 5>{} ), "!");
	static_assert(( size(multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}) == integral_constant<int, 5>{} ), "!");
	integral_constant<int, 5> five; five = 5;
#endif
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr){
	BOOST_REQUIRE( multi::extension_t<int>{5} == 5 ); // this is not a constexpr in cuda 10
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));
//	static_assert(( multi::extension_t{integral_constant<int, 5>{}, integral_constant<int, 12>{}}.contains(10) ), "!");
//	static_assert(( multi::extension_t{integral_constant<int, 5>{}, integral_constant<int, 12>{}}.contains(integral_constant<int, 10>{}) ), "!");


//	logic_assert( size(multi::range<int>{5, 5}) == 0 , "!");
//	static_assert( is_empty(multi::range<int>{5, 5}) , "!");
	
//	logic_assert( size(multi::range<int>{}) == 0 , "!");
//	logic_assert( empty(multi::range<int>{}) , "!");
	
	for(auto const& i : multi::range<int>{5, 12}) cout<< i <<' ';
	cout<<'\n';
	
//	static_assert(
//		empty(intersection(multi::range<int>{5, 12}, multi::range<int>{14, 16}))// == multi::range<int>{}
//	);
	cout<< intersection(multi::range<int>{5, 12}, multi::range<int>{14, 16}) <<'\n';
	cout<< intersection(multi::range<int>{5, 12}, multi::range<int>{14, 16}) <<'\n';

//	for(auto const& i : intersection(multi::range<int>{5, 12}, multi::range<int>{8, 16})) cout<< i <<' ';
//	cout <<'\n';
	
	multi::range<int> rr{5, 12};
	assert( rr.contains(6) );
	assert( not rr.contains(12) );
	for(auto it = rr.begin(); it != rr.end(); ++it) cout<< *it <<' ';
	cout<<'\n';
	for(auto it = rr.rbegin(); it != rr.rend(); ++it) cout<< *it <<' ';
	cout<<'\n';



//	cout<< *rr.rbegin() <<'\n';
//	for(auto it = rr.rbegin(); it != rr.rend(); ++it) cout<< *it <<' ';
//	cout <<'\n';
	
//	multi::extension<int> ei{5, 10}; 
#if 0
	std::iterator_traits<multi::range<multi::index>::const_iterator>::value_type p = multi::index{4};

	{
		multi::index_range ir{5, 10};
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		std::vector<multi::index_range::value_type> v(5);
		copy(begin(ir), end(ir), begin(v));
		assert(v[0] == 5);
		for(auto& i : ir) cout << i << ' ';
		cout << '\n';
		auto f = ir.find(6);
		cerr << "*f " << *f << '\n';
		assert(*f == 6);
		using std::find;
		auto f2 = find(ir.begin(), ir.end(), 12);
		assert(f2 == ir.end());
		auto f3 = find(ir.begin(), ir.end(), 2);
		assert(f3 == ir.end());
	}
/*	{
		multi::strided_index_range ir{6, 12, 2};
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		std::vector<multi::index_range::value_type> v(5);
		copy(begin(ir), end(ir), begin(v));
		assert( v[0] == 6 );
		assert( v[1] == 8 );
		for(auto& i : ir) cout << i <<' ';
		cout <<'\n';
	}*/
	{
		multi::index_range ir(5);
		cout << ir << " = {" << format(index_ % ", ", ir) << "}\n";
		assert(*begin(ir) == 5);
		assert(ir.front() == 5);
		assert(ir.back() == 5);
	}
	{
		multi::index_range ir; // partially formed
		ir = multi::index_range{8, 8};
		assert(ir.empty());
	}
	{
		multi::index_range ir = {};
		assert(ir.empty());
	}
	{
		multi::index_extension ie(5);
		cout << ie << " = {" << format(index_ % ", ", ie) << "}";
	}
#endif

}

