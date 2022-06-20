// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array_ref.hpp"

#include <boost/hana/integral_constant.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/serialization/nvp.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


#include<numeric>  // for accumulate

namespace hana = boost::hana;
namespace multi = boost::multi;

template<typename Integral, Integral const N>
struct integral_constant : private hana::integral_constant<Integral, N> {
//	using hana::integral_constant<Integral, n>::integral_constant;
	constexpr explicit operator Integral const&() const {
		return hana::integral_constant<Integral, N>::value;
	}
	integral_constant() = default;
	constexpr explicit integral_constant(Integral const& i) {
		assert(i == N);
	}
	constexpr auto operator==(Integral const& o) const {return static_cast<Integral const&>(*this)==o;}
	constexpr auto operator==(integral_constant const&/*other*/) {return std::true_type{};}
	template<Integral N2, typename = std::enable_if_t<(N2 != N)> >
	constexpr auto operator==(integral_constant<Integral, N2> const&/*other*/) {return std::false_type{};}
	template<Integral N2>
	friend constexpr auto operator+(integral_constant const&/*a*/, integral_constant<Integral, N2> const&/*b*/) {
		return integral_constant<Integral, hana::integral_constant<Integral, N>::value + N2>{};
	}
	template<Integral N2>
	friend constexpr auto operator-(integral_constant const&/*a*/, integral_constant<Integral, N2> const&/*b*/) {
		return integral_constant<Integral, hana::integral_constant<Integral, N>::value - N2>{};
	}
//	constexpr auto operator=(Integral other) -> integral_constant&{assert(other == n); return *this;}
	friend constexpr auto operator>=(Integral const& a, integral_constant const&/*self*/) {return a >= N;}
	friend constexpr auto operator< (Integral const& a, integral_constant const&/*self*/) {return a <  N;}
};

BOOST_AUTO_TEST_CASE(xml_serialization_index_range) {
	std::stringstream ss;
	multi::range<std::ptrdiff_t> const rg{5, 10};
	{
	    boost::archive::xml_oarchive oa(ss);
		oa<< ::boost::serialization::make_nvp("rg", rg);
	}
	{
		boost::archive::xml_iarchive ia(ss);
		multi::range<std::ptrdiff_t> rg2;
		ia>> ::boost::serialization::make_nvp("rg2", rg2);
		BOOST_REQUIRE( rg == rg2 );
	}
}

BOOST_AUTO_TEST_CASE(multi_range) {
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

BOOST_AUTO_TEST_CASE(multi_range_with_hana_literals) {
	static_assert(( integral_constant<int, 1234>{} == 1234 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == 1235 ), "!");
	static_assert(( (integral_constant<int, 1234>{} + integral_constant<int, 1>{})  == integral_constant<int, 1235>{} ), "!");
#if defined(__cpp_deduction_guides) and __cpp_deduction_guides
	static_assert(( multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}.size() == integral_constant<int, 5>{} ), "!");
	static_assert(( size(multi::range{integral_constant<int, 0>{}, integral_constant<int, 5>{}}) == integral_constant<int, 5>{} ), "!");
#endif
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr) {
	BOOST_REQUIRE( multi::extension_t<int>{5} == 5 ); // this is not a constexpr in cuda 10
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> rr{5, 12};
	BOOST_REQUIRE( rr.contains(6) );
	BOOST_REQUIRE( not rr.contains(12) );
}

BOOST_AUTO_TEST_CASE(multi_range2) {
	multi::index_extension x(10);

	BOOST_REQUIRE( *begin(x) == 0 );
	BOOST_REQUIRE( size(x) == 10 );
	BOOST_REQUIRE( x[0] == 0 );
	BOOST_REQUIRE( x[1] == 1 );
	BOOST_REQUIRE( x[9] == 9 );

	auto b = begin(x);
	BOOST_REQUIRE( b[0] == x[0] );
	BOOST_REQUIRE( b[1] == x[1] );

//	static_assert( ranges::forward_iterator< std::decay_t<decltype(b)> > , "!");

	BOOST_REQUIRE( std::accumulate( begin(x), end(x), 0) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

//  #if(__cpp_structured_bindings >= 201606)
//   {
//  	multi::iextensions<3> ies({{0, 3}, {0, 4}, {0, 5}});
//  	BOOST_REQUIRE( std::get<1>(ies).size() == 4 );
//  	auto [is, js, ks] = ies;
//  	BOOST_REQUIRE( is.size() == 3 );
//  }
//  #endif
}

