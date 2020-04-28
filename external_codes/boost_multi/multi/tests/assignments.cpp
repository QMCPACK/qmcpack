#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi assignments"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../tests/../array.hpp"

#include<boost/iterator/transform_iterator.hpp>
#include<boost/functional/hash.hpp>

#include<functional>
#include<iostream>
#include<vector>

namespace multi = boost::multi;
using std::cout;

multi::array_ref<double, 2> make_ref(double* p){return {p, {5, 7}};}

BOOST_AUTO_TEST_CASE(range_assignment){
{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x){return x+1;};
	std::vector<double> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7 );
}
{
	auto r = multi::make_range(5, 10);
	auto f = [](auto x){return x+1;};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 7 );
}
{
	auto r = multi::make_extension_t(10l);
	multi::array<double, 1> v(r.begin(), r.end());
	assert( r.size() == v.size() );
	assert( v[1] = 10 );
}
{
	auto r = multi::make_extension_t(10l);
	auto f = [](auto x){
		std::size_t seed = 1234;
	//	boost::hash_combine(seed, );
		seed ^= boost::hash<multi::index>{}(x) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		return seed/(double)std::numeric_limits<std::size_t>::max();
	};
	multi::array<double, 1> v(
		boost::make_transform_iterator(r.begin(), f), 
		boost::make_transform_iterator(r.end()  , f)
	);

	std::cerr << std::hash<std::size_t>{}(13) << std::endl;
	std::size_t seed = 12349l;
	//	boost::hash_combine(seed, );
//	seed ^= boost::hash<std::size_t>{}(13) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	boost::hash_combine(seed, 13);
	std::cerr << seed << std::endl;
	std::cerr << seed*1./std::numeric_limits<std::size_t>::max() << std::endl;

	assert( v.size() == r.size() );
	assert( v[1] >= 0. );
	assert( v[1] < 1.  );
	assert( std::all_of(begin(v), end(v), [](auto x){
		std::cout << x << std::endl;
		return x >= 0. and x < 1.;
	}) );
}
{
	multi::array<double, 1> v(10);
	auto r = extension(v);
	v.assign(r.begin(), r.end());
	assert( v[1] == 1 );
}
{
	multi::array<double, 1> v(10);
	auto r = extension(v);
	auto f = [](auto x){return x*2;};
	v.assign(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	assert( v[1] == 2 );
}
}

template<class Array> auto rearranged(Array&& arr){
	return std::forward<Array>(arr).unrotated().partitioned(2).transposed().rotated();
}

BOOST_AUTO_TEST_CASE(rearranged_assignment){
	multi::array<double, 4> tmp({14, 14, 7, 4});
	multi::array<double, 5> src({2, 14, 14, 7, 2}); src[0][1][2][3][1] = 99.;

	BOOST_REQUIRE( extensions(rearranged(tmp)) == extensions(src) );

//	BOOST_REQUIRE( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
//	tmp.unrotated().partitioned(2).transposed().rotated() = src;

	rearranged(tmp) = src;
	BOOST_REQUIRE( rearranged(tmp) == src );
	BOOST_REQUIRE( rearranged(tmp) == ~(tmp>>1).partitioned(2)<<1 );

	auto rearranged2 = [](auto&& arr){return std::forward<decltype(arr)>(arr).unrotated().partitioned(2).transposed().rotated();};
	rearranged2(tmp) = src;	
}

BOOST_AUTO_TEST_CASE(rvalue_assignments){
	using complex = std::complex<double>;

	std::vector<double> const v1(200, 99.);
	std::vector<complex> v2(200);
	auto linear1 = [&]{return multi::array_cptr<double, 1>(v1.data(), 200);};
	auto linear2 = [&]{return multi::array_ptr<complex, 1>(v2.data(), 200);};
	*linear2() = *linear1();

}

BOOST_AUTO_TEST_CASE(assignments){
	{
		std::vector<double> v(5*7, 99.);

		multi::array<double, 2> A{{5, 7}, 33.};
		multi::array_ref<double, 2>(v.data(), {5, 7}) = A;
		BOOST_REQUIRE( v[9] == 33. );
		BOOST_REQUIRE( not v.empty() );
		BOOST_REQUIRE( not is_empty(A) );

		multi::array<double, 1> V;
		BOOST_REQUIRE( V.empty() );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		multi::array_ref<double, 2> B{w.data(), {5, 7}};
		make_ref(v.data()) = std::move(B);
		make_ref(v.data()) = std::move(B).sliced(0,5);

		BOOST_REQUIRE( v[9] == 33. );
	}
	{
		std::vector<double> v(5*7, 99.), w(5*7, 33.);

		make_ref(v.data()) = make_ref(w.data());

		BOOST_REQUIRE( v[9] == 33. );
	}
}

