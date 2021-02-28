#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0x -lstdc++fs -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi interacting with Boost Interprocess"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include <boost/interprocess/managed_mapped_file.hpp>

#include<set>
#include<experimental/filesystem> // remove_all c++17

namespace bip = boost::interprocess;

using manager = bip::managed_mapped_file;

using std::experimental::filesystem::path;
using std::experimental::filesystem::remove;

template<class T> using mallocator = bip::allocator<T, manager::segment_manager>;
static auto get_allocator(manager& m){return m.get_segment_manager();}
static void mremove(path f){remove(f);}

std::set<std::string> candidates(manager& m){
	std::set<std::string> ret;
	for(auto it = get_allocator(m)->named_begin(); it != get_allocator(m)->named_end(); ++it)
		ret.insert(std::string(it->name(), it->name_length()));
	return ret;
}

#include "../array.hpp"

#include<numeric> // iota

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D> 
using marray = multi::array<T, D, mallocator<T>>;

BOOST_AUTO_TEST_CASE(multi_test_bip){

path file = "bip_mapped_file.bin";
mremove(file);
{
	manager m{bip::create_only, file.c_str(), 1 << 25};
	auto&& arr1d = 
		*m.construct<marray<int, 1>>("arr1d")(std::make_tuple(10), 99, get_allocator(m));
	auto&& arr2d = 
		*m.construct<marray<double, 2>>("arr2d")(std::make_tuple(1000, 1000), 0.0, get_allocator(m));
	auto&& arr3d = 
		*m.construct<marray<unsigned, 3>>("arr3d")(std::make_tuple(10, 10, 10), 0u, get_allocator(m));

	arr1d[3] = 33;
	arr2d[4][5] = 45.001;
	std::iota(arr3d[6][7].begin(), arr3d[6][7].end(), 100);

//	m.flush(); // this produces uninitialized access in icpc 19.1 and might not be necessary
}
{
	manager m{bip::open_only, file.c_str()};

	auto s = candidates(m);
	BOOST_REQUIRE( s.find("arr1d") != s.end() );
	BOOST_REQUIRE( s.find("arr2d") != s.end() );
	BOOST_REQUIRE( s.find("arr3d") != s.end() );

	auto&& arr1d = *m.find<marray<int, 1>>("arr1d").first; 
	BOOST_REQUIRE(std::addressof(arr1d));
	
	auto&& arr2d = *m.find<marray<double, 2>>("arr2d").first; 
	BOOST_REQUIRE(std::addressof(arr2d));
	
	auto&& arr3d = *m.find<marray<unsigned, 3>>("arr3d").first; 
	BOOST_REQUIRE(std::addressof(arr3d));

	BOOST_REQUIRE( arr1d[5] == 99 );
	BOOST_REQUIRE( arr1d[3] == 33 );

	BOOST_REQUIRE( arr2d[7][8] == 0. );
	BOOST_REQUIRE( arr2d[4][5] == 45.001 );

	BOOST_REQUIRE( arr3d[6][7][3] == 103 );

	m.destroy<marray<int, 1>>("arr1d");
	m.destroy<marray<double, 2>>("arr2d");
	m.destroy<marray<unsigned, 3>>("arr3d");
}
mremove(file);

}

