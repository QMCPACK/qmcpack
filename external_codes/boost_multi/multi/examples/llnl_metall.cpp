#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX -std=c++17 -I~/https/github.com/LLNL/metall.git/include/ $0 -o $0x -lstdc++fs&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#include<cassert>
#include<numeric> // iota
#include<iostream>

#include<metall/metall.hpp>

#include "../../multi/array.hpp"

template<class T> using mallocator = metall::manager::allocator_type<T>;

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D> 
using marray = multi::array<T, D, mallocator<T>>;

using std::tuple;

int main(){
std::filesystem::path dir = "llnl_metall_mapped_file.bin/";
remove_all(dir);
{
	metall::manager m{metall::create_only, dir.c_str(), 1<<25};
	auto&& arr1d = 
		*m.construct<marray<int     , 1>>("arr1d")(tuple{10}        , 99 , m.get_allocator());
	auto&& arr2d = 
		*m.construct<marray<double  , 2>>("arr2d")(tuple{1000, 1000}, 1.0, m.get_allocator());
	auto&& arr3d = 
		*m.construct<marray<unsigned, 3>>("arr3d")(tuple{10, 10, 10}, 1u , m.get_allocator());
	auto&& arr3d_cpy = 
		*m.construct<marray<unsigned, 3>>("arr3d_cpy")(tuple{0, 0, 0}, m.get_allocator());

	assert( arr1d[3] == 99 );
	assert( arr2d[4][5] == 1.0 );
	assert( arr3d[2][3][4] == 1u );

	arr1d[3] = 33;
	arr2d[4][5] = 45.001;
	std::iota(arr3d[6][7].begin(), arr3d[6][7].end(), 100);

	arr3d_cpy = arr3d;
	assert( arr3d_cpy[6][7][8] == arr3d[6][7][8] );
	m.flush();
}
{
	metall::manager m{metall::open_only, dir.c_str()};

	auto&& arr1d =
		*m.find<marray<int     , 1>>("arr1d").first; assert(std::addressof(arr1d));
	auto&& arr2d =
		*m.find<marray<double  , 2>>("arr2d").first; assert(std::addressof(arr2d));
	auto&& arr3d =
		*m.find<marray<unsigned, 3>>("arr3d").first; assert(std::addressof(arr3d));
	auto&& arr3d_cpy =
		*m.find<marray<unsigned, 3>>("arr3d_cpy").first; assert(std::addressof(arr3d));

	assert( arr1d[5] == 99 );
	assert( arr1d[3] == 33 );

	assert( arr2d[7][8] == 1.0 );
	assert( arr2d[4][5] == 45.001 );

	assert( arr3d[6][7][3] == 103 );
	assert( arr3d_cpy == arr3d );

	m.destroy<marray<int     , 1>>("arr1d");
	m.destroy<marray<double  , 2>>("arr2d");
	m.destroy<marray<unsigned, 3>>("arr3d");
	m.destroy<marray<unsigned, 3>>("arr3d_cpy");
}
remove_all(dir);
}

