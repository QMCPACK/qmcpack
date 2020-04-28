#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lpthread -lstdc++fs&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#include<cassert>
#include<numeric> // iota
#include<iostream>
#include<experimental/filesystem> // remove_all c++17

#include <boost/interprocess/managed_mapped_file.hpp>

#if defined(__clang__) && defined(__CUDA__)
namespace std{
	template<class T, class Long, class ULong, int N>
	struct pointer_traits<boost::interprocess::offset_ptr<T, Long, ULong, N>>{
		using pointer = boost::interprocess::offset_ptr<T, Long, ULong, N>;
		using element_type = T;
		using difference_type = std::ptrdiff_t;
		template<class U> using rebind = typename boost::interprocess::offset_ptr<U, Long, ULong, N>;
	//	static pointer pointer_to(element_type& e) noexcept{return std::addressof(e);}
	};
}
#endif

namespace bip = boost::interprocess;

using manager = bip::managed_mapped_file;

template<class T> using mallocator = bip::allocator<T, manager::segment_manager>;
static auto get_allocator(manager& m){return m.get_segment_manager();}
static void mremove(std::experimental::filesystem::path f){std::experimental::filesystem::remove(f);}

//static std::string candidates(manager& m){
//	std::string ret = "  candidates are:\n";
//	for(auto it = get_allocator(m)->named_begin(); it != get_allocator(m)->named_end(); ++it)
//		ret += '\t'+std::string(it->name(), it->name_length()) +'\n';
//	return ret;
//}

#include "../tests/../array.hpp"

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D> 
using marray = multi::array<T, D, mallocator<T>>;

using std::tuple;

int main(){

std::experimental::filesystem::path file = "bip_mapped_file.bin";
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

	m.flush();
}
{
	manager m{bip::open_only, file.c_str()};

	auto&& arr1d = 
		*m.find<marray<int, 1>>("arr1d").first; assert(std::addressof(arr1d));
	auto&& arr2d = 
		*m.find<marray<double, 2>>("arr2d").first; assert(std::addressof(arr2d));
	auto&& arr3d = 
		*m.find<marray<unsigned, 3>>("arr3d").first; assert(std::addressof(arr3d));

	assert( arr1d[5] == 99 );
	assert( arr1d[3] == 33 );

	assert( arr2d[7][8] == 0. );
	assert( arr2d[4][5] == 45.001 );

	assert( arr3d[6][7][3] == 103 );

	m.destroy<marray<int, 1>>("arr1d");//	eliminate<marray<int, 1>>(m, "arr1d"); 
	m.destroy<marray<double, 2>>("arr2d");//	eliminate<marray<double, 2>>(m, "arr2d");
	m.destroy<marray<unsigned, 3>>("arr3d");//	eliminate<marray<unsigned, 3>>(m, "arr3d");
}
mremove(file);
}

