#ifdef COMPILATION_INSTRUCTIONS
c++ -std=c++17 `#-Wfatal-errors` -I$HOME/soft/metall/include $0 -o $0.x -lpthread -lrt -lstdc++fs && $0.x  $@ && rm -f $0.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2019
#include<cassert>
#include<numeric> // iota
#include<iostream>
#if 0 // experimental use of "metall" allocator
#include<metall/metall.hpp>
using namespace metall;
template<class T> using mallocator = manager::allocator_type<T>;
decltype(auto) get_allocator(manager& m){return m.get_allocator();}
void sync(manager& m){m.sync();}
void mremove(std::string const& s){for(auto fix:{"bin_directory","chunk_directory","named_object_directory","segment"}) std::filesystem::remove(s +"_"+ fix);}
std::string candidates(manager& m){return "";}
#endif
#if 1 // boost interprocess declarations
#include <boost/interprocess/managed_mapped_file.hpp>
using namespace boost::interprocess;
using manager = managed_mapped_file;
template<class T> using mallocator = allocator<T, manager::segment_manager>;
auto get_allocator(manager& m){return m.get_segment_manager();}
void sync(manager& m){m.flush();}
void mremove(char const* f){std::remove(f);}//shared_memory_object::remove(f);}
std::string candidates(manager& m){
	std::string ret = "  candidates are:\n";
	for(auto it = get_allocator(m)->named_begin(); it != get_allocator(m)->named_end(); ++it)
        ret += '\t'+std::string(it->name(), it->name_length()) +'\n';
    return ret;
}
#endif

#include "../../multi/array.hpp"

namespace multi = boost::multi;
template<class T, auto D> using marray = multi::array<T, D, mallocator<T>>;

#if 0 // convenience functions for boost interprocess
template<class T, class... Args>
auto create(manager& m, std::string const& s, Args&&... args, void* = 0)
->decltype(*m.construct<T>(s.c_str())(std::forward<Args>(args)...)){
	return *m.construct<T>((s +"."+ typeid(T).name()).c_str())(std::forward<Args>(args)...);}
template<class T, class... Args>
auto create(manager& m, std::string const& s, Args&&... args)
->decltype(*m.construct<T>(s.c_str())(std::forward<Args>(args)..., get_allocator(m))){
	return *m.construct<T>((s +"."+ typeid(T).name()).c_str())(std::forward<Args>(args)..., get_allocator(m));}

template<class T>
decltype(auto) emerge(manager& m, std::string const& name){
	auto p= m.find<T>((name +"."+ typeid(T).name()).c_str()).first;
	if(!p) throw std::runtime_error{"object named\n\t"+ name +"."+ typeid(T).name() +"\n  not found,\n  check name and <type> consistency\n" + candidates(m)};
	return *p;
}
template<class T> 
void eliminate(manager& m, std::string const& s){m.destroy<T>((s +"."+ typeid(T).name()).c_str());}
#endif

using std::tuple;

int main(){
mremove("mapped_file.bin");
{
	manager m{create_only, "mapped_file.bin", 1 << 25};
	auto&& arr1d = //create<marray<int, 1>>(m, "arr1d", tuple{10}, 99); 
		*m.construct<marray<int, 1>>("arr1d")(tuple{10}, 99, get_allocator(m));
	auto&& arr2d = //create<marray<double, 2>>(m, "arr2d", tuple{1000, 1000}, 0.); 
		*m.construct<marray<double, 2>>("arr2d")(tuple{1000, 1000}, 0.0, get_allocator(m));
	auto&& arr3d = //create<marray<unsigned, 3>>(m, "arr3d", tuple{10, 10, 10}, 0); 
		*m.construct<marray<unsigned, 3>>("arr3d")(tuple{10, 10, 10}, 0u, get_allocator(m));

	arr1d[3] = 33;
	arr2d[4][5] = 45.001;
	std::iota(arr3d[6][7].begin(), arr3d[6][7].end(), 100);

	sync(m);
}

{
	manager m{open_only, "mapped_file.bin"};

	auto&& arr1d = //emerge<marray<int, 1>>(m, "arr1d");
		*m.find<marray<int, 1>>("arr1d").first; assert(std::addressof(arr1d));
	auto&& arr2d = //emerge<marray<double, 2>>(m, "arr2d");
		*m.find<marray<double, 2>>("arr2d").first; assert(std::addressof(arr2d));
	auto&& arr3d = //emerge<marray<unsigned, 3>>(m, "arr3d");
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
mremove("mapped_file.bin");
}

