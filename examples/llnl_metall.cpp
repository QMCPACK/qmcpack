#ifdef COMPILATION  // clang-format off
${CXX:-c++} -std=c++17 $CXXFLAGS -I../include -I$HOME/metall/include $0 -o$0.$X&&$0.$X&&rm $0.$X;exit
#endif // clang-format on
// Copyright 2019-2023 Alfredo A. Correa

#include <cassert>
#include <iostream>
#include <numeric>  // for std::iota

#include <metall/metall.hpp>

#include <multi/array.hpp>

template<class T>
using mallocator = metall::manager::allocator_type<T>;

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D>
using marray = multi::array<T, D, mallocator<T>>;

using std::tuple;

int main() {
	std::filesystem::path const dir = "llnl_metall_mapped_file.bin/";
	remove_all(dir);
	{
		metall::manager m{metall::create_only, dir.c_str(), 1 << 25};

		auto&& arr1d = *m.construct<marray<int, 1>>("arr1d")(std::tuple{10}, 5, m.get_allocator());
		auto&& arr2d = *m.construct<marray<double, 2>>("arr2d")(std::tuple{10, 10}, 1.0, m.get_allocator());
		auto&& arr3d = *m.construct<marray<unsigned, 3>>("arr3d")(std::tuple{10, 10, 10}, 1u, m.get_allocator());

		auto&& arr3d_copy = *m.construct<marray<unsigned, 3>>("arr3d_copy")(arr3d, m.get_allocator());

		assert(arr1d[3] == 5);
		assert(arr2d[4][5] == 1.0);
		assert(arr3d[2][3][4] == 1u);

		arr1d[3]    = 33;
		arr2d[4][5] = 45.001;
		std::iota(arr3d[6][7].begin(), arr3d[6][7].end(), 100);

		assert(arr3d_copy[6][7][8] == 1u);

		auto&& arr3d_assign = *m.construct<marray<unsigned, 3>>("arr3d_assign")(m.get_allocator());
		arr3d_assign        = arr3d;

		assert(arr3d_assign == arr3d);

		assert(arr3d_assign[6][7][8] == arr3d[6][7][8]);
		//  m.flush();
	}
	{
		metall::manager m{metall::open_only, dir.c_str()};

		auto const& arr1d = *m.find<marray<int, 1>>("arr1d").first;
		auto const& arr2d = *m.find<marray<double, 2>>("arr2d").first;
		auto const& arr3d = *m.find<marray<unsigned, 3>>("arr3d").first;

		auto const& arr3d_copy = *m.find<marray<unsigned, 3>>("arr3d_copy").first;
		assert(std::addressof(arr3d));

		auto const& arr3d_assign = *m.find<marray<unsigned, 3>>("arr3d_assign").first;
		assert(std::addressof(arr3d));

		assert(arr1d[5] == 5);
		assert(arr1d[3] == 33);

		assert(arr2d[7][8] == 1.0);
		assert(arr2d[4][5] == 45.001);

		assert(arr3d[6][7][3] == 103);

		assert(arr3d_assign == arr3d);

		m.destroy<marray<int, 1>>("arr1d");
		m.destroy<marray<double, 2>>("arr2d");
		m.destroy<marray<unsigned, 3>>("arr3d");
		m.destroy<marray<unsigned, 3>>("arr3d_copy");
		m.destroy<marray<unsigned, 3>>("arr3d_assign");
	}

	remove_all(dir);
}