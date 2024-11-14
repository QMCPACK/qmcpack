#ifdef COMPILATION  // clang-format off
${CXX:-c++} -std=c++17 $CXXFLAGS -I../include $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif  // clang-format on
// Copyright 2019-2024 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi interacting with Boost Interprocess"
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <boost/interprocess/managed_mapped_file.hpp>

#include <boost/interprocess/allocators/adaptive_pool.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>

#include <filesystem>
#include <set>

#include <scoped_allocator>

namespace bip = boost::interprocess;

using manager = bip::managed_mapped_file;

using std::filesystem::path;
using std::filesystem::remove;

template<class T> using mallocator = bip::allocator<T, manager::segment_manager>;
static auto get_allocator(manager& m) { return m.get_segment_manager(); }
static void mremove(path f) { remove(f); }

auto objects_directory(manager& m) {
	std::set<std::string> ret;
	std::transform(
		get_allocator(m)->named_begin(), get_allocator(m)->named_end(),
		std::inserter(ret, ret.end()),
		[](auto const& e) { return e.name(); }
	);
	return ret;
}

#include <multi/array.hpp>

#include <numeric>  // iota

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D>
using marray = multi::array<T, D, mallocator<T>>;

BOOST_AUTO_TEST_CASE(const multi_test_bip) {

	path const file = "bip_mapped_file.bin";
	{
		mremove(file);
		manager m{bip::create_only, file.c_str(), 1 << 25};  // objects with same name produce boost::interprocess_exception::library_error
		auto&&  arr1d = *m.construct<marray<int, 1>>("arr1d")(std::tuple{10}, 99, get_allocator(m));
		auto&&  arr2d = *m.construct<marray<double, 2>>("arr2d")(std::tuple{10, 10}, 0.0, get_allocator(m));
		auto&&  arr3d = *m.construct<marray<unsigned, 3>>("arr3d")(std::tuple{10, 10, 10}, 0u, get_allocator(m));

		arr1d[3]    = 33;
		arr2d[4][5] = 45.001;

		std::iota(arr3d[6][7].begin(), arr3d[6][7].end(), 100);

		auto const& arr3d_copy = *m.construct<marray<unsigned, 3>>("arr3d_copy")(arr3d, get_allocator(m));
		BOOST_REQUIRE( arr3d == arr3d_copy );

		//  m.flush(); // this produces uninitialized access in icpc 19.1 and might not be necessary
	}
	{
		manager m{bip::open_only, file.c_str()};

		auto const s = objects_directory(m);
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

		BOOST_REQUIRE( arr2d[7][8] ==  0.0    );
		BOOST_REQUIRE( arr2d[4][5] == 45.001 );

		BOOST_REQUIRE( arr3d[6][7][3] == 103 );

		auto&& arr3d_copy = *m.find<marray<unsigned, 3>>("arr3d_copy").first;
		BOOST_REQUIRE(std::addressof(arr3d_copy));

		BOOST_REQUIRE( arr3d == arr3d_copy );

		m.destroy<marray<int, 1>>("arr1d");
		m.destroy<marray<double, 2>>("arr2d");
		m.destroy<marray<unsigned, 3>>("arr3d");
		mremove(file);
	}
}

template<class T> using alloc = bip::adaptive_pool<
	T, bip::managed_shared_memory::segment_manager>;

BOOST_AUTO_TEST_CASE(const scoped_allocator_vector_of_arrays) {

	using bipc_row    = multi::array<int, 1, alloc<int>>;
	using bipc_matrix = std::vector<bipc_row, std::scoped_allocator_adaptor<alloc<bipc_row>>>;

	bip::shared_memory_object::remove("Demo");
	{
		bip::managed_shared_memory s{bip::create_only, "Demo", 65536};

		bipc_matrix v(s.get_segment_manager());

		v.emplace_back(multi::extensions_t<1>(3), 99.);
		std::iota(v[0].begin(), v[0].end(), 42);

		assert(v[0][1] == 43);
		bip::shared_memory_object::remove("Demo");
	}
}

BOOST_AUTO_TEST_CASE(const scoped_allocator_arrays_of_vector) {

	using bipc_row    = std::vector<int, alloc<int>>;
	using bipc_matrix = multi::array<bipc_row, 1, std::scoped_allocator_adaptor<alloc<bipc_row>>>;

	bip::shared_memory_object::remove("Demo");
	{
		bip::managed_shared_memory s{bip::create_only, "Demo", 65536};

		bipc_matrix      v(bipc_matrix::extensions_type(10), bipc_row{s.get_segment_manager()}, s.get_segment_manager());
		std::vector<int> row(3, 99);
		v[0].assign(row.begin(), row.end());

		BOOST_REQUIRE( v[0][1] == 99 );
		bip::shared_memory_object::remove("Demo");
	}
}

BOOST_AUTO_TEST_CASE(const scoped_allocator_arrays_of_array) {

	using bipc_row    = multi::array<int, 1, alloc<int>>;
	using bipc_matrix = multi::array<bipc_row, 1, std::scoped_allocator_adaptor<alloc<bipc_row>>>;

	bip::shared_memory_object::remove("Demo");
	{
		bip::managed_shared_memory s{bip::create_only, "Demo", 165536};

		bipc_matrix          v(bipc_matrix::extensions_type(10), bipc_row{bipc_matrix::extensions_type(3), 5, s.get_segment_manager()}, s.get_segment_manager());
		multi::array<int, 1> row = {97, 98, 99};
		std::copy(row.begin(), row.end(), v[0].begin());

		BOOST_REQUIRE( v[0][1] == 98 );
		BOOST_REQUIRE( v[1][1] ==  5 );

		v.reextent(bipc_matrix::extensions_type(12), bipc_row{bipc_matrix::extensions_type(3), 5, s.get_segment_manager()});

		bip::shared_memory_object::remove("Demo");
	}
}
