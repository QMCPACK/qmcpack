// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi allocators"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<fstream>

// #include "../detail/serialization.hpp"

#if 0  // 1 to test cereal
	#include <cereal/cereal.hpp>

	#include <cereal/archives/binary.hpp>
	#include <cereal/archives/xml.hpp>
	#include <cereal/types/string.hpp>

	using XOArchive = cereal::XMLOutputArchive;
	using XIArchive = cereal::XMLInputArchive;

	using BOArchive = cereal::BinaryOutputArchive;
	using BIArchive = cereal::BinaryInputArchive;

	using cereal::make_nvp;
#else
	#include <boost/archive/binary_iarchive.hpp>
	#include <boost/archive/binary_oarchive.hpp>

	#include <boost/archive/xml_iarchive.hpp>
	#include <boost/archive/xml_oarchive.hpp>
	#include <boost/serialization/string.hpp>

	using XOArchive = boost::archive::xml_oarchive;
	using XIArchive = boost::archive::xml_iarchive;

	using BOArchive = boost::archive::binary_oarchive;
	using BIArchive = boost::archive::binary_iarchive;

	using boost::serialization::make_nvp;
//	using boost::serialization::make_array;
#endif

#include <numeric>
#include <string>

#include <boost/multi_array.hpp>

namespace multi = boost::multi;

struct array {
	using input_archive  = boost::archive::xml_iarchive;  // cereal::JSONInputArchive ;
	using output_archive = boost::archive::xml_oarchive;  // cereal::JSONOutputArchive;

	template<class Array, class IStream>
	static auto load(IStream&& is) -> Array {
		using boost::serialization::make_nvp;  // cereal::make_nvp;  //
		Array value{};
		input_archive{is} >> make_nvp("value", value);
		return value;
	}

	template<class Array, class OStream>
	static void save(OStream&& os, Array const& value) {
		using boost::serialization::make_nvp;  // cereal::make_nvp;  //
		output_archive{os} << make_nvp("value", value);
	}
};

BOOST_AUTO_TEST_CASE(json) {
	namespace multi = boost::multi;
	multi::array<std::string, 2> A = {{"00", "01"}, {"10", "11"}};
	array::save(std::ofstream{"file"}, A);

	auto B = array::load<multi::array<std::string, 2>>(std::ifstream{"file"});
	BOOST_REQUIRE(A == B);
}

BOOST_AUTO_TEST_CASE(extensions_serialization) {
	multi::array<double, 2> arr({10, 10});
	auto const x = arr.extensions();
	std::stringstream ss;
	{
		XOArchive xoa{ss};
		xoa<<                                   make_nvp("x", x);
	//	xoa<< multi::archive_traits<XOArchive>::make_nvp("x", x);
	//	xoa<<                                          AR_NVP(x);
	//	xoa<<                                      CEREAL_NVP(x);
	//	xoa<<                                                 x ;
	}
	{
		multi::extensions_t<2> y;
		{
			XIArchive xia{ss};
			xia>>                                   make_nvp("x", y);
		//	xia>> multi::archive_traits<XIArchive>::make_nvp("x", y);
		//	xia>>                           cereal::make_nvp("x", y);
		//	xia>>                                                 y ;
		}
		BOOST_REQUIRE(x == y);
	}
}

BOOST_AUTO_TEST_CASE(carray_serialization) {
	double const A[3][3] = {{0., 1., 2.}, {3., 4., 5.}, {6., 7., 8.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
	std::stringstream ss;
	{
		{
			XOArchive xoa{ss};
			xoa<<                                   make_nvp("A", A);
		//	xoa<<                                          AR_NVP(A);
		//	xoa<<                                      CEREAL_NVP(A);
		//	xoa<<                                                 A ;
		//	xoa<< multi::archive_traits<XOArchive>::make_nvp("A", A);
		}
		std::ofstream ofs{"serialization_A.xml"};
		ofs<< ss.str();
	}
	{
		double B[3][3];  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
		XIArchive xia{ss};
		xia>>                                   make_nvp("A", B);
	//	xia>>                           cereal::make_nvp("A", B);
	//	xia>>                                                 B ;
	//	xia>> multi::archive_traits<XIArchive>::make_nvp("A", B);
		BOOST_REQUIRE( B[1][2] == 5. );  // NOLINT(clang-analyzer-core.UndefinedBinaryOperatorResult) is it?
		BOOST_REQUIRE( A[1][2] == B[1][2] );  // NOLINT(clang-analyzer-core.UndefinedBinaryOperatorResult) is it?
	}
}

BOOST_AUTO_TEST_CASE(array_serialization) {
	multi::array<double, 2> arr({10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		XOArchive xoa{ss};
		xoa<<      make_nvp("arr", arr);
	//	xoa<<               AR_NVP(arr);
	//	xoa<<           CEREAL_NVP(arr);
	//	xoa<<                      arr ;
	//	xoa<< CEREAL_NVP(arr);
	//	xoa<< multi::archive_traits<XOArchive>::make_nvp("arr", arr);
	}
	{
		multi::array<double, 2> arr2;
		{
			XIArchive xia{ss};
			xia>>                                   make_nvp("arr", arr2);
		//	xia>>                           cereal::make_nvp("arr", arr2);
		//	xia>>                                                   arr2 ;
		//	xia>> multi::archive_traits<XIArchive>::make_nvp("arr", arr2);
		}
		BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
		BOOST_REQUIRE( arr2 == arr );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_string) {
	multi::array<std::string, 2> arr({10, 10});
	auto const x = extensions(arr);
	for(auto i : std::get<0>(x) ) {
		for(auto j : std::get<1>(x) ) {
			arr[i][j] = std::to_string(i) + std::to_string(j);
		}
	}

	std::stringstream ss{};
	{
		XOArchive xoa{ss};
		xoa<<                                   make_nvp("arr", arr);
	//	xoa<<                                            AR_NVP(arr) ;
	//	xoa<<                           BOOST_SERIALIZATION_NVP(arr) ;
	//	xoa<<                                        CEREAL_NVP(arr) ;
	//	xoa<<                                                   arr ;
	//	xoa<< multi::archive_traits<XOArchive>::make_nvp("arr", arr);
	}
	{
		multi::array<std::string, 2> arr2{};
		{
			XIArchive xia{ss};
			xia>>                                   make_nvp("arr", arr2);
		//	xia>>                           cereal::make_nvp("arr", arr2);
		//	xia>>                                                   arr2 ;
		//	xia>> multi::archive_traits<XIArchive>::make_nvp("arr", arr2);
		}
		BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
		BOOST_REQUIRE( arr2 == arr );
	}
}

//#if not defined(__NVCC__)  // some code contained here doesn't compile with nvcc 11.0,11.1 and 11.2
BOOST_AUTO_TEST_CASE(array_serialization_binary) {
	multi::array<double, 2> arr({10, 10}, 0.);
	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		BOArchive boa(ss);
		boa<< arr;
	}
	{
		multi::array<double, 2> arr2{};
		{
			BIArchive bia{ss};
			bia>> arr2;
		}
		BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
		BOOST_REQUIRE( arr2 == arr );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_string_binary) {
	multi::array<std::string, 2> arr({10, 10});
	auto const x = extensions(arr);
	for(auto i : std::get<0>(x) ) {
		for(auto j : std::get<1>(x) ) {
			arr[i][j] = std::to_string(i) + std::to_string(j);
		}
	}

	std::stringstream ss{};
	{
		BOArchive boa{ss};
		boa<< arr;
	//	boa<< multi::archive_traits<BOArchive>::make_nvp("arr", arr);
	}
	{
		multi::array<std::string, 2> arr2{};
		{
			BIArchive bia{ss};
			bia>> arr2;
			//	bia>> multi::archive_traits<BIArchive>::make_nvp("arr", arr2);
		}
		BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
		BOOST_REQUIRE( arr2 == arr );
	}
}

//#if not defined(__NVCC__)  // some code contained here doesn't compile with nvcc 11.0,11.1 and 11.2
BOOST_AUTO_TEST_CASE(vector) {
	std::vector<double> v(100); std::iota(begin(v), end(v), 10.);

	std::stringstream ss;
	{
		XOArchive xoa{ss};
		xoa<< make_nvp("v_data", multi::archive_traits<XOArchive>::make_array(v.data(), v.size()));
	//	xoa<< make_nvp("v_data",                                   make_array(v.data(), v.size()));
	//	xoa<< make_nvp("v_data",             boost::serialization::make_array(v.data(), v.size()));
	}
	{
		std::vector<double> w(100);
		XIArchive xia{ss};
		xia>> make_nvp("v_data", multi::archive_traits<XIArchive>::make_array(w.data(), w.size()));
	//	xia>> make_nvp("v_data",                                   make_array(w.data(), w.size()));
	//	xia>> make_nvp("v_data",             boost::serialization::make_array(w.data(), w.size()));
		BOOST_REQUIRE( v == w );
	}
}

BOOST_AUTO_TEST_CASE(vector_binary) {
	std::vector<double> v(100); std::iota(begin(v), end(v), 10.);

	std::stringstream ss{};
	{
		BOArchive xoa{ss};
		xoa<< make_nvp("v_data", multi::archive_traits<XOArchive>::make_array(v.data(), v.size()));
	//	xoa<< make_nvp("v_data",                                   make_array(v.data(), v.size()));
	//	xoa<< make_nvp("v_data",             boost::serialization::make_array(v.data(), v.size()));
	}
	{
		std::vector<double> w(100);
		BIArchive xia{ss};
		xia>> make_nvp("v_data", multi::archive_traits<XIArchive>::make_array(w.data(), w.size()));
	//	xia>> make_nvp("v_data",                                   make_array(w.data(), w.size()));
	//	xia>> make_nvp("v_data",             boost::serialization::make_array(w.data(), w.size()));
		BOOST_REQUIRE( v == w );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_3D) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		XOArchive xoa{ss};
		xoa<<      make_nvp("arr", arr);
	//	xoa<<               AR_NVP(arr);
	//	xoa<<           CEREAL_NVP(arr);
	//	xoa<<                      arr ;
	//	xoa<< multi::archive_traits<XOArchive>::make_nvp("arr", arr);
	}
	{
		multi::array<double, 3> arr2{};
		{
			XIArchive xia{ss};
			xia>>                                   make_nvp("arr", arr2);
		//	xia>>                           cereal::make_nvp("arr", arr2);
		//	xia>>                                                   arr2 ;
		//	xia>> multi::archive_traits<XIArchive>::make_nvp("arr", arr2);
		}
		BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
		BOOST_REQUIRE( arr2 == arr );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_3D_inplace) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	XOArchive{ss}<< make_nvp("arr", arr);

	multi::array<double, 3> arr2{};
	XIArchive{ss}>> make_nvp("arr", arr2);

	BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
	BOOST_REQUIRE( arr2 == arr );
}

BOOST_AUTO_TEST_CASE(array_serialization_2D_inplace_file) {
	multi::array<double, 2> arr({2, 2}, 99.);

	{
		std::ofstream ofs{"file.xml"};
		XOArchive{ofs}<< make_nvp("arr", arr);
	}  // flush the file stream

	multi::array<double, 2> arr2{};
	std::ifstream ifs{"file.xml"};
	XIArchive{ifs}>> make_nvp("arr", arr2);

	BOOST_REQUIRE( extensions(arr2) == extensions(arr) );
	BOOST_REQUIRE( arr2 == arr );
}

#if not defined(__NVCC__)  // some code contained here doesn't compile with nvcc 11.0,11.1 and 11.2
BOOST_AUTO_TEST_CASE(array_serialization_3D_part_binary_lvalue) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		BOArchive boa{ss};
		auto&& arr2 = arr[2];
		boa& arr2;
	}
	{
		BOOST_REQUIRE( arr[3] != arr[2] );
		{
			BIArchive bia{ss};
			auto&& arr3 = arr[3];
			bia& arr3;
		}
		BOOST_REQUIRE( arr[3] == arr[2] );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_3D_part_xml_lvalue) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		XOArchive boa{ss};
		auto&& arr2 = arr[2];
		boa<< multi::archive_traits<XOArchive>::make_nvp("arr2", arr2);
	}
	{
		BOOST_REQUIRE( arr[3] != arr[2] );
		{
			XIArchive bia{ss};
			auto&& arr3 = arr[3];
			bia>> multi::archive_traits<XOArchive>::make_nvp("arr2", arr3);
		}
		BOOST_REQUIRE( arr[3] == arr[2] );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_3D_part_binary) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		BOArchive boa{ss};
		boa& arr[2];
	}
	{
		BOOST_REQUIRE( arr[3] != arr[2] );
		{
			BIArchive bia{ss};
			bia& arr[3];
		}
		BOOST_REQUIRE( arr[3] == arr[2] );
	}
}

BOOST_AUTO_TEST_CASE(array_serialization_3D_part_xml) {
	multi::array<double, 3> arr({10, 10, 10}, 0.);

	BOOST_REQUIRE(( arr.extension() == boost::multi::index_range{0, 10} ));

	std::iota(arr.data_elements(), arr.data_elements() + arr.num_elements(), 1000.);

	std::stringstream ss{};
	{
		XOArchive boa{ss};
		boa<< multi::archive_traits<XOArchive>::make_nvp("arr2", arr[2]);
	}
	{
		BOOST_REQUIRE( arr[3] != arr[2] );
		{
			XIArchive bia{ss};
			bia>> multi::archive_traits<XOArchive>::make_nvp("arr2", arr[3]);
		}
		BOOST_REQUIRE( arr[3] == arr[2] );
	}
}
#endif
