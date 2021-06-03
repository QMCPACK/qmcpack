#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework -lstdc++fs -lboost_serialization -lboost_iostreams&&$0x&&rm $0x;exit
#endif
// © Alfredo Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi serialization"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../adaptors/serialization/xml_archive.hpp"
//#include<boost/archive/xml_oarchive.hpp>
//#include<boost/archive/xml_iarchive.hpp>

#include<boost/archive/text_oarchive.hpp>
#include<boost/archive/text_iarchive.hpp>
#include<boost/archive/binary_oarchive.hpp>
#include<boost/archive/binary_iarchive.hpp>
#include<boost/serialization/nvp.hpp>
#include<boost/serialization/binary_object.hpp>
#include<boost/serialization/complex.hpp>

#include "../array.hpp"

//#include "../adaptors/cuda.hpp"

#include<boost/iostreams/filtering_stream.hpp>
#include<boost/iostreams/filter/gzip.hpp>

#include<fstream>
#include<experimental/filesystem>

#include<chrono>
#include<iostream>
#include<random>

#include<boost/serialization/vector.hpp>

#include<boost/multi_array.hpp>

namespace multi = boost::multi;
namespace fs = std::experimental::filesystem;

struct watch : private std::chrono::high_resolution_clock{
	std::string name_;
	time_point  start_;
	mutable bool engaged = true;
	watch(std::string name = "") : name_{name}, start_{now()}{}
	auto operator*() const{engaged = false; return std::chrono::duration<double>(now() - start_).count();}
	~watch(){
		if(engaged){
			auto count = operator*();
			std::cerr<< name_ <<": "<< count <<" sec"<<std::endl;
		}
	}
};

BOOST_AUTO_TEST_CASE(multi_serialization_static_small_xml){
	multi::static_array<double, 2> d2D({10, 10});
	std::mt19937 eng{std::random_device{}()};
	auto gen = [&](){return std::uniform_real_distribution<>{}(eng);};
	std::for_each(begin(d2D), end(d2D), [&](auto&& r){std::generate(begin(r), end(r), gen);});
	auto const name = "serialization-static-small.xml";
	{
		std::ofstream ofs{name}; assert(ofs);
		boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
	}
	{
		std::ifstream ifs{name}; assert(ifs); decltype(d2D) d2D_copy(extensions(d2D), 99.);
		boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
		BOOST_REQUIRE( d2D_copy == d2D );
	}
	std::cout<< fs::file_size(name) <<'\n';
	fs::remove(name);
}

BOOST_AUTO_TEST_CASE(multi_serialization_small_xml){
	multi::array<double, 2> d2D({10, 10});
	std::mt19937 e{std::random_device{}()};
//	auto g = std::bind(std::uniform_real_distribution<>{}, e);//
	auto g = [&](){return std::uniform_real_distribution<>{}(e);};
	std::for_each(begin(d2D), end(d2D), [&](auto&& r){std::generate(begin(r), end(r), g);});
	auto const name = "serialization-small.xml";
	{
		std::ofstream ofs{name}; assert(ofs);
		boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
	}
	{
		std::ifstream ifs{name}; assert(ifs); decltype(d2D) d2D_copy(extensions(d2D));
		boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
		BOOST_REQUIRE( d2D_copy == d2D );
	}
	{
		std::ofstream ofs{"serialization-small-part.xml"}; assert(ofs);
		auto&& a = d2D({0, 5}, {0, 5});
		boost::archive::xml_oarchive{ofs} << boost::serialization::make_nvp("d2D_part", a);//BOOST_SERIALIZATION_NVP(d2D);
		fs::remove("serialization-small-part.xml");
	}
	std::cout<< fs::file_size(name) <<'\n';
	fs::remove(name);
}


BOOST_AUTO_TEST_CASE(multi_serialization_static_large_xml){
	watch w("static_large_xml");
	multi::static_array<double, 2> d2D({1000, 1000});
	auto gen = [e=std::mt19937{std::random_device{}()}]() mutable{return std::uniform_real_distribution<>{}(e);};
	std::for_each(begin(d2D), end(d2D), [&](auto&& r){std::generate(begin(r), end(r), gen);});
	auto const name = "serialization-static-large.xml";
	{
		std::ofstream ofs{name}; assert(ofs);
		boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
	}
	{
		std::ifstream ifs{name}; assert(ifs); decltype(d2D) d2D_copy(extensions(d2D));
		boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
		BOOST_REQUIRE( d2D_copy == d2D );
	}
	std::cout<< fs::file_size(name) <<'\n';
	fs::remove(name);
}

BOOST_AUTO_TEST_CASE(multi_serialization_static_small){
	{
		multi::static_array<double, 0> d0D = 12.;
		std::ofstream ofs{"serialization-static_0D.xml"}; assert(ofs);
		boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d0D);
		fs::remove("serialization-static_0D.xml");
	}
	{
		multi::array<double, 2> d2D = {
			{150., 16., 17., 18., 19.}, 
			{  5.,  5.,  5.,  5.,  5.}, 
			{100., 11., 12., 13., 14.}, 
			{ 50.,  6.,  7.,  8.,  9.}  
		};
		auto gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable{return d(e);};
		std::for_each(
			begin(d2D), end(d2D), 
			[&](auto&& r){std::generate(begin(r), end(r), gen);}
		);
		auto name = "serialization-small-double2D.xml";
		[&, _ = watch("xml write double")]{
			std::ofstream ofs{"serialization-small-double2D.xml"}; assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		}();
		std::cerr<<"size "<< double(fs::file_size(name))/1e6 <<"MB\n";
		fs::remove("serialization-small-double2D.xml");
	}
	{
		multi::array<double, 2> d2D = {
			{150., 16., 17., 18., 19.}, 
			{  5.,  5.,  5.,  5.,  5.}, 
			{100., 11., 12., 13., 14.}, 
			{ 50.,  6.,  7.,  8.,  9.}  
		};
		d2D.reextent({2000, 2000});
		auto gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable{return d(e);};
		std::for_each(
			begin(d2D), end(d2D), 
			[&](auto&& r){std::generate(begin(r), end(r), gen);}
		);
		[&, _ = watch("xml write double")]{
			std::ofstream ofs{"serialization-double.xml"}; assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		}();
		std::cerr<<"size "<< double(fs::file_size("serialization-double.xml"))/1e6 <<"MB\n";
		fs::remove("serialization-double.xml");
	}

	using complex = std::complex<float>;

	auto const d2D = []{multi::array<complex, 2> _({10000, 1000});
		auto gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable{return std::complex<double>{d(e), d(e)};};
		std::for_each(begin(_), end(_), [&](auto&& r){std::generate(begin(r), end(r), gen);}); return _;
	}();
	auto size = sizeof(double)*d2D.num_elements();
	using std::cerr;
	std::cout<<"data size (in memory) "<< size <<std::endl;
	{
		fs::path file{"serialization.bin"};
		auto count = [&, w=watch("binary write")]{
			std::ofstream ofs{file}; assert(ofs);
			boost::archive::binary_oarchive{ofs} << d2D;
			return *w;
		}();
		cerr<<"size  "<< double(file_size(file))/1e6 <<"MB\n";
		cerr<<"speed "<< double(size)/1e6/count <<"MB/s\n";
		std::decay_t<decltype(d2D)> d2D_cpy;
		auto count_load = [&, w=watch("binary load")]{
			std::ifstream ifs{file}; assert(ifs);
			boost::archive::binary_iarchive{ifs} >> d2D_cpy;
			return *w;
		}();
		std::cerr<<"load speed "<< double(file_size(file))/1e6/count_load <<"MB/s\n";
		BOOST_REQUIRE( d2D == d2D_cpy );
		fs::remove(file);
	}
	{
		using std::cout;
		fs::path file{"serialization.xml"};
		cout<< file << std::endl;
		auto count = [&, w = watch("xml write base64")]{
			std::ofstream ofs{file}; assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
			return *w;
		}();
		cout<<"data size "<< double(size)/1e6 << "MB\n";
		cout<<"file size "<< double(file_size(file))/1e6 <<"MB\n";
		cout<<"save speed "<< double(size)/1e6/count <<"MB/s"<< std::endl;
		multi::array<complex, 2> d2D_cpy;
		auto count2 = [&, w = watch("xml load base64")]{
			std::ifstream ifs{file}; assert(ifs);
			boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_cpy);
			return *w;
		}();
		cout<<"load speed "<< double(size)/1e6/count2 <<"MB/s"<< std::endl;
		BOOST_REQUIRE( d2D_cpy == d2D );
		fs::remove(file);
	}
	return;
#if 0
	{
		multi::cuda::managed::array<complex, 2> cud2D({2000, 2000});
		[&, _=watch("cuda binary write")]{
			std::ofstream ofs{"serialization.bin"}; assert(ofs);
			boost::archive::binary_oarchive{ofs} << cud2D;
		}();
		std::cerr<<"size "<< (fs::file_size("serialization.bin")/1e6) <<"MB\n";
	}
#endif
	{
		[&, _ = watch("text write")]{
			std::ofstream ofs{"serialization.txt"}; assert(ofs);
			boost::archive::text_oarchive{ofs} << d2D;
		}();
		std::cerr<<"size "<< double(fs::file_size("serialization.txt"))/1e6 <<"MB\n";
		fs::remove("serialization.txt");
	}
	{
		multi::array<complex, 2> d2D_copy;//(extensions(d2D), 9999.);
		[&, _ = watch("text read")]{
			std::ifstream ifs{"serialization.txt"}; assert(ifs);
			boost::archive::text_iarchive{ifs} >> d2D_copy;
		}();
		BOOST_REQUIRE( d2D_copy == d2D );
		fs::remove("serialization.txt");
	}
	{
		multi::array<complex, 2> d2D_copy;//(extensions(d2D), 9999.);
		[&, _=watch("binary read")]{
			std::ifstream ifs{"serialization.bin"}; assert(ifs);
			boost::archive::binary_iarchive{ifs} >> d2D_copy;
		}();
		BOOST_REQUIRE( d2D_copy == d2D );
		fs::remove("serialization.bin");
	}
	{
		[&, _=watch("binary compressed write")]{
			std::ofstream ofs{"serialization_compressed.bin.gz"};
			{
				boost::iostreams::filtering_stream<boost::iostreams::output> f;
				f.push(boost::iostreams::gzip_compressor());
				f.push(ofs);
				boost::archive::binary_oarchive{f} << d2D;
			}
		}();
		std::cerr<<"size "<< double(fs::file_size("serialization.bin.gz"))/1e6 <<"MB\n";
		fs::remove("serialization.bin.gz");
	}
	{
		[&, _ = watch("compressed xml write")]{
			std::ofstream ofs{"serialization.xml.gz"}; assert(ofs);
			{
				boost::iostreams::filtering_stream<boost::iostreams::output> f;
				f.push(boost::iostreams::gzip_compressor());
				f.push(ofs);
				boost::archive::xml_oarchive{f} << BOOST_SERIALIZATION_NVP(d2D);
			}
		}();
		std::cerr<<"size "<< double(fs::file_size("serialization.xml.gz"))/1e6 <<"MB\n";
		fs::remove("serialization.xml.gz");
	}
	{
		multi::array<complex, 2> d2D_copy;//(extensions(d2D), 9999.);
		[&, _ = watch("xml read")]{
			std::ifstream ifs{"serialization.xml"}; assert(ifs);
			boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
		}();
		BOOST_REQUIRE( d2D_copy == d2D );
		fs::remove("serialization.xml");
	}
}

BOOST_AUTO_TEST_CASE(test_utility_serialization_2d){
	double carr[3][10] = {
		{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
		{10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
		{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
	};
	multi::array_ref<double, 2> marr(&carr[0][0], {3, 10});
	boost::multi_array_ref<double, 2> Marr(&carr[0][0], boost::extents[3][10]);

	namespace arxiv = boost::archive;
	{
		std::ofstream ofs{"utility_serialization_marr.xml"}; assert(ofs);
		arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(marr);
		fs::remove("utility_serialization_marr.xml");
	}
	{
		std::ofstream ofs{"utility_serialization_marr_as_value.xml"}; assert(ofs);
		multi::array<double, 2> const& marr_value = decay(marr);//static_cast<multi::array<double, 2> const&>(marr);
		BOOST_REQUIRE( marr_value.data_elements() == marr.data_elements() );
		arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(marr_value);
		fs::remove("utility_serialization_marr_as_value.xml");
	}
	{
		std::ofstream ofs{"utility_serialization_carr.xml"}; assert(ofs);
		arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(carr);
		fs::remove("utility_serialization_carr.xml");
	}
	{
//		std::ofstream ofs{"utility_serialization_Marr.xml"}; assert(ofs);
//		arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(Marr);
	}
}

