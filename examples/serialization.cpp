#if 0 && defined(COMPILATION)
${CXX:-c++} $0 -o $0x -I../include -lstdc++fs -lboost_serialization -lboost_iostreams&& $0x&& rm $0x;
exit
#endif
// Copyright 2018-2026 Alfredo A. Correa

#include <multi/array.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/multi_array.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

namespace multi = boost::multi;

namespace fs = std::filesystem;

struct watch : private std::chrono::high_resolution_clock {
	std::string  name_;
	time_point   start_  = std::chrono::high_resolution_clock::now();
	mutable bool engaged = true;

	watch() = default;
	watch(std::string_view name) : name_{name} {}

	auto operator*() const {
		engaged = false;
		return std::chrono::duration<double>(now() - start_).count();
	}
	auto operator=(watch const&) = delete;
	~watch() {
		if(engaged) {
			auto count = operator*();
			std::cerr << name_ << ": " << count << " sec" << std::endl;
		}
	}
};

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)

	BOOST_AUTO_TEST_CASE(print_xml) {
		multi::array<std::string, 2> const A{
			{"w", "x"},
			{"y", "z"},
		};
		boost::archive::xml_oarchive(std::cout, boost::archive::no_header)
			<< boost::make_nvp("A", A());
	}

	BOOST_AUTO_TEST_CASE(print_xml) {
		multi::array<std::string, 2> const A{
			{"w", "x"},
			{"y", "z"},
		};
		auto&& AA = A();
		boost::archive::xml_oarchive(std::cout, boost::archive::no_header)
			<< boost::make_nvp("A", AA);
	}

	BOOST_AUTO_TEST_CASE(multi_serialization_static_small_xml) {
		multi::dynamic_array<double, 2> d2D({10, 10});

		std::mt19937_64 eng(std::random_device{}());

		auto gen = [&]() {
			return std::uniform_real_distribution<>{}(eng);
		};

		std::for_each(begin(d2D), end(d2D), [&](auto&& r) { std::generate(begin(r), end(r), gen); });
		std::string const filename = "serialization-static-small.xml";
		{
			std::ofstream ofs{filename};
			assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		}
		{
			std::ifstream ifs{filename};
			assert(ifs);
			decltype(d2D) d2D_copy(extensions(d2D), 99.0);
			boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
			BOOST_TEST( d2D_copy == d2D );
		}
		std::cout << fs::file_size(filename) << '\n';
		fs::remove(filename);

		{
			std::ostringstream oss;
			{
				boost::archive::text_oarchive xoa{oss};

				std::vector<int> v = {1, 2, 3};
				std::for_each(v.begin(), v.end(), [&xoa](auto const& e) { xoa << e; });
				// std::accumulate(v.begin(), v.end(), &xoa, [](boost::archive::text_oarchive* x, int e) {return &(*x << BOOST_SERIALIZATION_NVP(e));});
			}
			std::cout << oss.str() << std::endl;
		}
	}

	BOOST_AUTO_TEST_CASE(multi_serialization_small_xml) {
		multi::array<double, 2> d2D({10, 10});
		std::mt19937_64         e(std::random_device{}());

		//  auto g = std::bind(std::uniform_real_distribution<>{}, e);//
		auto g = [&]() {
			return std::uniform_real_distribution<>{}(e);
		};

		std::for_each(begin(d2D), end(d2D), [&](auto&& row) { std::generate(begin(row), end(row), g); });
		std::string const filename = "serialization-small.xml";
		{
			std::ofstream ofs{filename};
			assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		}
		{
			std::ifstream ifs{filename};
			assert(ifs);
			decltype(d2D) d2D_copy(extensions(d2D));
			boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
			BOOST_TEST( d2D_copy == d2D );
		}
		{
			std::ofstream ofs{"serialization-small-part.xml"};
			assert(ofs);
			auto&& a = d2D({0, 5}, {0, 5});
			boost::archive::xml_oarchive{ofs} << boost::serialization::make_nvp("d2D_part", a);  // BOOST_SERIALIZATION_NVP(d2D);
			fs::remove("serialization-small-part.xml");
		}
		std::cout << fs::file_size(filename) << '\n';
		fs::remove(filename);
	}

	BOOST_AUTO_TEST_CASE(multi_serialization_static_large_xml) {

		multi::dynamic_array<double, 2> d2D({1000, 1000});

		auto gen = [e = std::mt19937_64(std::random_device{}())]() mutable {
			return std::uniform_real_distribution<>{}(e);
		};
		std::for_each(begin(d2D), end(d2D), [&](auto&& row) { std::generate(begin(row), end(row), gen); });

		watch w("static_large_xml");

		std::string const filename = "serialization-static-large.xml";
		{
			std::ofstream ofs{filename};
			assert(ofs);
			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		}
		{
			std::ifstream ifs{filename};
			assert(ifs);
			decltype(d2D) d2D_copy(extensions(d2D));
			boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
			BOOST_TEST( d2D_copy == d2D );
		}
		std::cout << fs::file_size(filename) << '\n';
		fs::remove(filename);
	}

	BOOST_AUTO_TEST_CASE(multi_serialization_static_small) {
		{
			multi::dynamic_array<double, 0> d0D{12.0};

			std::ofstream ofs{"serialization-static_0D.xml"};
			assert(ofs);

			boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d0D);
			fs::remove("serialization-static_0D.xml");
		}
		{
			multi::array<double, 2> d2D = {
				{150.0, 16.0, 17.0, 18.0, 19.0},
				{  5.0,  5.0,  5.0,  5.0,  5.0},
				{100.0, 11.0, 12.0, 13.0, 14.0},
				{ 50.0,  6.0,  7.0,  8.0,  9.0},
			};
			auto gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable {
				return d(e);
			};
			std::for_each(
				begin(d2D), end(d2D),
				[&](auto&& row) { std::generate(begin(row), end(row), gen); }
			);
			std::string const filename = "serialization-small-double2D.xml";
			[&, _ = watch("xml write double")] {
				std::ofstream ofs{filename};
				assert(ofs);
				boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
			}();
			std::cerr << "size " << double(fs::file_size(filename)) / 1e6 << "MB\n";
			fs::remove(filename);
		}
		{
			multi::array<double, 2> d2D = {
				{150.0, 16.0, 17.0, 18.0, 19.0},
				{  5.0,  5.0,  5.0,  5.0,  5.0},
				{100.0, 11.0, 12.0, 13.0, 14.0},
				{ 50.0,  6.0,  7.0,  8.0,  9.0},
			};
			d2D.reextent({2000, 2000});
			auto gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable {
				return d(e);
			};
			std::for_each(
				begin(d2D), end(d2D),
				[&](auto&& r) { std::generate(begin(r), end(r), gen); }
			);
			[&, _ = watch("xml write double")] {
				std::ofstream ofs{"serialization-double.xml"};
				assert(ofs);
				boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
			}();
			std::cerr << "size " << double(fs::file_size("serialization-double.xml")) / 1e6 << "MB\n";
			fs::remove("serialization-double.xml");
		}

		using complex = std::complex<float>;

		auto const d2D = [] {
			multi::array<complex, 2> _({10000, 1000});
			auto                     gen = [d = std::uniform_real_distribution<double>{-1, 1}, e = std::mt19937{std::random_device{}()}]() mutable {
                return std::complex<double>{d(e), d(e)};
			};
			std::for_each(begin(_), end(_), [&](auto&& r) { std::generate(begin(r), end(r), gen); });
			return _;
		}();
		auto size = sizeof(double) * d2D.num_elements();
		using std::cerr;
		std::cout << "data size (in memory) " << size << std::endl;
		{
			fs::path file{"serialization.bin"};
			auto     count = [&, w = watch("binary write")] {
                std::ofstream ofs{file};
                assert(ofs);
                boost::archive::binary_oarchive{ofs} << d2D;
                return *w;
			}();
			cerr << "size  " << double(file_size(file)) / 1e6 << "MB\n";
			cerr << "speed " << double(size) / 1e6 / count << "MB/s\n";
			std::decay_t<decltype(d2D)> d2D_cpy;

			auto count_load = [&, w = watch("binary load")] {
				std::ifstream ifs{file};
				assert(ifs);
				boost::archive::binary_iarchive{ifs} >> d2D_cpy;
				return *w;
			}();

			std::cerr << "load speed " << double(file_size(file)) / 1e6 / count_load << "MB/s\n";
			BOOST_TEST( d2D == d2D_cpy );
			fs::remove(file);
		}
		// {
		//  using std::cout;
		//  fs::path file{"serialization.xml"};
		//  cout << file << std::endl;
		//  auto count = [&, w = watch("xml write base64")] {
		//      std::ofstream ofs{file};
		//      assert(ofs);
		//      boost::archive::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(d2D);
		//      return *w;
		//  }();
		//  cout << "data size " << double(size) / 1e6 << "MB\n";
		//  cout << "file size " << double(file_size(file)) / 1e6 << "MB\n";
		//  cout << "save speed " << double(size) / 1e6 / count << "MB/s" << std::endl;
		//  multi::array<complex, 2> d2D_cpy;

		//  auto count2 = [&, w = watch("xml load base64")] {
		//      std::ifstream ifs{file};
		//      assert(ifs);
		//      boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_cpy);
		//      return *w;
		//  }();

		//  cout << "load speed " << double(size) / 1e6 / count2 << "MB/s" << std::endl;
		//  BOOST_TEST( d2D_cpy == d2D );
		//  fs::remove(file);
		// }
		// #if 0
		//  {
		//      multi::cuda::managed::array<complex, 2> cud2D({2000, 2000});
		//      [&, _=watch("cuda binary write")]{
		//          std::ofstream ofs{"serialization.bin"}; assert(ofs);
		//          boost::archive::binary_oarchive{ofs} << cud2D;
		//      }();
		//      std::cerr<<"size "<< (fs::file_size("serialization.bin")/1e6) <<"MB\n";
		//  }
		// #endif
		// {
		//  [&, _ = watch("text write")] {
		//      std::ofstream ofs("serialization.txt");
		//      assert(ofs);
		//      boost::archive::text_oarchive{ofs} << d2D;
		//      assert(ofs);
		//  }();
		//  std::cerr << "size " << double(fs::file_size("serialization.txt")) / 1e6 << "MB\n";
		//  fs::remove("serialization.txt");
		// }
		// {
		//  multi::array<complex, 2> d2D_copy;  //(extensions(d2D), 9999.0);
		//  [&, _ = watch("text read")] {
		//      std::ifstream ifs("serialization.txt");
		//      assert(ifs);
		//      boost::archive::text_iarchive{ifs} >> d2D_copy;
		//      assert(ifs);
		//  }();
		//  BOOST_TEST( d2D_copy == d2D );
		//  fs::remove("serialization.txt");
		// }
		//  {
		//      multi::array<complex, 2> d2D_copy;  //(extensions(d2D), 9999.0);
		//      [&, _ = watch("binary read")] {
		//          std::ifstream ifs{"serialization.bin"};
		//          assert(ifs);
		//          boost::archive::binary_iarchive{ifs} >> d2D_copy;
		//      }();
		//      BOOST_TEST( d2D_copy == d2D );
		//      fs::remove("serialization.bin");
		//  }
		//  {
		//      [&, _ = watch("binary compressed write")] {
		//          std::ofstream ofs{"serialization_compressed.bin.gz"};
		//          {
		//              boost::iostreams::filtering_stream<boost::iostreams::output> f;
		//              f.push(boost::iostreams::gzip_compressor());
		//              f.push(ofs);
		//              boost::archive::binary_oarchive{f} << d2D;
		//          }
		//      }();
		//      std::cerr << "size " << double(fs::file_size("serialization.bin.gz")) / 1e6 << "MB\n";
		//      fs::remove("serialization.bin.gz");
		//  }
		{
			[&, _ = watch("compressed xml write")] {
				std::ofstream ofs{"serialization.xml.gz"};
				assert(ofs);
				{
					boost::iostreams::filtering_stream<boost::iostreams::output> f;
					f.push(boost::iostreams::gzip_compressor());
					f.push(ofs);
					boost::archive::xml_oarchive{f} << BOOST_SERIALIZATION_NVP(d2D);
				}
			}();
			std::cerr << "size " << double(fs::file_size("serialization.xml.gz")) / 1e6 << "MB\n";
			fs::remove("serialization.xml.gz");
		}
		//      {
		//          multi::array<complex, 2> d2D_copy;  //(extensions(d2D), 9999.);
		//          [&, _ = watch("xml read")] {
		//              std::ifstream ifs{"serialization.xml"};
		//              assert(ifs);
		//              boost::archive::xml_iarchive{ifs} >> BOOST_SERIALIZATION_NVP(d2D_copy);
		//          }();
		//          BOOST_TEST( d2D_copy == d2D );
		//          fs::remove("serialization.xml");
		//      }
	}

	BOOST_AUTO_TEST_CASE(test_utility_serialization_2d) {
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy types
		double carr[3][10] = {
			{ 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0},
			{10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0},
			{20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0},
		};
		multi::array_ref<double, 2> marr(&carr[0][0], {3, 10});
		//		boost::multi_array_ref<double, 2> Marr(&carr[0][0], boost::extents[3][10]);

		namespace arxiv = boost::archive;
		{
			std::ofstream ofs{"utility_serialization_marr.xml"};
			assert(ofs);
			arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(marr);
			fs::remove("utility_serialization_marr.xml");
		}
		{
			std::ofstream ofs{"utility_serialization_carr.xml"};
			assert(ofs);
			arxiv::xml_oarchive{ofs} << BOOST_SERIALIZATION_NVP(carr);
			fs::remove("utility_serialization_carr.xml");
		}
	}

	return boost::report_errors();
}
