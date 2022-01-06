#include <benchmark/benchmark.h>
#include "../array.hpp"

//#include <cereal/archives/json.hpp>
//#include <cereal/archives/xml.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

//#include "/home/correaa/prj/alf/boost/archive/yml/yml_iarchive.hpp"
//#include "/home/correaa/prj/alf/boost/archive/yml/yml_oarchive.hpp"
//#include <cereal/archives/portable_binary.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

//#include <cereal/archives/binary.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include<fstream>
#include<sstream>

#include<random>

constexpr auto N = 64;

namespace multi = boost::multi;

template<class Ar>
void BM_oserialization(benchmark::State& st) {
	auto const A = [] {
	    std::random_device rd;
	    std::mt19937 mt(rd());
	    std::uniform_real_distribution<double> dist(-1.0, +1.0);

		multi::array<double, 4> A({N, N, N, N});
		std::generate(begin(elements(A)), end(elements(A)), [&]{return dist(mt);});
		return A;
	}();

	std::cerr<< A.num_elements()*sizeof(double) <<std::endl;

    benchmark::ClobberMemory();
	for(auto _  : st) {
		std::ofstream fs{"file"};
		Ar xa{fs};
		xa<< multi::archive_traits<Ar>::make_nvp("A", A);
		fs.flush();
		benchmark::DoNotOptimize(A);
	    benchmark::ClobberMemory();
	}
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(double));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

template<class Ar>
void BM_iserialization(benchmark::State& st) {
	multi::array<double, 4> A({N, N, N, N});

    benchmark::ClobberMemory();
	for(auto _  : st) {
		std::ifstream fs{"file"};
		Ar xa{fs};
		xa>> multi::archive_traits<Ar>::make_nvp("A", A);
		benchmark::DoNotOptimize(A);
	    benchmark::ClobberMemory();
	}
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(double));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

BENCHMARK_TEMPLATE(BM_oserialization, boost::archive::xml_oarchive   );
BENCHMARK_TEMPLATE(BM_iserialization, boost::archive::xml_iarchive   );

//BENCHMARK_TEMPLATE(BM_oserialization, boost::archive::yml_oarchive   );
//BENCHMARK_TEMPLATE(BM_iserialization, boost::archive::yml_iarchive   );

//BENCHMARK_TEMPLATE(BM_oserialization, boost::archive::text_oarchive  );
//BENCHMARK_TEMPLATE(BM_iserialization, boost::archive::text_iarchive  );

//BENCHMARK_TEMPLATE(BM_oserialization, boost::archive::binary_oarchive);
//BENCHMARK_TEMPLATE(BM_iserialization, boost::archive::binary_iarchive);

//BENCHMARK_TEMPLATE(BM_oserialization, cereal::JSONOutputArchive  );
//BENCHMARK_TEMPLATE(BM_iserialization, cereal::JSONInputArchive   );

//BENCHMARK_TEMPLATE(BM_oserialization, cereal::XMLOutputArchive  );
//BENCHMARK_TEMPLATE(BM_iserialization, cereal::XMLInputArchive   );

//BENCHMARK_TEMPLATE(BM_oserialization, cereal::PortableBinaryOutputArchive  );
//BENCHMARK_TEMPLATE(BM_iserialization, cereal::PortableBinaryInputArchive   );

//BENCHMARK_TEMPLATE(BM_oserialization, cereal::BinaryOutputArchive  );
//BENCHMARK_TEMPLATE(BM_iserialization, cereal::BinaryInputArchive   );

template<class Ar>
void BM_gzip_oserialization(benchmark::State& st) {
	auto const A = []{
	    std::random_device rd;
	    std::mt19937 mt(rd());
	    std::uniform_real_distribution<double> dist(-1.0, +1.0);

		multi::array<double, 4> A({N, N, N, N});
		std::generate(begin(elements(A)), end(elements(A)), [&]{return dist(mt);});

		return A;
	}();

    benchmark::ClobberMemory();
	for(auto _  : st) {
		std::ofstream ofs{"file.gz"};
		boost::iostreams::filtering_ostream out;
		out.push(boost::iostreams::gzip_compressor{boost::iostreams::zlib::best_speed});
		out.push(ofs);
		Ar xa(out);
		xa<< multi::archive_traits<Ar>::make_nvp("A", A);
		benchmark::DoNotOptimize(A);
	    benchmark::ClobberMemory();
	}
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(double));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

template<class Ar>
void BM_gzip_iserialization(benchmark::State& st) {
	multi::array<double, 4> A({N, N, N, N});

    benchmark::ClobberMemory();
	for(auto _  : st) {
		std::ifstream ifs{"file.gz"};
		boost::iostreams::filtering_istream in;
		in.push(boost::iostreams::gzip_decompressor{});
		in.push(ifs);
		Ar xa(in);
		xa>> multi::archive_traits<Ar>::make_nvp("A", A);
		benchmark::DoNotOptimize(A);
	    benchmark::ClobberMemory();
	}
	st.SetBytesProcessed(st.iterations()*A.num_elements()*sizeof(double));
	st.SetItemsProcessed(st.iterations()*A.num_elements());
}

//BENCHMARK_TEMPLATE(BM_gzip_oserialization, cereal::XMLOutputArchive);
//BENCHMARK_TEMPLATE(BM_gzip_iserialization, cereal::XMLInputArchive);

BENCHMARK_TEMPLATE(BM_gzip_oserialization, boost::archive::xml_oarchive);
BENCHMARK_TEMPLATE(BM_gzip_iserialization, boost::archive::xml_iarchive);

BENCHMARK_MAIN();
