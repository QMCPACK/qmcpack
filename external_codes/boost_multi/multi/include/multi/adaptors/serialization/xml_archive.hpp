#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2021
#include "../../utility.hpp"

#include<boost/archive/xml_iarchive.hpp>
#include<boost/archive/xml_oarchive.hpp>

namespace boost{
namespace multi{

template<>
struct archive_traits<boost::archive::xml_oarchive>{
	template<class T> static auto make_nvp(char const* name, T& value) -> decltype(auto){
		return boost::serialization::make_nvp(name, value);
	} 
};
template<>
struct archive_traits<boost::archive::xml_iarchive>{
	template<class T> static auto make_nvp(char const* name, T& value) -> decltype(auto){
		return boost::serialization::make_nvp(name, value);
	}
};

} // end namespace multi
} // end namespace boost

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi adaptor serialization xml_archive"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_adaptor_serialization_xml_archive){
	BOOST_REQUIRE(true);
}
#endif

