/* -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*- */
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_VERSION_HPP
#define BOOST_MPI3_VERSION_HPP

#include <mpi.h>

#include "detail/call.hpp"

#include <cassert>
#include <iostream>
#include <tuple>  // tie

#define BOOST_MPI3_MAJOR_VERSION 0   // NOLINT(cppcoreguidelines-macro-usage,modernize-macro-to-enum)
#define BOOST_MPI3_MINOR_VERSION 81  // NOLINT(cppcoreguidelines-macro-usage,modernize-macro-to-enum)
#define BOOST_MPI3_PATCH_VERSION 0   // NOLINT(cppcoreguidelines-macro-usage,modernize-macro-to-enum)

#define BOOST_MPI3_VERSION_STRING "Boost.MPI3/0.81"  // NOLINT(cppcoreguidelines-macro-usage)

#define BOOST_MPI3_VERSION (BOOST_MPI3_MAJOR_VERSION * 100 + BOOST_MPI3_MINOR_VERSION * 10)

namespace boost {
namespace mpi3 {

struct version_t {
	int major;  // NOLINT(misc-non-private-member-variables-in-classes)
	int minor;  // NOLINT(misc-non-private-member-variables-in-classes)
//	version_t() = default;
//	constexpr explicit version_t(int major, int minor = 0) : major{major}, minor{minor} {}  // NOLINT(bugprone-easily-swappable-parameters)
	friend std::ostream& operator<<(std::ostream& os, version_t const& self) {
		return os << self.major << '.' << self.minor;
	}
	constexpr bool operator<(version_t const& other) const {
		return std::tie(major, minor) < std::tie(other.major, other.minor);
	}
	constexpr bool operator>(version_t const& o) const { return o < *this; }

	constexpr bool operator==(version_t const& o) const { return not operator<(o) and not operator>(o); }
	constexpr bool operator!=(version_t const& o) const { return not operator==(o);}

	constexpr bool operator>=(version_t const& o) const { return operator>(o) or operator==(o); }
	constexpr bool operator<=(version_t const& o) const { return operator<(o) or operator==(o); }
};

inline constexpr auto Version() -> version_t { return {MPI_VERSION, MPI_SUBVERSION}; }  // NOLINT(readability-identifier-naming)
static constexpr auto VERSION = version_t{MPI_VERSION, MPI_SUBVERSION};

inline auto version() {
	version_t ret{};
	MPI_(Get_version)(&ret.major, &ret.minor);
	return ret;
}

inline auto library_version() -> std::string {
	std::array<char, MPI_MAX_LIBRARY_VERSION_STRING> mpi_lib_ver{};
	int  len = 0;
	MPI_(Get_library_version)(mpi_lib_ver.data(), &len);
	return {mpi_lib_ver.data(), static_cast<std::string::size_type>(len)};
}

inline auto library_version_short() -> std::string {
	std::string ret = library_version();
	{
		auto found = ret.find('\n');
		if(found != std::string::npos) {ret = std::string(ret.c_str(), found);}
	}
	{
		auto found = ret.find(',');
		if(found != std::string::npos) {ret = std::string(ret.c_str(), found);}
	}
	return ret;
}

}  // namespace mpi3
}  // namespace boost

//#ifdef _TEST_BOOST_MPI3_VERSION

//#include "../mpi3/main.hpp"
//#include<cassert>

// namespace mpi3 = boost::mpi3;
// using std::cout;

// int mpi3::main(int, char*[], mpi3::communicator world){

//	assert(( mpi3::version() == mpi3::Version() ));
//	assert(( mpi3::version() == mpi3::version_t{MPI_VERSION, MPI_SUBVERSION} ));
//	assert(( mpi3::version() == mpi3::version_t{3, 1} ));
//	assert(( mpi3::version() <  mpi3::version_t{3, 2} ));
//	assert(( mpi3::version() >  mpi3::version_t{3, 0} ));
//	if(world.rank() == 0){
//		cout
//			<<"mpi Version                : "<< mpi3::Version()               <<'\n'
//			<<"mpi version                : "<< mpi3::version()               <<'\n'
//			<<"mpi3 library version       : "<< mpi3::library_version()       <<'\n'
//			<<"mpi3 library version short : "<< mpi3::library_version_short() <<'\n'
//			<<"mpi3 wrapper version       : "<< BOOST_MPI3_VERSION            <<'\n'
//			<<"mpi3 wrapper version string: "<< BOOST_MPI3_VERSION_STRING     <<'\n'
//		;
//	}
//	assert( BOOST_MPI3_VERSION >= 071 );
//	return 0;

//}
//#endif
#endif
