#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*- */
mpicxx.mpich -D_TEST_BOOST_MPI3_VERSION -lboost_serialization -xc++ $0 -o $0x&&mpirun -np 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_VERSION_HPP
#define BOOST_MPI3_VERSION_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "detail/call.hpp"

#include<cassert>
#include<tuple> // tie
#include<iostream>

#define BOOST_MPI3_MAJOR_VERSION 0
#define BOOST_MPI3_MINOR_VERSION 71
#define BOOST_MPI3_PATCH_VERSION 0
#define BOOST_MPI3_VERSION_STRING "Boost.MPI3/0.71"

#define BOOST_MPI3_VERSION 0*100 + BOOST_MPI3_MINOR_VERSION*10

namespace boost{
namespace mpi3{

struct version_t{
	int major;
	int minor;
	version_t() = default;
	constexpr version_t(int major, int minor = 0) : major{major}, minor{minor}{}
	friend std::ostream& operator<<(std::ostream& os, version_t const& self){
		return os << self.major <<'.'<< self.minor;
	}
	constexpr bool operator<(version_t const& other) const{
		return std::tie(major, minor) < std::tie(other.major, other.minor);
	}
	constexpr bool operator>(version_t const& o) const{return o < *this;}
	constexpr bool operator==(version_t const& o) const{return not operator<(o) and not operator>(o);}
	constexpr bool operator>=(version_t const& o) const{return operator>(o) or operator==(o);}
	constexpr bool operator<=(version_t const& o) const{return operator<(o) or operator==(o);}
};

constexpr version_t Version(){return {MPI_VERSION, MPI_SUBVERSION};}

version_t version(){
	version_t ret;
	MPI_(Get_version)(&ret.major, &ret.minor);
	assert( ret == Version() );
	return ret;
}

std::string library_version(){
	int len;
	char mpi_lib_ver[MPI_MAX_LIBRARY_VERSION_STRING];
	MPI_(Get_library_version)(mpi_lib_ver, &len);
	return std::string(mpi_lib_ver, len);
}

std::string library_version_short(){
	std::string ret = library_version();
	{
		auto found = ret.find('\n');
		if(found != std::string::npos) ret = std::string(ret.c_str(), found);
	}
	{
		auto found = ret.find(',');
		if(found != std::string::npos) ret = std::string(ret.c_str(), found);
	}	
	return ret;
}

}}

#ifdef _TEST_BOOST_MPI3_VERSION

#include "../mpi3/main.hpp"
#include<cassert>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	assert(( mpi3::version() == mpi3::Version() ));
	assert(( mpi3::version() == mpi3::version_t{MPI_VERSION, MPI_SUBVERSION} ));
	assert(( mpi3::version() == mpi3::version_t{3, 1} ));
	assert(( mpi3::version() <  mpi3::version_t{3, 2} ));
	assert(( mpi3::version() >  mpi3::version_t{3, 0} ));
	if(world.rank() == 0){
		cout 
			<<"mpi Version                : "<< mpi3::Version()               <<'\n'
			<<"mpi version                : "<< mpi3::version()               <<'\n'
			<<"mpi3 library version       : "<< mpi3::library_version()       <<'\n'
			<<"mpi3 library version short : "<< mpi3::library_version_short() <<'\n'
			<<"mpi3 wrapper version       : "<< BOOST_MPI3_VERSION            <<'\n'
			<<"mpi3 wrapper version string: "<< BOOST_MPI3_VERSION_STRING     <<'\n'
		;
	}
	assert( BOOST_MPI3_VERSION >= 071 );
	return 0;

}
#endif
#endif

