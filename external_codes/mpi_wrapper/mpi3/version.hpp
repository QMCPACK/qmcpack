#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_VERSION -lboost_serialization $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.cpp $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef BOOST_MPI3_VERSION_HPP
#define BOOST_MPI3_VERSION_HPP
#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<boost/operators.hpp>

#include<cassert>
#include<tuple> // tie
#include<iostream>

namespace boost{
namespace mpi3{

struct version_t{// : private boost::totally_ordered<version_t>{
	int major;
	int minor;
	constexpr version_t() : major{}, minor{}{}
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

constexpr version_t Version(){
	version_t ret{MPI_VERSION, MPI_SUBVERSION};
	return ret;
}

version_t version(){
	version_t ret;
	MPI_Get_version(&ret.major, &ret.minor);
	assert( ret == Version() );
	return ret;
}

std::string library_version(){
	int len;
    char mpi_lib_ver[MPI_MAX_LIBRARY_VERSION_STRING];
    MPI_Get_library_version(mpi_lib_ver, &len);
    return std::string(mpi_lib_ver, len);
}
std::string library_version_short(){
	std::string ret = library_version();
	auto found = ret.find('\n');
	if(found != std::string::npos) return std::string(ret.c_str(), found);
	return ret;
}

}}

#ifdef _TEST_BOOST_MPI3_VERSION

#include "../mpi3/main.hpp"

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
			<<"size "<< world.size() <<'\n'
			<<"mpi version "<< mpi3::version() <<'\n'
		;
	}
	cout << mpi3::library_version_short() <<'\n';
	return 0;
}
#endif
#endif

