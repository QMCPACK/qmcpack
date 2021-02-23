#if COMPILATION// -*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4; -*-
$CXXX `mpicxx -showme:compile|sed 's/-pthread/ /g'` -std=c++14 $0 -o $0x `mpicxx -showme:link|sed 's/-pthread/ /g'`&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#ifndef BOOST_MPI3_CORE_HPP
#define BOOST_MPI3_CORE_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<stdexcept>

namespace boost{
namespace mpi3{

inline bool initialized(){
	int flag = -1;
	int s = MPI_Initialized(&flag); 
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot probe initialization"};
	return flag;
}

inline bool finalized(){
	int flag = -1;
	int s = MPI_Finalized(&flag);
	if(s != MPI_SUCCESS) throw std::runtime_error{"cannot probe finalization"};
	return flag;
}

}}

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MPI3_ENVIRONMENT

#include<cassert>

namespace mpi3 = boost::mpi3;

int main(){
	assert(not mpi3::initialized() );
}

#endif
#endif

