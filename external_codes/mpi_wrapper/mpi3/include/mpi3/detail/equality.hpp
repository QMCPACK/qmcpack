/* -*- indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4;-*- */
//#if COMPILATION_INSTRUCTIONS
//mpic++ -D_TEST_BOOST_MPI3_DETAIL_EQUALITY -xc++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
//#endif
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MPI3_DETAIL_EQUALITY_HPP
#define BOOST_MPI3_DETAIL_EQUALITY_HPP

// #define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<type_traits> // underglying_type

namespace boost{
namespace mpi3{
namespace detail{

static_assert( sizeof(int)>=sizeof(MPI_IDENT), "standard requires");

enum equality : int{
	identical = MPI_IDENT,     // same (same address)
	congruent = MPI_CONGRUENT, // equal (in the value sense)
	similar   = MPI_SIMILAR,   // same processes in different order (permutation)
	unequal   = MPI_UNEQUAL    // not equal in any sense
};

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_DETAIL_EQUALITY

//namespace mpi3 = boost::mpi3;

//int main(){

//	mpi3::detail::equality q = mpi3::detail::identical;
//	q = mpi3::detail::congruent;
//	q = mpi3::detail::similar;
//	q = mpi3::detail::unequal;

//}

//#endif
#endif
