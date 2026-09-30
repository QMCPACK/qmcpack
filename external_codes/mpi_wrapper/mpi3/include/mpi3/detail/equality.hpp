// Copyright 2018-2025 Alfredo A. Correa

#ifndef BOOST_MPI3_DETAIL_EQUALITY_HPP
#define BOOST_MPI3_DETAIL_EQUALITY_HPP

#include <mpi3/detail/mpi_impl.h>

#include<type_traits> // underglying_type

namespace boost{
namespace mpi3{
namespace detail{

static_assert( sizeof(int)>=sizeof(MPI_IDENT), "standard requires");

enum class equality : int {  // NOLINT(performance-enum-size)
	identical = MPI_IDENT,     // same (same address)
	congruent = MPI_CONGRUENT, // equal (in the value sense)
	similar   = MPI_SIMILAR,   // same processes in different order (permutation)
	unequal   = MPI_UNEQUAL    // not equal in any sense
};

static_assert(sizeof(equality) >= sizeof(MPI_IDENT));

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

#endif
