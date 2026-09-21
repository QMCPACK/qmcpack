// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/detail/mpi_impl.h>

#include<type_traits> // make_signed_t

namespace boost {
namespace mpi3 {

using size_t = MPI_Aint;
using address = std::make_signed_t<size_t>;

}  // end namespace mpi3
}  // end namespace boost
