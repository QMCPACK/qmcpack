//#define OMPI_SKIP_MPICXX 1 // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<type_traits> // make_signed_t

namespace boost {
namespace mpi3 {

using size_t = MPI_Aint;
using address = std::make_signed_t<size_t>;

}  // end namespace mpi3
}  // end namespace boost
