#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_MPI3_DETAIL_BUFFER $0x.cpp -o $0x.x && time mpirun -n 3 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_DETAIL_BUFFER_HPP
#define MPI3_DETAIL_BUFFER_HPP

//#include "../../mpi3/communicator_fwd.hpp"
//#include "../../mpi3/detail/iterator.hpp"
#include "../../mpi3/detail/datatype.hpp" // packed
#include "../../mpi3/vector.hpp"


namespace boost{
namespace mpi3{
namespace detail{

struct buffer : mpi3::uvector<detail::packed>{
	int pos = 0;
	buffer() = default;
	buffer(std::size_t r){ reserve(r); }
	buffer(buffer const&) = delete;
};

}}}

#ifdef _TEST_MPI3_BUFFER

#include "../mpi3/communicator.hpp"
#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
}

#endif
#endif

