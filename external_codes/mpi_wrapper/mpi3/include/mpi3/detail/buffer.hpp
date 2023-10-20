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
	int pos = 0;  // NOLINT(misc-non-private-member-variables-in-classes) TODO(correaa) : make private
	buffer() = default;
	explicit buffer(std::size_t r) {reserve(r);}
	buffer(buffer const&) = delete;
	buffer(buffer&&) = delete;
	buffer& operator=(buffer const&) = delete;
	buffer& operator=(buffer&&) = delete;
	~buffer() = default;
};

} // end namespace detail
} // end namespace mpi3
} // end namespace boost

//#ifdef _TEST_MPI3_BUFFER

//#include "../mpi3/communicator.hpp"
//#include "../mpi3/main.hpp"

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int, char*[], mpi3::communicator world){
//}

//#endif
#endif

