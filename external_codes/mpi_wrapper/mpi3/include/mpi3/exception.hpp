#ifndef BOOST_MPI3_EXCEPTION
#define BOOST_MPI3_EXCEPTION

#include<stdexcept>

namespace boost {
namespace mpi3 {

struct exception : std::runtime_error {
	using runtime_error::runtime_error;
//	~exception() override = default;
};

struct invalid_communicator : exception {
	using exception::exception;
//	~invalid_communicator() override = default;
};

}  // end namespace mpi3
}  // end namespace boost
#endif
