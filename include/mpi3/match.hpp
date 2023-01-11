/* -*- indent-tabs-mode: t -*- */
#ifndef MPI3_MATCH_HPP
#define MPI3_MATCH_HPP

#include "../mpi3/message.hpp"
#include "../mpi3/status.hpp"

#include<mpi.h>

namespace boost{
namespace mpi3{

struct match : public message, public status {  // NOLINT(fuchsia-multiple-inheritance)
	friend class communicator;
	template<class It>
	auto receive(It dest){
		return receive_n(
			dest, 
			count<typename std::iterator_traits<It>::value_type>()
		);
	}
};

}  // end namespace mpi3
}  // end namespace boost
#endif
