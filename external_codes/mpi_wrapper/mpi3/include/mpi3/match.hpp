#ifndef BMPI3_MATCH_HPP
#define BMPI3_MATCH_HPP

#include "../mpi3/message.hpp"
#include "../mpi3/status.hpp"

#include<mpi.h>

namespace boost{
namespace mpi3{

#if not defined(EXAMPI)
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
#endif

}  // end namespace mpi3
}  // end namespace boost
#endif
