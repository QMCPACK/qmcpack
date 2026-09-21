// Copyright 2018-2025 Alfredo A. Correa

#ifndef BMPI3_MATCH_HPP
#define BMPI3_MATCH_HPP

#include <mpi3/message.hpp>
#include <mpi3/status.hpp>

#include <mpi3/detail/mpi_impl.h>

namespace boost{
namespace mpi3{

#ifndef EXAMPI
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
