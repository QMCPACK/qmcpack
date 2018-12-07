#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 `#-Wfatal-errors` -D_TEST_MPI3_MESSAGE $0x.cpp -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_MATCH_HPP
#define MPI3_MATCH_HPP

#include "../mpi3/message.hpp"
#include "../mpi3/status.hpp"

#include<mpi.h>

namespace boost{
namespace mpi3{

struct match : public message, public status{
	friend class communicator;
	template<class It>
	auto receive(It dest){
		return receive_n(
			dest, 
			count<typename std::iterator_traits<It>::value_type>()
		);
	}
};

}}

#ifdef _TEST_MPI3_MATCH

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
}

#endif
#endif

