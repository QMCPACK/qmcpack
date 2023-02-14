#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++17 -Wall -Wextra -Wpedantic -D_TEST_MPI3_RMA_MEMORY $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// (C) Copyright 2019 Alfredo A. Correa
#ifndef MPI3_RMA_MEMORY_HPP
#define MPI3_RMA_MEMORY_HPP

#include "../../mpi3/window.hpp"

namespace boost{
namespace mpi3{
namespace rma{


class mutex{
//	window<void>* w_;
	mutex(mutex const&) = delete;
};


}}}

#ifdef _TEST_MPI3_RMA_MEMORY

#include "../../mpi3/main.hpp"
#include "../../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3; 

int mpi3::main(int, char*[], mpi3::communicator world){



	return 0;
}

#endif


#endif
