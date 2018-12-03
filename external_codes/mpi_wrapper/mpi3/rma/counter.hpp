#if COMPILATION_INSTRUCTIONS
(echo "#include<"$0">" > $0x.cpp) && mpicxx -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_COUNTER $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_RMA_COUNTER_HPP
#define BOOST_MPI3_RMA_COUNTER_HPP


#include "../../mpi3/window.hpp"

namespace boost{
namespace mpi3{

template<class T>
class counter{
	mpi3::window w_;
	counter(mpi3::communicator& c, T num){
		int lnum = num / c.size();
		int lleft = num % c.size();
		if(rank < lleft) ++lnum;
		
	}
};

}}

#ifdef _TEST_BOOST_MPI3_COUNTER

#include "alf/boost/mpi3/main.hpp"

using std::cout;

int boost::mpi3::main(int argc, char* argv[], boost::mpi3::communicator& world){
	
}

#endif
#endif

