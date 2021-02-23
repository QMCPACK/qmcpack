#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_MPI3_FUTURE $0x.cpp -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef MPI3_FUTURE_HPP
#define MPI3_FUTURE_HPP

#include "../mpi3/communicator.hpp"

namespace boost{
namespace mpi3{

template<class T>
struct future{
	template<class F>
	future(communicator& comm, F&& f, Args&&... args){
		if(not comm.root()){
			int r = comm.rank();
			comm.send_n(&r, 1, 0);
		}else{
			int r = -1;
			std::vector<request> reqs;
			req
			auto req = comm.ireceive_n(&r, 1);
		}
	}
};

// futures factory, todo implement launch policy
template<class F, class... Args> auto async(
	communicator& c, F&& f, Args&&... args
) -> future<decltype(F(std::forward(args)...)>{
	return {c, std::forward<F>(f), std::forward<Args>(args)...};
}

}}

#ifdef _TEST_MPI3_FUTURE

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){
	return 0;
}

#endif
#endif

