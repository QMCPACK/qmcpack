#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x&&mpirun --oversubscribe -n 8 $0x&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

#include<list>
#include<vector>

namespace mpi3 = boost::mpi3;

struct projector{
	projector(mpi3::communicator& comm) : comm_{comm}{}
	projector(projector const&) = default;
private:
	mutable mpi3::communicator comm_;
};

int mpi3::main(int, char*[], mpi3::communicator world){

	{
		std::list<mpi3::communicator> v;
		v.emplace_back(world);
		v.emplace_back(world);
	}
//	{ // doesn't compile, communicator is not copiable
//		std::vector<mpi3::communicator> v = {world, world};
//		v.emplace_back(world);
//		v.emplace_back(world);
//	}
	{
		std::vector<projector> v = {world, world};
		v.emplace_back(world);
		v.emplace_back(world);
	}

	return 0;
}

