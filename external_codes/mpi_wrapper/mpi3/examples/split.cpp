#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 -Wall `#-Wfatal-errors` $0 -o $0x.x -lboost_serialization && time mpirun -np 16 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator& world){

	{
	//	int color = world.rank() / 4;
	//	mpi3::communicator row = world.split(color, world.rank());
		mpi3::communicator row = world / 4; // shortcut
		cout << "world " << world.rank() << '/' << world.size() << " ---> row " << row.rank() << '/' << row.size() << '\n';
	}
	{
	//	int color = world.rank() % 4;
	//	mpi3::communicator col = world.split(color, world.rank());
		mpi3::communicator col = world % 4; // shortcut
		cout << "world " << world.rank() << '/' << world.size() << " ---> col " << col.rank() << '/' << col.size() << '\n';
	}

	return 0;
}

