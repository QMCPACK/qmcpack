#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra $0 -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/environment.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

#define N 8

int main(){

	mpi3::environment env;
	mpi3::communicator world = env.world();

	int table[N][N];
	
	assert(world.size() == N);

	if(world.rank() == 0)
		for(int i = 0; i != N; ++i)
			for(int j = 0; j != N; ++j)
				table[i][j] = i + j;

	{
		int row[N];
		world.scatter(&table[0][0], &table[0][0] + N*N, &row[0], 0);
		for(int i = 0; i != N; ++i) assert(row[i] == i + world.rank());
	}
	{
		int row[N];
		world.scatter_from(&row[0], &row[0] + N, &table[0][0], 0);
		for(int i = 0; i != N; ++i) assert(row[i] == i + world.rank());
	}

	return 0;
}

