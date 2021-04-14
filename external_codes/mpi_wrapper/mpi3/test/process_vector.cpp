#if COMPILATION_INSTRUCTIONS
mpicxx -O3 -std=c++14 $0 -o $0.x -lboost_serialization -lboost_container&&time mpirun -np 2 $0.x&&rm -f $0.x;exit
#endif

#include "../main.hpp"
#include "../process.hpp"

#include <boost/serialization/vector.hpp>

struct long_long{
	long long value;
	long_long& operator=(long long v){value = v; return *this;}
};

template<class Archive>
void serialize(Archive& ar, long_long& l, unsigned = 0){
	ar & l.value;
}

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	assert( world.size() == 2 );
	long long size = 100000000;
	switch(world.rank()){
	case 0:
		{
			std::vector<long_long> v(size); std::iota(v.begin(), v.end(), 0.);
		//	assert(std::accumulate(v.begin(), v.end(), 0.) == size*(size-1)/2 );
			world[1] << v;
		}
		break;
	case 1:
		{
			std::vector<long_long> w; world[0] >> w;
			assert( w.size() == size );
			assert( w[45].value = 45 );
		//	assert(std::accumulate(w.begin(), w.end(), 0.) == size*(size-1)/2 );
		}
		break;
	}
	return 0;

}

