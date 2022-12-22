#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -I$HOME/include -Wall -Wextra -Wpedantic `#-Wfatal-errors` $0 -o $0x.x && time mpirun -n 8 $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../mpi3/environment.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"
#include<numeric> // iota
#include<list>

//#include<boost/range.hpp>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	std::cerr << sizeof(MPI_Offset) << " " << sizeof(int) << std::endl;
	#if 0
	if(0){
		if(auto duo = (world < 2)){
			assert( duo.size() == 2 );
			std::vector<double> v = {1, 2, 2};
			std::vector<int> ns = {1, 2}; if(!duo.root()) ns.clear();
			std::vector<std::vector<double>::iterator> its = {v.begin(), v.begin()+1}; if(!world.root()) its.clear();
			std::vector<double> w(duo.rank()+1);
			duo.scatterv_n(begin(its), begin(ns), begin(w));
			switch(duo.rank()){
				case 0: assert(( w == std::vector<double>{1} )); break;
				case 1: std::cerr <<"w2:" << w[0] << ", " << w[1] << std::endl; break;
			}
		}
	}
	#endif
	#if 1
	{
		if(auto quintet = (world < 5)){
			std::vector<std::vector<double>> vs = { {1}, {2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {5, 5, 5, 5, 5} }; if(!quintet.root()) vs.clear();
			std::vector<double> w( quintet.rank() + 1 );
			auto e = quintet.scatterv(vs, begin(w)); assert( e == end(w) );
			assert( w == std::vector<double>(quintet.rank() + 1, quintet.rank() + 1) );
		}
	}
	#endif

	return 0;
}

