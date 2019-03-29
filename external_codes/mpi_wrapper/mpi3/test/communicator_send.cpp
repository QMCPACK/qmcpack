#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -O3 -Wall -Wextra $0 -o $0x.x -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` && time mpirun -n 3 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

#include<complex>
#include<string>
#include<list>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() > 1 );

	switch(world.rank()){
		case 0: {
			std::list<int> b = {3, 4, 5};
			world.send(cbegin(b), cend(b), 1);
		}; break;
		case 1: {
			std::vector<int> b2(3);
			world.receive(begin(b2), 0);
			assert( b2[1] == 4. );
		}; break;
	}
	switch(world.rank()){
		case 0: {
			std::vector<std::string> b = {"hola", "blah", "chau"};
			world.send(cbegin(b), cend(b), 1);
		}; break;
		case 1: {
			std::list<std::string> b2(3);
			world.receive(begin(b2));//, 0);
			assert( *begin(b2) == "hola" and *rbegin(b2) == "chau" );
		}; break;
	}
	switch(world.rank()){
		case 0: {
			std::istringstream iss{"1 2 3"};
			world.send(std::istream_iterator<int>{iss}, {}, 1);
		}; break;
		case 1: {
			std::vector<int> out(3);
			world.receive(begin(out), 0);
			assert((out == std::vector<int>{1,2,3}));
		} break;
	}

	return 0;
}

