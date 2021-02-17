#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -O3 -Wall -Wextra $0 -o $0x.x -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY `#-lboost_serialization` -pthread && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

#include<complex>
#include<future> // async
#include<string>
#include<list>

namespace mpi3 = boost::mpi3;
using std::cout;

using namespace std::chrono_literals;

int mpi3::main(int, char*[], mpi3::communicator world){
	assert( world.size() > 1 );

	{
		switch(world.rank()){
			case 0: {
				std::list<int> b = {3, 4, 5};
				world.send(cbegin(b), cend(b), 1);
			}; break;
			case 1: {
				std::vector<int> b2(3);
				auto e = world.receive(begin(b2), 0);
				assert( e == end(b2) );
				assert( b2[1] == 4. );
			}; break;
		}
	}
	{
		switch(world.rank()){
			case 0: {
				std::list<int> b = {3, 4, 5};
				std::this_thread::sleep_for(10s);
				auto req = std::async([&](){return world.send(cbegin(b), cend(b), 1);});
			}; break;
			case 1: {
				std::vector<int> b2(3);
				auto req = std::async([&](){return world.receive(begin(b2), 0);});
				std::this_thread::sleep_for(2s);
				std::cout << "rank " << world.rank() << " has the req" << std::endl;
				assert( req.get() == end(b2) );
				assert( b2[1] == 4. );
			}; break;
		}
	}
	return 0;
}

