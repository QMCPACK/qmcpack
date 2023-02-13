// © Alfredo Correa 2018-2021

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world)-> int try {
	assert( world.size() == 6 );

	mpi3::communicator fifth = world/5;

	cout << "I am rank " << world.rank() << " in " << world.name() << ", ";

	if(fifth){cout<<"I am also   "<< fifth.rank() <<" in "<< fifth.name() <<'\n';}
	else     {cout<<"I am not in "<< fifth.name() <<'\n';}

	return 0;
} catch(...) {return 1;}
