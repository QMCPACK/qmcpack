// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/ostream.hpp"

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {

	mpi3::ostream wout(world);
	
	mpi3::communicator third = world/3; // or other division
	mpi3::communicator leaders = world.keep(third.root()); // same as world.split(third.root()?0:mpi3::undefined);
 
	wout << "I am 'world' rank "<<world.rank(); 
	if(third){
		wout<<" and 'third':"<< third.name() <<"'s rank "<<third.rank() <<" with color attribute "<< mpi3::any_cast<int>(third.attribute("color"));
	}else{
		wout<<" and not in 'third'";
	}
	if(leaders){
		wout<<" and 'leader:'"<< leaders.name() <<"'s rank "<< leaders.rank() <<" with color attribute "<< mpi3::any_cast<int>(third.attribute("color"));
	}else{
		wout <<" and not in 'leader'";
	}
	wout << std::endl;

	return 0;
} catch(...) {return 1;}
