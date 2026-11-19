// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#ifndef EXAMPI
#include <mpi3/ostream.hpp>
#endif

#include <iostream>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	mpi3::communicator third = world/3; // or other division
	mpi3::communicator const leaders = world.keep(third.root()); // same as world.split(third.root()?0:mpi3::undefined);

#ifndef EXAMPI
	mpi3::ostream wout(world);
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
	wout << '\n' << std::flush;
#endif
}
