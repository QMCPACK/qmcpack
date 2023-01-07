#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

#include<list>
#include<vector>

namespace mpi3 = boost::mpi3;

void f(mpi3::communicator      & /*comm*/) {}
void g(mpi3::communicator        /*comm*/) {}
void h(mpi3::communicator const& /*comm*/) {}

auto ovrld(mpi3::communicator        /*comm*/) -> std::string {return "by value";}
auto ovrld(mpi3::communicator      & /*comm*/) -> std::string {return "by reference";}
//auto ovrld(mpi3::communicator const& /*comm*/) -> std::string {return "by const reference";}

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {

	f(world);
	//  f(world.duplicate());  // doesn't work, good

	g(world.duplicate());
	g(world);  // works, implicit calls to duplicate

	h(world);
	h(world.duplicate());

//  assert( ovrld(world) == "by ???" );  // ambiguous, not much to do, overload by reference can never called
//  assert( ovrld(std::ref(world)) == "by ???" ); // ambiguous
	assert( ovrld(world.duplicate()) == "by value" );
	assert( ovrld(mpi3::communicator{world}) == "by value" );

	return 0;
} catch(...) {return 1;}

