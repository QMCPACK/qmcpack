// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <string>

namespace mpi3 = boost::mpi3;

namespace  {

void f(mpi3::communicator      & /*comm*/) {}
void g(mpi3::communicator        /*comm*/) {}
void h(mpi3::communicator const& /*comm*/) {}

auto ovrld(mpi3::communicator        /*comm*/) -> std::string {return "by value";}
// auto ovrld(mpi3::communicator      & /*comm*/) -> std::string {return "by reference";}
// auto ovrld(mpi3::communicator const& /*comm*/) -> std::string {return "by const reference";}

}  // end namespace

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	f(world);
	//  f(world.duplicate());  // doesn't work, good

	g(world.duplicate());
	g(world);  // works, implicit calls to duplicate

	h(world);
	h(world.duplicate());

//  BOOST_TEST( ovrld(world) == "by ???" );  // ambiguous, not much to do, overload by reference can never called
//  BOOST_TEST( ovrld(std::ref(world)) == "by ???" ); // ambiguous
	BOOST_TEST( ovrld(world.duplicate()) == "by value" );
	BOOST_TEST( ovrld(mpi3::communicator{world}) == "by value" );

	return boost::report_errors();
} catch(...) {
	return 1;
}
