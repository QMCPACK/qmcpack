// Copyright 2018-2026 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include<list>
#include<vector>

namespace mpi3 = boost::mpi3;

struct projector {
	explicit projector(mpi3::communicator& comm) : comm_{comm} {}
#ifdef _MSC_VER
	projector(projector const& other) : comm_{other.comm_} {}  // mutable member: const& gives non-const access. Necessary for MSVC
#endif

 private:
	mutable mpi3::communicator comm_;
};

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	{
		std::list<mpi3::communicator> v;
		v.emplace_back(world);
		v.emplace_back(world);
	}
// #if 0
//  { // doesn't compile, communicator is not copiable
//      std::vector<mpi3::communicator> v = {world, world};
//      v.emplace_back(world);
//      v.emplace_back(world);
//  }
// #endif
	{
		std::vector<projector> v = {projector{world}, projector{world}};
		v.emplace_back(world);
		v.emplace_back(world);
	}
}
