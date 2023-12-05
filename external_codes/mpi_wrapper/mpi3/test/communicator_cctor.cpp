// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"

#include<complex>
#include<future>
#include<string>

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {

	static_assert(sizeof(MPI_Comm) == sizeof(mpi3::communicator) );

{
	mpi3::communicator& w2 = mpi3::grip_communicator(world.handle());
	assert(  w2 ==  world );
	assert( &w2 == &world );

	assert(  mpi3::grip_communicator(MPI_COMM_WORLD) ==  world );
	assert( &mpi3::grip_communicator(MPI_COMM_WORLD) != &world );
}

//	assert( reinterpret_cast<mpi3::communicator&>(MPI_COMM_WORLD) == world );

{
	mpi3::communicator mty;
	assert( mty.empty() );
//  assert( mty.size() == 0 );
	mpi3::communicator mty2 = mty;
	assert( mty2.empty() );
//  assert( mty2.size() == 0 );
}
{
	std::vector<mpi3::communicator> comms;
	comms.emplace_back(world);
	comms.emplace_back(world);
	comms.emplace_back(world);

	std::vector<mpi3::communicator> comms2;
	comms2.reserve(3);
//  for(auto& e:comms) {comms2.emplace_back(e);}  // ok, but old style
//  std::copy(comms.begin(), comms.end(), std::back_inserter(comms2));  // doesn't work because it calls cctor
	std::transform(comms.begin(), comms.end(), std::back_inserter(comms2), [](auto&& e) {return e;});  // calls dup ctor
}
{
	std::size_t const NTHREADS = 10;
	std::vector<std::future<int>> fs;
	for(int i=0; i != NTHREADS; ++i) {  // NOLINT(altera-unroll-loops)
#if 0 // this is problematic because copy (mpi_comm_dup) is not thread safe
		fs.emplace_back(std::async([&world](){
			auto comm = world; // hangs here
			return 5;
		}));
#endif
		fs.emplace_back(std::async([thread_comm=world](){
		//	auto comm2 = thread_comm; // works, just to test
			return 5;
		}));
#if 0 // more correct
		fs.emplace_back(std::async([world](){
			auto comm2 = world; // works, just to test
			return 5;
		}));
#endif
		std::cout << "created thread" << std::endl;
	}
	for(std::size_t i = 0; i != NTHREADS; ++i) {  // NOLINT(altera-unroll-loops)
		auto five = fs[i].get();
		assert( five == 5 );
		std::cout << "joined thread" << std::endl;
	}
}
	return 0;
} catch(...) {return 1;}
