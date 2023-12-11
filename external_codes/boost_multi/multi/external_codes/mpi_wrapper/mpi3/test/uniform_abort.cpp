// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>

#include <chrono>
#include <stdexcept>  // std::runtime_error
#include <thread>

namespace mpi3 = boost::mpi3;

// failures

void uniform_fail(mpi3::communicator const& comm) {  // cppcheck-suppress [unusedFunction,unmatchedSuppression]
	using namespace std::chrono_literals;
	std::this_thread::sleep_for(comm.rank() * 1s);

	std::cout << "uniform_fail in n = " << comm.rank() << " is about to fail" << std::endl;
	throw std::logic_error{"global but unsynchronized error"};
}

void nonuniform_fail(mpi3::communicator const& comm) {  // cppcheck-suppress [unusedFunction,unmatchedSuppression]
	using namespace std::chrono_literals;
	std::this_thread::sleep_for(comm.rank() * 1s);

	if(comm.rank() > 1) {
		std::cout << "nonuniform_fail in n = " << comm.rank() << " is about to fail" << std::endl;
		throw std::logic_error{"nonglobal error"};
	}
}

// handlers

void unconditional_abort(mpi3::communicator const& comm) {  // cppcheck-suppress unusedFunction
	std::cout << "not essential message: aborting from rank " << comm.rank() << std::endl;
	comm.abort();
}

void barriered_abort(mpi3::communicator& comm) {  // cppcheck-suppress unusedFunction
	comm.barrier();
	std::cout << "not essential message: aborting from rank " << comm.rank() << std::endl;
	comm.abort();
}

// template<class Duration>
// void abort_after(mpi3::communicator& comm, Duration d) {
// 	auto const t0 = mpi3::wall_time();
// 	while((mpi3::wall_time() - t0) < d) {  // NOLINT(altera-unroll-loops) spin loop
// 	}
// 	std::cout << "not essential message: aborting from rank " << comm.rank() << " after others join" << std::endl;
// 	comm.abort();
// }

// template<class Duration>
// void timedout_abort(mpi3::communicator& comm, Duration d) {
// 	auto       rbarrier = comm.ibarrier();
// 	auto const t0       = mpi3::wall_time();
// 	while(not rbarrier.completed() and (mpi3::wall_time() - t0) < d) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch) spin loop
// 	}

// 	if(not rbarrier.completed()) {
// 		std::cout << "non essential message: aborting from rank " << comm.rank() << " after timeout" << std::endl;
// 	} else {
// 		std::cout << "not essential message: aborting from rank " << comm.rank() << " after others join" << std::endl;
// 	}

// 	comm.abort();
// }

// template<class Duration>
// void timedout_throw(mpi3::communicator& comm, Duration d) {
// 	auto       rbarrier = comm.ibarrier();
// 	auto const t0       = mpi3::wall_time();
// 	while(not rbarrier.completed() and (mpi3::wall_time() - t0) < d) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops) spin loop
// 	}

// 	if(rbarrier.completed()) {
// 		std::cout << "non essential message: throwing from rank " << comm.rank() << " before timeout" << std::endl;
// 		throw;  // cppcheck-suppress [rethrowNoCurrentException,unmatchedSuppression] ; experimental line
// 	}
// 	std::terminate();
// }

template<class Duration>
[[noreturn]] void mpi3_timed_terminate(Duration d, mpi3::communicator& comm = mpi3::environment::get_world_instance()) {
	std::cout << "terminate called" << std::endl;
	auto       rbarrier = comm.ibarrier();
	auto const t0       = mpi3::wall_time();
	while(not rbarrier.completed() and (mpi3::wall_time() - t0) < d) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops) spin loop
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(1s);
	}

	if(rbarrier.completed()) {
		if(comm.root()) {
			std::cout << "not essential message: terminate from rank " << comm.rank() << " after others joined before timeout of " << std::chrono::duration_cast<std::chrono::seconds>(d).count() << " seconds" << std::endl;
			comm.abort(911);
		}
	} else {
		std::cout << "non essential message: terminate from rank " << comm.rank() << " after timeout of " << std::chrono::duration_cast<std::chrono::seconds>(d).count() << " seconds, not all processes failed within that time." << std::endl;
	}

	comm.abort(911);
	// never call std::terminate from here
	std::abort();  // necessary to avoid error for returning in a [[noreturn]] function
}

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int {  // NOLINT(bugprone-exception-escape) part of the test is that there is no `try` here
	assert(world.size() == 4);

// unconditional abort
#if 0
	// (-) prints only one message, (+) program terminates immediately
	try {
		uniform_fail(world);
	} catch(std::logic_error&) {
		unconditional_abort(world);
	}
#endif

#if 0
	// (-) prints only one message, (+) program terminates immediately
	try {
		nonuniform_fail(world);  // non-uniform error
	} catch(std::logic_error& e) {
		unconditional_abort(world);
	}
#endif

// barriered abort
#if 0
	// (+) prints all available messages, (+) program terminates immediately
	try {
		uniform_fail(world);
	} catch(std::logic_error& e) {
		barriered_abort(world);
	}
#endif

#if 0
	// (+) prints all available messages, (-) it DEADLOCKS (here or later)
	try {
		nonuniform_fail(world);
	} catch(std::logic_error& e) {
		barriered_abort(world);
	}
#endif

// abort after hard sleep
#if 0
	// (+) prints all available messages, (~) program terminates after hard timeout
	try {
		uniform_fail(world);  // non-uniform error
	} catch(std::logic_error&) {
	    using namespace std::chrono_literals;
		abort_after(world, 20s);
	}
#endif

#if 0
	// (+) prints all available messages, (~) program terminates after hard timeout
	try {
		nonuniform_fail(world);  // non-uniform error
	} catch(std::logic_error&) {
	    using namespace std::chrono_literals;
		abort_after(world, 20s);
	}
#endif

// timedout_abort
#if 0
	// (+) prints all available messages, (+) program terminates very quickly
	try {
		uniform_fail(world);
	} catch(std::logic_error&) {
	    using namespace std::chrono_literals;
		timedout_abort(world, 20s);
	}
#endif

#if 0
	// (+) prints all available messages, (~) program terminates after timeout
	try {
		nonuniform_fail(world);
	} catch(std::logic_error&) {
	    using namespace std::chrono_literals;
		timedout_abort(world, 20s);
	}
#endif

// timedout_terminate
#if 1
	// (+) prints all available messages, (+) program terminates very quickly
	std::set_terminate([] {
		using namespace std::chrono_literals;
		mpi3_timed_terminate(20s);
	});

	try {
		uniform_fail(world);
	} catch(std::logic_error&) {
		throw;
	}
#endif

#if 0
	// (+) prints all available messages, (~) program terminates after timeout
	std::set_terminate([]{
	    using namespace std::chrono_literals;
		mpi3_timed_terminate(20s);
	});

	try {
		nonuniform_fail(world);
	} catch(std::logic_error&) {
		std::terminate();
	}
#endif

	// I am putting a collective here to produce an deadlock if some abort strategy leaks processes
	{
		int n     = 1;
		int total = 0;
		world.all_reduce_n(&n, 1, &total);
		assert(total == world.size());
	}

	return 0;
}
// catch(...) {return 1;}
