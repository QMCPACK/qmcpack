// Copyright 2021-2025 Alfredo A. Correa

// Test for separate compilation / library usage.

#include "../communicator.hpp"

namespace mpi3 = boost::mpi3;

void do_broadcast(mpi3::communicator &c) {  // NOLINT(misc-use-internal-linkage) declared and used elsewhere // cppcheck-suppress unusedFunction
	int a = 2;
	c.broadcast_value(a);
}

