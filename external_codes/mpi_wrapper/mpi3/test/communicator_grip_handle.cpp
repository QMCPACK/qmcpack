// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"

#include<iostream>
#include<numeric>
#include<vector>

namespace bmpi3 = boost::mpi3;

int main(int argc, char** argv) try {

	MPI_Init(&argc, &argv);
	MPI_Comm W{};
	MPI_Comm_dup(MPI_COMM_WORLD, &W);
	{
		bmpi3::communicator& w = bmpi3::grip_communicator(W);
		// cppcheck-suppress[assertWithSideEffect,unmatchedSuppression]
		assert(w.handle() == W);

		std::vector<double> const xsend(10,  5.);
		std::vector<double>       xrecv(10, -1.);

		switch( w.rank() ) {
			case 0: {
				w.receive(begin(xrecv), end(xrecv), 1);
				assert(xrecv[5] ==  5.);
				break;
			}
			case 1: {
				w.send(begin(xsend), end(xsend), 0);
				assert(xrecv[5] == -1.);
				break;
			}
		}
	}
	MPI_Comm_free(&W);
	MPI_Finalize();

	return 0;
} catch(...) {return 1;}
