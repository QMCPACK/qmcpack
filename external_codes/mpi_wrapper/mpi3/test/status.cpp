// Copyright 2019-2024 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/request.hpp>

#include <cstddef>
#include <vector>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	std::vector<std::size_t> bufsizes = {1, 100, 10000};
	mpi3::communicator&      comm     = world;

	int const dest = comm.size() - 1;
	for(std::vector<int>::size_type cs = 0; cs != bufsizes.size(); ++cs) {  // NOLINT(altera-unroll-loops) TODO(correaa) use algorithm
		if(comm.rank() == 0) {
			auto const n = bufsizes[cs];
			comm.send_n(&n, 1, dest);
			std::vector<char> buf(n);

			mpi3::request const req = comm.isend(buf.begin(), buf.end(), dest, static_cast<int>(cs + n + 1UL));
			//  req.cancel();
			//  mpi3::status s =
			//  req.wait();
			//  if(not s.cancelled()) cout << "failed to cancel request\n";
			//  else
			//  n = 0;
		} else if(comm.rank() == dest) {
			std::size_t n;  // NOLINT(cppcoreguidelines-init-variables) delayed init
			comm.receive_n(&n, 1, 0);
			if(n > 0) {
				std::vector<char> btemp(n);
				comm.receive_n(btemp.data(), n, 0);
			}
		}
		comm.barrier();
	}
	return 0;
}
