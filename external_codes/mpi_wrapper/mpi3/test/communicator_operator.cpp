#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

#include <chrono> //literals
#include <thread> //sleep_for

namespace mpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals;


auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world) -> int try {

	std::vector<double> inbuf(100);
	std::vector<double> outbuf(100);

	switch(world.rank()) {
		case 0: {
			iota(begin(outbuf), end(outbuf), 0.0);
			std::this_thread::sleep_for(2s);
			cout <<"world["<< world.rank() <<"] about to isent"<< std::endl;
			mpi3::request r = world.isend(begin(outbuf), end(outbuf), 1);
			cout <<"comm["<< world.rank() <<"] isent"<< std::endl;
		//	r.wait();
		}; break;
		case 1: {
			cout <<"comm["<< world.rank() <<"] about to ireceive"<< std::endl;
			mpi3::request r;//= world.ireceive_n(inbuf.begin(), inbuf.size(), 0);
			MPI_Irecv(
				inbuf.data(), static_cast<int>(inbuf.size()),
				detail::basic_datatype<double>{},
				MPI_ANY_SOURCE, MPI_ANY_TAG, world.get(), &r.impl_
			);
			cout <<"comm["<< world.rank() <<"] ireceived"<< std::endl;
			MPI_Wait(&r.impl_, MPI_STATUS_IGNORE);  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast) for macro
		//	r.wait();
		}; break;
		default: break;
	}
	cout <<"comm["<< world.rank() <<"] completed op"<< std::endl;

	if(world.rank() == 1) {assert( inbuf[9] == 9. );}

	return 0;
} catch(...) {return 1;}
