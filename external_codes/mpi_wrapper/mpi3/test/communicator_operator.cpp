#if COMPILATION_INSTRUCTIONS
mpic++ $0 -o $0x&&time mpirun -n 2 $0x&&rm $0x;exit
#endif

#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

#include <chrono> //literals
#include <thread> //sleep_for

namespace mpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals;


int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<double> inbuf(100);
	std::vector<double> outbuf(100);

	mpi3::request r;
	switch(world.rank()){
		case 0: {
			iota(begin(outbuf), end(outbuf), 0.0);
			std::this_thread::sleep_for(2s);
			cout <<"world["<< world.rank() <<"] about to isent"<< std::endl;
			r = world.isend(begin(outbuf), end(outbuf), 1);
			cout <<"comm["<< world.rank() <<"] isent"<< std::endl;
			r.wait();
		}; break;
		case 1: {
			cout <<"comm["<< world.rank() <<"] about to ireceive"<< std::endl;
			r = world.ireceive_n(inbuf.begin(), inbuf.size(), 0);
			cout <<"comm["<< world.rank() <<"] ireceived"<< std::endl;
			r.wait();
		}
	}
	cout <<"comm["<< world.rank() <<"] completed op"<< std::endl;

	if(world.rank() == 1) assert( inbuf[9] == 9. );

	return 0;
}

