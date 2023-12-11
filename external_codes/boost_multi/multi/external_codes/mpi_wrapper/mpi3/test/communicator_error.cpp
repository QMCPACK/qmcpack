#if COMPILATION_INSTRUCTIONS
mpic++ -O3 -std=c++14 -Wall -Wextra `#-Wfatal-errors` $0 -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::error_code ec = mpi3::error::invalid_buffer_pointer;
	assert( ec == mpi3::error::invalid_buffer_pointer);
	assert( ec != std::io_errc::stream );	
	
	try{
		world.broadcast_n((int*)nullptr, 0, -1);
	}catch(std::system_error const& e){
		if(world.root())
			cout 
				<<"code --> "<< e.code() <<'\n'
				<<"message --> "<< e.code().message() <<'\n'
				<<"what --> "<< e.what() <<'\n'
			;
	}
	cout <<"finish\n";
	
	return 0;

}

