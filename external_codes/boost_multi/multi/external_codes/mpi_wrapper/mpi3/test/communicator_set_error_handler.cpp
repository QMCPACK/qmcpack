#include "../../mpi3/main.hpp"
#include "../../mpi3/communicator.hpp"
#include "../../mpi3/error_handler.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

auto mpi3::main(int/*argc*/, char**/*argv*/, mpi3::communicator world)->int try{

	world.set_error_handler(mpi3::error_handler::code); // default, internal function returns codes
	double d = 5.;
	try{
		world.send_n(&d, 1, 100);
	}catch(...){
		cout << "catched exception" << std::endl;
		return 0;
	}

//	world.set_error_handler(mpi3::error_handler::fatal); // fail immediately 
//	world.send(&d, &d + 1, 100);

	return 1;
} catch(...) {return 911;}
