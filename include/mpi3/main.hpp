// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_MAIN_HPP
#define BOOST_MPI3_MAIN_HPP

#include "../mpi3/communicator.hpp"
#include "../mpi3/environment.hpp"
#include "../mpi3/exception.hpp"
#include "../mpi3/timed_terminate.hpp"

#include<chrono>

namespace boost {
namespace mpi3 {

static int main(int /*argc*/, char** /*argv*/, boost::mpi3::communicator /*world*/); // if you include this file you should define `::boost::mpi3::main`NOLINT(bugprone-exception-escape)

}  // end namespace mpi3
}  // end namespace boost

// cppcheck-suppress syntaxError ; bug cppcheck 2.3
auto main(int argc, char** argv) -> int /*try*/ {  // NOLINT(misc-definitions-in-headers,bugprone-exception-escape) : if you include this file you shouldn't have your own `::main`, you should define `boost::mpi3::main(int argc, char** argv, boost::mpi3::communicator world)` instead
	boost::mpi3::environment const env{argc, argv};
	std::set_terminate([]{
	    using namespace std::chrono_literals;
		boost::mpi3::timed_terminate(3s);
	});
//	try {
		int const ret = boost::mpi3::main(argc, argv, /*env.*/ boost::mpi3::environment::get_world_instance());
		boost::mpi3::environment::get_world_instance().barrier();
		return ret;
//	} catch(std::exception& e) {
//		if(boost::mpi3::environment::get_world_instance().root()) {std::cerr<<"exception message: "<< e.what() <<"\n\n\n"<<std::endl;}
//		return 1;
//	} catch(...) {
//		if(boost::mpi3::environment::get_world_instance().root()) {std::cerr<<"unknown exception message"<<std::endl;}
//		return 1;
//	}
}
// catch(...) {
//	std::cerr<<"unknown error in MPI pogram"<<std::endl;
//	return 1;
//}
#endif
