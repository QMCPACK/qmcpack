#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -Wall -Wextra -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_WALL_CLOCK $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
//  (C) Copyright Alfredo A. Correa 2018.
#ifndef BOOST_MPI3_WALL_CLOCK_HPP
#define BOOST_MPI3_WALL_CLOCK_HPP

#define OMPI_SKIP_MPICXX 1 // workaround for https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<chrono>

namespace boost{
namespace mpi3{

inline auto wall_time(){return std::chrono::duration<double>(MPI_Wtime());}
inline auto wall_tick(){return std::chrono::duration<double>(MPI_Wtick());}

struct wall_clock{
	using rep = double;
	using period = std::ratio<1>; // one second
	using duration = std::chrono::duration<double>;
	using time_point = std::chrono::time_point<wall_clock>;
	static time_point now() noexcept{return time_point{wall_time()};};
	static duration tick() noexcept{return duration{wall_tick()};}
};

template<class Duration = std::chrono::duration<double>>
void wall_sleep_for(Duration d){
	auto then = wall_clock::now();
	while(wall_clock::now() - then < d){}
}

}}

#ifdef _TEST_BOOST_MPI3_WALL_CLOCK

#include "../mpi3/environment.hpp"

#include<iostream>
#include<thread> // this_tread::sleep_for

namespace mpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals; // 2s

int main(int argc, char* argv[]){
	mpi3::environment::initialize(argc, argv); // same as MPI_Init(...);
	assert(mpi3::environment::is_initialized());

	auto then = mpi3::wall_time();
	std::this_thread::sleep_for(2s);
	cout << (mpi3::wall_time() - then).count() <<'\n'; 

	mpi3::environment::finalize(); // same as MPI_Finalize()
	assert(mpi3::environment::is_finalized());
	return 0;
}
#endif
#endif

