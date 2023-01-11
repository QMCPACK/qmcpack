/* -*- indent-tabs-mode: t -*- */
//#if COMPILATION
//mpic++ -x c++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
//#endif
// Copyright 2018-2022 Alfredo A. Correa

#ifndef BOOST_MPI3_WALL_CLOCK_HPP
#define BOOST_MPI3_WALL_CLOCK_HPP

#include<mpi.h>

#include "communicator.hpp"

#include<chrono>

namespace boost {
namespace mpi3 {

inline auto wall_time(){return std::chrono::duration<double>(MPI_Wtime());}
inline auto wall_tick(){return std::chrono::duration<double>(MPI_Wtick());}

struct wall_clock {
	using rep = double;
	using period = std::ratio<1>; // one second
	using duration = std::chrono::duration<double>;
	using time_point = std::chrono::time_point<wall_clock>;

	static time_point now() noexcept{return time_point{wall_time()};};
	static duration tick() noexcept{return duration{wall_tick()};}
};

template<class Duration = std::chrono::duration<double>>
void wall_sleep_for(Duration d) {
	auto then = wall_clock::now();
	// spin now
	while(wall_clock::now() - then < d) {}  // NOLINT(altera-unroll-loops) this is time loop
}

class wall_timer {
	mpi3::communicator comm_;
	std::string title_;
	wall_clock::time_point start_;

 public:
	explicit wall_timer(mpi3::communicator comm, std::string title = "")
	: comm_{std::move(comm)}, title_{std::move(title)}, start_{wall_clock::now()} {}

	wall_timer(wall_timer const&) = delete;
	wall_timer(wall_timer     &&) = delete;

	wall_timer& operator=(wall_timer const&) = delete;
	wall_timer& operator=(wall_timer     &&) = delete;

	~wall_timer() {  // NOLINT(bugprone-exception-escape) TODO(correaa) may be it should be able to throw
		auto const diff = wall_clock::now() - start_;
		auto const min = comm_.min(diff.count());  // cppcheck-suppress unreadVariable ; bug in cppcheck 2.3
		auto const max = comm_.max(diff.count());  // cppcheck-suppress unreadVariable ; bug in cppcheck 2.3
		auto const total = (comm_ += diff.count());
		auto const avg = total/comm_.size();
		auto const speed_up = max / total;
		if(comm_.root()) {std::cerr<<"# "<< title_ <<" timing "<< min <<"["<< avg <<"]"<< max <<" sec, speed up = x"<< speed_up <<std::endl;}
	}
};

}  // end namespace mpi3
}  // end namespace boost

//#if not __INCLUDE_LEVEL__

//#include "../mpi3/environment.hpp"

//#include<iostream>
//#include<thread> // this_tread::sleep_for

//namespace bmpi3 = boost::mpi3;
//using std::cout;
//using namespace std::chrono_literals; // 2s

//void f(bmpi3::communicator c){
//	bmpi3::wall_timer f_watch{c, "OptionalTitle"};
//	std::this_thread::sleep_for(std::chrono::seconds(c.rank()+1));
//// prints "# OptionalTitle timing 1.00007[2.50007]4.00006 sec"
//}

//int main(int argc, char* argv[]){
//	bmpi3::environment::initialize(argc, argv); // same as MPI_Init(...);
//	assert(bmpi3::environment::is_initialized());

//	auto then = bmpi3::wall_time();
//	std::this_thread::sleep_for(2s);
//	cout<< (bmpi3::wall_time() - then).count() <<'\n'; 

//	{
//		auto world = bmpi3::environment::get_world_instance();
//		f(world);
//	}

//	bmpi3::environment::finalize(); // same as MPI_Finalize()
//	assert(bmpi3::environment::is_finalized());
//}
//#endif
#endif

