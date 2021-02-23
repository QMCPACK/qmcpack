#if COMPILATION /* -*- indent-tabs-mode: t -*- */
mpic++ -x c++ $0 -o $0x&&mpirun -n 4 $0x&&rm $0x;exit
#endif
// © Copyright Alfredo A. Correa 2018-2019s

#ifndef BOOST_MPI3_WALL_CLOCK_HPP
#define BOOST_MPI3_WALL_CLOCK_HPP

#define OMPI_SKIP_MPICXX 1 // workaround for https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include "communicator.hpp"

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

struct wall_timer{
	mpi3::communicator comm_;
	wall_clock::time_point start_;
	std::string title_;
	wall_timer(wall_timer const&) = delete;
	wall_timer(mpi3::communicator comm, std::string title = "") : comm_{std::move(comm)}, start_{wall_clock::now()}, title_{title}{}
	~wall_timer(){
		auto diff = wall_clock::now() - start_;
		auto min = comm_.min(diff.count());
		auto max = comm_.max(diff.count());
		auto avg = (comm_ += diff.count())/comm_.size();
		if(comm_.root()) std::cerr<<"# "<< title_ <<" timing "<< min <<"["<< avg <<"]"<< max <<" sec"<<std::endl;
	}
};

}}

#if not __INCLUDE_LEVEL__

#include "../mpi3/environment.hpp"

#include<iostream>
#include<thread> // this_tread::sleep_for

namespace bmpi3 = boost::mpi3;
using std::cout;
using namespace std::chrono_literals; // 2s

void f(bmpi3::communicator c){
	bmpi3::wall_timer f_watch{c, "OptionalTitle"};
	std::this_thread::sleep_for(std::chrono::seconds(c.rank()+1));
// prints "# OptionalTitle timing 1.00007[2.50007]4.00006 sec"
}

int main(int argc, char* argv[]){
	bmpi3::environment::initialize(argc, argv); // same as MPI_Init(...);
	assert(bmpi3::environment::is_initialized());

	auto then = bmpi3::wall_time();
	std::this_thread::sleep_for(2s);
	cout<< (bmpi3::wall_time() - then).count() <<'\n'; 

	{
		auto world = bmpi3::environment::get_world_instance();
		f(world);
	}

	bmpi3::environment::finalize(); // same as MPI_Finalize()
	assert(bmpi3::environment::is_finalized());
}
#endif
#endif

