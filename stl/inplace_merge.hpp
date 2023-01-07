#if COMPILATION_INSTRUCTIONS
(echo "#include<"$0">" > $0x.cpp) && mpicxx -O3 -std=c++17 -lboost_timer -Wfatal-errors -D_TEST_BOOST_MPI3_STL_INPLACE_MERGE $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_STL_INPLACE_MERGE_HPP
#define BOOST_MPI3_STL_INPLACE_MERGE_HPP

#include "../../mpi3/shared_window.hpp"
#include "../../mpi3/communicator.hpp"
#include<algorithm>
#include "../../mpi3/mutex.hpp"
namespace boost{
namespace mpi3{

//boost::mpi3::mutex m(world);

template<class RandomAccess>
void inplace_merge(
	communicator const& comm, 
	RandomAccess first, RandomAccess middle, RandomAccess last
){
	boost::mpi3::mutex m(communicator::world);
	m.lock();
	std::cout << "in comm " << comm.name() << "-" << comm.rank() << " of " << comm.size() << " mergin " << std::distance(first, middle) << "|" << std::distance(middle, last) << "=" << std::distance(first, last) << std::endl;
	m.unlock();
	if (comm.size() == 1) return std::inplace_merge(first, middle, last);
	if (first == middle or middle == last) return;
	if (last - first == 2) {
		using std::iter_swap;
		if (*middle < *first){
			if(comm.rank() == 0){iter_swap(first, middle);}
		}
		return;
	}
	auto first_cut = first;
	auto second_cut = middle;
	using std::lower_bound;
	if(middle - first > last - middle){
		first_cut += (middle - first)/2;
		second_cut = lower_bound(middle, last, *first_cut);
	}else{
		second_cut += (last - middle)/2;
		first_cut = lower_bound(first, middle, *second_cut);
	}
	using std::rotate;

	RandomAccess new_middle;
	int new_middle_int;
	if(comm.rank() == 0){
		new_middle = rotate(first_cut, middle, second_cut);
		new_middle_int = std::distance(first, new_middle);
	}
	comm.broadcast_value(new_middle_int);
	new_middle = first + new_middle_int;

	using boost::mpi3::inplace_merge;
	communicator newcomm(comm/2);
	if(comm.rank() <= comm.size()/2) inplace_merge(newcomm, first, first_cut, new_middle);
	else inplace_merge(newcomm, new_middle, second_cut, last);
	comm.barrier();
}

}}

static int level = 0;



#ifdef _TEST_BOOST_MPI3_STL_INPLACE_MERGE
#include<iostream>
#include<algorithm>
#include<random>

#include "alf/boost/mpi3/main.hpp"
#include "alf/boost/mpi3/mutex.hpp"
#include "alf/boost/mpi3/shm/vector.hpp"
#include<cassert>
#include <boost/timer/timer.hpp>
#include "alf/boost/mpi3/wall_clock.hpp"

int rand(int lower, int upper){
//	static std::random_device rd;
	static std::mt19937 rng;//(rd());
	static std::uniform_int_distribution<int> uni(lower, upper); 
	return uni(rng);
}
int rand(int upper = RAND_MAX){return rand(0, upper);}

namespace mpi3 = boost::mpi3;

int mpi3::main(int argc, char* argv[], boost::mpi3::communicator const& world){

	int N = 640000000;
	int MIDDLE = N/3 + 4;
	boost::mpi3::shm::vector<double> v(N, world);
	world.barrier();
	if(world.rank() == 0){
		double acc = 0;
		int i = 0;
		for(; i != MIDDLE; ++i){
			acc += rand(1, 10);
			v[i] = acc;
		}
		acc = 0;
		for(; i != N; ++i){
			acc += rand(1, 10);
			v[i] = acc;
		}
		for(int j = 0; j != N; ++j){
		//	if(j == MIDDLE) std::cout << "| ";
		//	std::cout << v[j] << " ";
		}
		std::vector<double> vv(N);
	//	std::copy(v.begin(), v.end(), vv.begin());
		{
			boost::timer::auto_cpu_timer t(std::cerr, "std::inplace_merge %w seconds\n");
			std_inplace_merge(vv.begin(), vv.begin() + MIDDLE, vv.end());
			assert(std::is_sorted(vv.begin(), vv.end()));
		}
	}
	world.barrier();
	return 0;
//	auto t1 = boost::mpi3::wall_clock::now();
	boost::mpi3::inplace_merge(world, v.begin(), v.begin() + MIDDLE, v.end());
	world.barrier();
//	auto t2 = boost::mpi3::wall_clock::now();
//	if(world.master()) std::cout << "time " << t2 - t1 << " seconds " << std::endl;
	if(world.rank() == 0){
		assert(std::is_sorted(v.begin(), v.end()));
	}
}

#endif
#endif

