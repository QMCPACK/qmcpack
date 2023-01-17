/* -*- indent-tabs-mode: t -*- */

#ifndef BOOST_MPI3_VECTOR_HPP
#define BOOST_MPI3_VECTOR_HPP

#include "../mpi3/allocator.hpp"

#include<mpi.h>

#include<vector>

namespace boost{
namespace mpi3{

template<class T>
using vector = std::vector<T, mpi3::allocator<T>>;

template<class T>
using uvector = std::vector<T, mpi3::uallocator<T>>;

//template<class T>
//struct uvector_iterator : std::vector<T, mpi3::uallocator<T>>{};

/*
template<class T>
struct uvector : std::vector<T, mpi3::uallocator<T>>{
	using std::vector<T, mpi3::uallocator<T>>::vector;
	uvector(std::vector<T> const& other) 
	: std::vector<T, mpi3::uallocator<T>>(reinterpret_cast<std::vector<T, mpi3::uallocator<T>>&>(other)){}
	uvector(std::vector<T>&& other) 
	: std::vector<T, mpi3::uallocator<T>>(std::move(reinterpret_cast<std::vector<T, mpi3::uallocator<T>>&>(other))){}
	operator std::vector<T>&&() &&{
		return std::move(reinterpret_cast<std::vector<T>&>(*this));
	}
	operator std::vector<T>&() &{
		return reinterpret_cast<std::vector<T>&>(*this);
	}
	operator std::vector<T> const&() const&{
		return reinterpret_cast<std::vector<T> const&>(*this);
	}
	
};*/

}  // namespace mpi3
}  // namespace boost

//#ifdef _TEST_BOOST_MPI3_VECTOR

//#include <boost/timer/timer.hpp>

//#include "../mpi3/main.hpp"
//#include "../mpi3/detail/iterator.hpp"

//using std::cout;
//namespace mpi3 = boost::mpi3;

//int mpi3::main(int, char*[], mpi3::communicator world){

//	mpi3::vector<long long> v(100);
//	std::iota(v.begin(), v.end(), 0);
//	assert( std::accumulate(v.begin(), v.end(), 0) == (v.size()*(v.size() - 1))/2 );

//	mpi3::uvector<std::size_t> uv(100);
////	assert( std::accumulate(uv.begin(), uv.end(), 0) == 0 );

//	std::vector<std::size_t> v2(uv.begin(), uv.end());
//	assert(v2.size() == 100);
//	assert(uv.size() == 100);

//	std::vector<std::size_t> v3(uv.begin(), uv.end()); uv.clear();
//	assert(v3.size() == 100);
//	assert(uv.size() == 0);

//	using boost::timer::cpu_timer;
//	using boost::timer::cpu_times;

//	auto test_size = 100000000; // 100MB
//	cout << "test size " << test_size << " chars \n";
//	using vector = std::vector<char, mpi3::allocator<char>>;
//	cout << "std::vector<char, mpi3::allocator<char>> (pinned, initialized)\n";
//	{
//		cpu_timer t;
//		t.start();
//		cpu_times allocation_start = t.elapsed();
//		vector v(test_size);
//		cpu_times allocation_end = t.elapsed();
//		cpu_times firstuse_start = t.elapsed();
//		for(std::size_t i = 0; i != test_size; ++i) v[i] = 'l';
//		cpu_times firstuse_end = t.elapsed();
//		cpu_times seconduse_start = t.elapsed();
//		for(std::size_t i = 0; i != test_size; ++i) v[i] = 'z';
//		cpu_times seconduse_end = t.elapsed();
//		cout 
//			<< "\tallocation (+initialization if any) " << (allocation_end.wall - allocation_start.wall)/1.e6 << " milliseconds \n"
//			<< "\tfirst use " << (firstuse_end.wall - firstuse_start.wall)/1.e6 << " milliseconds \n"
//			<< "\tsecond use " << (seconduse_end.wall - seconduse_start.wall)/1.e6 << " milliseconds \n"
//		;
//	}
///*
//output:
//test size 100000000 chars (mpic++ -O3)
//std::vector<char> (no pinned, initialized)
//	allocation (+initialization if any) 33.3729 milliseconds 
//	first use 4.78287 milliseconds 
//	second use 6.60647 milliseconds 

//test size 100000000 chars 
//std::vector<char, mpi3::allocator<char>> (pinned, initialized)
//	allocation (+initialization if any) 0.006804 milliseconds 
//	first use 31.1551 milliseconds 
//	second use 4.82273 milliseconds 

//std::vector<char, mpi3::uallocator<char>> (pinned, uninitialized)
//	allocation (+initialization if any) 0.007946 milliseconds 
//	first use 30.7246 milliseconds 
//	second use 4.92651 milliseconds 
//*/
//	return 0;
//}

//#endif
#endif

