/* -*- indent-tabs-mode: t -*- */

#ifndef BOOST_MPI3_ALLOCATOR_HPP
#define BOOST_MPI3_ALLOCATOR_HPP

#include "../mpi3/address.hpp"

#include<mpi.h>

#include<limits>
#include<vector>

namespace boost{
namespace mpi3{

struct /*__attribute__((aligned(0)))*/ bad_alloc : std::bad_alloc{using std::bad_alloc::bad_alloc;};

inline void* malloc(mpi3::size_t size) {
	void* ret;  // NOLINT(cppcoreguidelines-init-variables) delayed init
	int const s = MPI_Alloc_mem(size, MPI_INFO_NULL, &ret);
	if(s != MPI_SUCCESS) {return nullptr;}  //s throw bad_alloc();//"cannot allocate " + std::to_string(size) + " bytes");
	return ret;
}

inline void free(void* ptr){
	MPI_(Free_mem)(ptr);
}

template<class T>
struct /*__attribute__((aligned(0)))*/ allocator{
	using size_type = mpi3::size_t;
	using value_type = T;
	using pointer = T*;

	allocator() = default;

	// cppcheck-suppress noExplicitConstructor
	template<class U> allocator(allocator<U> const&/*other*/) {}  // NOLINT(google-explicit-constructor,hicpp-explicit-conversions) : allocator convention

	auto allocate(size_type n) {
		if(void* ptr = mpi3::malloc(n * static_cast<size_type>(sizeof(T)))) {
			return static_cast<pointer>(ptr);
		}
		throw bad_alloc();
	}
	void             deallocate(pointer p, std::size_t /*size*/) { mpi3::free(p); }
	static size_type max_size() { return std::numeric_limits<size_type>::max(); }
};

template<typename T>
struct uallocator : allocator<T>{
	template<class U> void construct(U* /*unused*/){
		static_assert(
			std::is_trivially_destructible<T>{}, 
			"uallocator cannot be used with non trivial types"
		);
	}
};

template< class T1, class T2 > constexpr
bool operator==(allocator<T1> const&/*self*/, uallocator<T2> const&/*other*/) {  // TODO(correaa) check that both(?) are trivial?
	return true;
}

template< class T1, class T2 > constexpr
bool operator==(uallocator<T1> const&/*self*/, allocator<T2> const&/*other*/) {
	return true;
}

template <class T>
constexpr std::add_const_t<T>& as_const(T& t) noexcept{return t;}

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_ALLOCATOR

//#include "../mpi3/main.hpp"

//#include<cassert>
//#include<vector>

//#include<boost/container/flat_set.hpp>

//namespace mpi3 = boost::mpi3;
//using std::cout;

//int mpi3::main(int argc, char* argv[], mpi3::communicator world){

//	std::vector<double, mpi3::allocator<double>> v(1000000);
//	std::vector<double, mpi3::uallocator<double>> uv(1000000);
//	std::iota(v.begin(), v.end(), 0.);
//	using boost::mpi3::data;
//	assert( data(uv.begin()) == &*uv.begin() );
//	assert( std::accumulate(v.begin(), v.end(), 0.) == (v.size()*(v.size() - 1))/2 );
//	return 0;
//	
//	{
//		boost::container::flat_set<double, std::less<double>, mpi3::allocator<double> > fs;
//		fs.insert(5.);
//		fs.insert(3.);
//		auto it = fs.begin(); 
//		assert(*it == 3.); 
//		++it; 
//		assert(*it == 5.); 
//	}
//	{
//		boost::container::flat_set<int, std::less<int>, std::allocator_traits<mpi3::allocator<double>>::rebind_alloc<int>> fs;
//		fs.insert(5);
//		fs.insert(3);
//		auto it = fs.begin(); 
//		assert(*it == 3); 
//		++it; 
//		assert(*it == 5);
//	} 
//	{
//		boost::container::flat_set<std::pair<double, double>, std::less<std::pair<double, double>>, mpi3::allocator<std::pair<double, double>>> fsp;
//		fsp.insert({1.,2.});
//		fsp.insert({3.,4.});
//		auto it = fsp.begin(); 
//		assert(*it == std::make_pair(1.,2.)); 
//		++it; 
//		assert(*it == std::make_pair(3.,4.)); 
//	}
//	return 0;
//}

//#endif
#endif
