#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_ADDRESS -lboost_serialization $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_ADDRESS_HPP
#define BOOST_MPI3_ADDRESS_HPP

#include "../mpi3/detail/call.hpp"

#define OMPI_SKIP_MPICXX 1 // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<stdexcept>

namespace boost{
namespace mpi3{

using address = MPI_Aint;
using size_t = MPI_Aint;

inline address get_address(void const* location){
	address ret; 
	// this function requires an initialized environment, TODO should be a (static?) member of environment?
	int status = MPI_Get_address(location, &ret); // MPI_Address is deprecated
	if(status != MPI_SUCCESS) throw std::runtime_error{"error taking address"};
	return ret;
}

template<class T>
address addressof(T const& t){return get_address(std::addressof(t));}

}}

#ifdef _TEST_BOOST_MPI3_ADDRESS

#include "../mpi3/main.hpp"

using std::cout;
namespace mpi3 = boost::mpi3;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<int> v(10);
	mpi3::address a1 = mpi3::addressof(v[0]);
	mpi3::address a2 = mpi3::addressof(v[1]);
	assert( a2 - a1 == sizeof(int) );

	return 0;
}

#endif
#endif

