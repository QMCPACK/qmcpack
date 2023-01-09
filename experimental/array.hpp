#if COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_WINDOW $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_ARRAY_HPP
#define BOOST_MPI3_ARRAY_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>
#include "../mpi3/window.hpp"
#include "../mpi3/mutex.cpp"

namespace boost{
namespace mpi3{

template<class T>
class array{
	mpi3::window<> win_;
	mpi3::mutex m_; //mpi3::window lock_;
	int dim1_;
	int dim2_;
	int chunk2_;
	array(mpi3::communicator& comm, int dim1, int dim2) : dim1(dim1), dim2(dim2), m_(comm){
		chunk2_ = dim2_/comm.size();
		assert( dim2_ % comm.size() == 0 );
		mpi3::size_t local_size = dim1_*chunk2_;
		win_ = comm.make_window<T>(local_size);
	}
	void set(int i, int j, T const& t){
		int rank = (j - 1)/chunk2_;
		int jfirst = rank*chunk2_+1;
		int jlast = (rank+1)*chunk2_;
		if(jlast > j) jlast = j;
		
		win_.
	}
	~array(){
		T* ptr = win_.base();
		if(ptr) mpi3::free(ptr);
	}
};

}}

#ifdef _TEST_BOOST_MPI3_ARRAY

#include "../../main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){
	return 0;
}

#endif
#endif

