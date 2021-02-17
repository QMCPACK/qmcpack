#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wfatal-errors -D_TEST_BOOST_MPI3_DYNAMIC_WINDOW $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_DYNAMIC_WINDOW_HPP
#define BOOST_MPI3_DYNAMIC_WINDOW_HPP

#include<mpi.h>
#include "../mpi3/window.hpp"

namespace boost{
namespace mpi3{

template<class T = void>
struct dynamic_window : window<T>{
	protected:
	dynamic_window() : window<T>{}{}
	public:
	dynamic_window(communicator& comm){
		int s = MPI_Win_create_dynamic(MPI_INFO_NULL, comm.get(), &(this->impl_));
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot create dynamic window");
	}
	template<class TT = char>
	void attach_n(TT* base, mpi3::size_t n){MPI_Win_attach(this->impl_, base, n*sizeof(TT));}
	void detach(void* base){MPI_Win_detach(this->impl_, base);}
/*	void* bottom() const{
		void* base; int flag;
		int s = MPI_Win_get_attr(impl_, MPI_BOTTOM, &base, &flag);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot get base");
		assert(flag);
		return base;
	}*/

};

//window communicator::make_dynamic_window(){
//	return dynamic_window(t, n, *this);
//}

}}

#ifdef _TEST_BOOST_MPI3_DYNAMIC_WINDOW

#include "../mpi3/main.hpp"

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator world){

	mpi3::dynamic_window<int> dw(world);
	std::vector<int> a(1000, world.rank());
	dw.attach_n(a.data(), a.size());
	world.barrier();
	dw.lock(0);
	if(world.rank() == 0){
		std::vector<int> v(1000);
		dw.get_n(v.data(), v.size(), 3, 0);
		assert(a[0] == 3);
	}
	dw.unlock(0);
	world.barrier();
	dw.detach(a.data());

	return 0;
}

#endif
#endif

