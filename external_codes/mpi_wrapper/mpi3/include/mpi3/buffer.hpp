#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall -Wfatal-errors -D_TEST_BOOST_MPI3_BUFFER $0x.cpp -o $0x.x && time mpirun -n 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_BUFFER_HPP
#define BOOST_MPI3_BUFFER_HPP

#define OMPI_SKIP_MPICXX 1 // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<stdexcept>
#include<cassert>

namespace boost{
namespace mpi3{

template<class T>
void buffer_attach_n(T* data, std::size_t n){
	static_assert(sizeof(T)%sizeof(char) == 0, "");
	int status = MPI_Buffer_attach((void*)data, n*sizeof(T)/sizeof(char));
	if(status != MPI_SUCCESS) throw std::runtime_error("cannot attach buffer");
}
template<class T>
void buffer_attach(T* first, T* last){
	buffer_attach_n(first, std::distance(first, last));
}

template<class T, class Size>
void attach_n(T* data, Size n){return buffer_attach_n(data, n);}

template<class T>
void attach(T* first, T* last){return buffer_attach(first, last);}

template<class C>
void attach(C& c){return buffer_attach(c);}

std::pair<char*, int> buffer_detach(){
	char* buffer = 0;
	int size = -1;
	int s = MPI_Buffer_detach(&buffer, &size);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot buffer detach"); 
	return {buffer, size};
}
std::pair<char*, int> detach(){
	return buffer_detach();
}

template<class T = char>
struct scoped_buffer_attach{
	scoped_buffer_attach(T* data, int size){ 
		buffer_attach_n(data, size);
	}
	~scoped_buffer_attach(){
		buffer_detach();
	}
};

template<class T = char>
struct scoped_buffer{
	T* buf_;
	int size_;
	scoped_buffer(int size) : size_(size){
		buf_ = new T[size];
		buffer_attach_n(buf_, size_*sizeof(T));
	}
	~scoped_buffer(){
		std::pair<char*, int> check = buffer_detach();
		assert(check.first == (char*)(buf_));
	//	assert(check.second/sizeof(T) == size_);
		delete[] buf_;
	}
};

}}

#ifdef _TEST_BOOST_MPI3_BUFFER

#include "../mpi3/main.hpp"
#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int, char*[], mpi3::communicator world){

	std::vector<int> a(10);

	mpi3::scoped_buffer<int> buf(2000);

	for(int j = 0; j !=10; ++j){
		auto r = world.bsend_init_n(a.data(), a.size(), 0, 27 + j);
		for(int i = 0; i != 10; ++i) a[i] = (world.rank() + 10*j)*world.size() + i;
		r.start();
	//	r.wait(); // requests wait automatically on destruction (they don't start automatically)
	}
	
	std::vector<int> b(10);
	
	if(world.root())
		for(int i = 0; i != world.size(); ++i)
			for(int j = 0; j != 10; ++j){
				world.receive(b.data(), i, 27 + j);
				for(int k = 0; k != 10; ++k)
					if(b[k] != (i + 10*j)*world.size() + k) assert(0);
			}

	return 0;
}

#endif
#endif

