/* -*- indent-tabs-mode: t -*- */
//#if COMPILATION_INSTRUCTIONS 
//(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall `#-Wfatal-errors` -D_TEST_BOOST_MPI3_GENERALIZED_REQUEST $0x.cpp -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
//#endif
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MPI3_GENERALIZED_REQUEST_HPP
#define BOOST_MPI3_GENERALIZED_REQUEST_HPP

#include "../mpi3/request.hpp"
#include "../mpi3/status.hpp"

#include<mpi.h>

#include<stdexcept>

namespace boost {
namespace mpi3 {

struct default_request {
	static status query() {
		status ret{};
		ret.set_source(MPI_UNDEFINED);
		ret.set_tag(MPI_UNDEFINED);
		ret.set_cancelled();
		ret.set_elements<char>(0);
		return ret;
	}
	void free() {}
	void cancel(int /*complete*/) {}
};

struct generalized_request : request{
	template<class F>
	static int query_fn(void *extra_state, MPI_Status *status) {
		try {
			*status = reinterpret_cast<F*>(extra_state)->query().impl_;  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa)
		} catch(...) {
			return MPI_ERR_UNKNOWN;
		}
		return MPI_SUCCESS;
	}
	template<class F>
	static int free_fn(void* extra_state) {
		try {
			reinterpret_cast<F*>(extra_state)->free();  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa)
		} catch(...) {return MPI_ERR_UNKNOWN;}
		return MPI_SUCCESS;
	}
	template<class F>
	static int cancel_fn(void* extra_state, int complete) {
		try {
			reinterpret_cast<F*>(extra_state)->cancel(complete);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa)
		} catch(...) {return MPI_ERR_UNKNOWN;}
		return MPI_SUCCESS;
	}
	template<class F>
	explicit generalized_request(F& f){
		int s = MPI_Grequest_start(
			&query_fn<F>, //MPI_Grequest_query_function *query_fn,
  			&free_fn<F>, //MPI_Grequest_free_function *free_fn,
  			&cancel_fn<F>, //MPI_Grequest_cancel_function *cancel_fn,
  			std::addressof(f),	//	void *extra_state,
  			&impl_ //MPI_Request *request
		);
		if(s != MPI_SUCCESS) {throw std::runtime_error("cannot create generalized request");}
	}
	void complete() { MPI_(Grequest_complete)(impl_); }
};

}  // end namespace mpi3
}  // end namespace boost

//#ifdef _TEST_BOOST_MPI3_GENERALIZED_REQUEST

//#include "../mpi3/main.hpp"
//#include<iostream>

//using std::cout;
//namespace mpi3 = boost::mpi3;

//struct custom_request{
//	int counter = 0;
//	boost::mpi3::status query(){
//		counter -= 1;
//		boost::mpi3::status ret;
//		ret.set_source(MPI_UNDEFINED);
//		ret.set_tag(MPI_UNDEFINED);
//		ret.set_cancelled();
//		ret.set_elements<char>(0);
//		return ret;
//	}
//	void free(){}
//	void cancel(int complete){}
//};

//struct print_request{
////	double* first;
////	boost::mpi3::size_type n;
////	int counter = 0;
//	boost::mpi3::status query(){
//	//	counter -= 1;
//		std::cout << "query" << std::endl;
//		boost::mpi3::status ret;
//		ret.set_source(MPI_UNDEFINED);
//		ret.set_tag(MPI_UNDEFINED);
//		ret.set_cancelled();
//		ret.set_elements<char>(0);
//		return ret;
//	}
//	void free(){
//		std::cout << "free" << std::endl;
//	}
//	void cancel(int complete){
//		std::cout << "cancel " << complete << std::endl;
//	}
//};

//int mpi3::main(int, char*[], mpi3::communicator world){

//	{
//		custom_request c{};
//		mpi3::generalized_request gr(c);
//		assert(not gr.completed());
//		gr.complete();
//		gr.wait();
//	}
//	{
//		custom_request c{1};
//		mpi3::generalized_request gr(c);
//		gr.complete();
//		gr.wait();
//		assert(c.counter == 0);
//	}
//	{
//		assert( world.size() == 2);

//	//	int right = (world.rank() + 1) % world.size();
//	//	int left = world.rank() - 1;
//	//	if(left < 0) left = world.size() - 1;

//		using T = double;
//		std::vector<T> buffer(10); std::iota(buffer.begin(), buffer.end(), 0);
//		std::vector<T> buffer2(10);
//	//	mpi3::request r1 = world.ireceive(buffer2.begin(), buffer2.end(), left, 123);
//		print_request pr{};
//		mpi3::generalized_request r1(pr);
//	//	world.send(buffer.begin(), buffer.end(), right, 123);
//		std::cout << "middle" << std::endl;
//		r1.wait();
//	//	assert( buffer == buffer2 );
//	}
//	return 0;
//}

//#endif
#endif

