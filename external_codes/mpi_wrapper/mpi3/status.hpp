#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 `#-Wfatal-errors` -D_TEST_BOOST_MPI3_STATUS $0x.cpp -o $0x.x && time mpirun -n 2 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_STATUS_HPP
#define BOOST_MPI3_STATUS_HPP

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>
#include "../mpi3/detail/datatype.hpp"

#include<stdexcept>

namespace boost{
namespace mpi3{

struct status{
	MPI_Status impl_;
	status() = default;
	~status(){
	//	if(impl_ != MPI_STATUS_NULL) 
	//	MPI_Status_free(&impl_);
	}
	template<class T>// = char>
	int count() const{ //entries received of datatype T
		int ret = -1;
		int s = MPI_Get_count(&impl_, detail::basic_datatype<T>{}, &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot count");
		return ret;
	}
	template<class T>
	int elements() const{
		int ret = -1;
		int s = MPI_Get_elements(&impl_, detail::basic_datatype<T>{}, &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot elements");
		return ret;
	}
	template<class T>
	void set_elements(int count){
		int s = MPI_Status_set_elements(&impl_, detail::basic_datatype<T>{}, count);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot set elements");
	}
	int source() const{return impl_.MPI_SOURCE;}
	void set_source(int s){impl_.MPI_SOURCE = s;}
//	void source(int s){set_source(s);}

	int tag() const{return impl_.MPI_TAG;}
	void set_tag(int t){impl_.MPI_TAG = t;}
//	void tag(int t){set_tag(t);}
	void set_cancelled(bool flag = true){
		int s = MPI_Status_set_cancelled(&impl_, flag);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot set cancelled");
	}
	bool cancelled() const{
		int ret = -1;
		int s = MPI_Test_cancelled(&impl_, &ret);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot test cancelled");
		return ret;
	}
};

}}

#ifdef _TEST_BOOST_MPI3_STATUS

#include "../mpi3/main.hpp"
#include "../mpi3/request.hpp"

#include<iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int mpi3::main(int argc, char* argv[], mpi3::communicator& world){
	std::vector<int> bufsizes = {1, 100, 10000, 1000000};
	mpi3::communicator& comm = world;

	int dest = comm.size() - 1;
	for(int cs = 0; cs != bufsizes.size(); ++cs){
		if(comm.rank() == 0){
			int n = bufsizes[cs];
			std::vector<char> buf(n);
			mpi3::request req = comm.isend(buf.begin(), buf.end(), dest, cs + n + 1);
			req.cancel();
		//	mpi3::status s = 
			req.wait();
		//	if(not s.cancelled()) cout << "failed to cancel request\n";
		//	else 
			n = 0;
			comm.send_n(&n, 1, dest, 123);
			n = cs + n + 1;
			comm.send_n(&n, 1, dest, 123);
		}else if(comm.rank() == dest){
			int n = -1;
			int tag = 0;
			comm.receive_n(&n, 1, 0, 123);
			comm.receive_n(&tag, 1, 0, 123);
			if(n > 0){
				std::vector<char> btemp(n);
				comm.receive_n(btemp.data(), n, 0, tag); 
			}
		}
		comm.barrier();
	}
}

#endif
#endif

