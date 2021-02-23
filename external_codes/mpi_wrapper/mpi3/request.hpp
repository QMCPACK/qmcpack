#if COMPILATION_INSTRUCTIONS /* -*- indent-tabs-mode: t -*- */
(echo "#include\""$0"\"" > $0x.cpp) && mpic++ -O3 -std=c++14 -Wall `#-Wfatal-errors` -D_TEST_BOOST_MPI3_REQUEST $0x.cpp -o $0x.x && time mpirun -np 4 $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MPI3_REQUEST_HPP
#define BOOST_MPI3_REQUEST_HPP

//#include "../mpi3/communicator.hpp"
#include "../mpi3/detail/iterator.hpp" // detail::data
#include "../mpi3/status.hpp"

#define OMPI_SKIP_MPICXX 1  // https://github.com/open-mpi/ompi/issues/5157
#include<mpi.h>

#include<vector>
#include<stdexcept>

namespace boost{
namespace mpi3{

struct request{
	MPI_Request impl_ = MPI_REQUEST_NULL;
	request() = default;
	request(request const& other) = delete;// : impl_(other.impl_), owner_(false){}
//private:
//	template<
//		class ContIt,
//		class value_type = typename std::iterator_traits<ContIt>::value_type,
//		class datatype = detail::datatype<value_type>
//	>
public:
	request(request&& other) : impl_(other.impl_){other.impl_ = MPI_REQUEST_NULL;}// = default;
	request& operator=(request const&) = delete;
	request& operator=(request&& other){
		request(std::move(other)).swap(*this);
		return *this;
	}
	bool completed() const{
		int ret = -1;
		MPI_Request_get_status(impl_, &ret, MPI_STATUS_IGNORE);
		return ret;
	}
	status get_status() const{
		status ret;
		int ignore = -1;
		MPI_Request_get_status(impl_, &ignore, &ret.impl_);
		return ret;
	}
	void swap(request& other){std::swap(impl_, other.impl_);}
	void cancel(){MPI_Cancel(&impl_);}
	bool valid() const{return impl_ != MPI_REQUEST_NULL;}
	~request(){
		wait();
		if(impl_ != MPI_REQUEST_NULL) MPI_Request_free(&impl_);
	}
	void wait(){
		int s = MPI_Wait(&impl_, MPI_STATUS_IGNORE);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot wait on request");
	}
	status get(){
		status ret;
		int s = MPI_Wait(&impl_, &ret.impl_);
		if(s != MPI_SUCCESS) throw std::runtime_error("cannot wait on request");
		return ret;
	}
	void start(){
		int status = MPI_Start(&impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot start request");
	}
	status test() const{return get_status();}
};

inline std::vector<status> test_some(std::vector<request> const& requests){
	int outcount = -1;
	std::vector<int> ignore(requests.size());
	std::vector<status> ret(requests.size()); 
	int s = MPI_Testsome(requests.size(), const_cast<MPI_Request*>(&(requests.data()->impl_)), &outcount, ignore.data(), &(ret.data()->impl_));
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot test some");
	return ret;
}

inline std::vector<int> completed_some(std::vector<request> const& requests){
	int outcount = -1;
	std::vector<int> ret(requests.size());
	int s = MPI_Testsome(requests.size(), const_cast<MPI_Request*>(&(requests.data()->impl_)), &outcount, ret.data(), MPI_STATUSES_IGNORE);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot completed some");
	ret.resize(outcount);
	return ret;
}

#if 0
struct send_request : request{
	using request::request;
	template<
		class ContIt,
		class value_type = typename std::iterator_traits<ContIt>::value_type,
		class datatype = detail::datatype<value_type>
	>
	send_request(ContIt first, ContIt last, int dest, int tag, communicator& comm){
		int status = MPI_Send_init(boost::mpi3::detail::data(first), std::distance(first, last), datatype{}, dest, tag, comm.impl_, &impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot init send");
	}
};

struct receive_request : request{
	using request::request;
	template<
		class ContIt,
		class value_type = typename std::iterator_traits<ContIt>::value_type,
		class datatype = detail::datatype<value_type>
	>
	receive_request(ContIt first, ContIt last, int source, int tag, communicator& comm){
		int status = MPI_Recv_init(boost::mpi3::detail::data(first), std::distance(first, last), datatype{}, source, tag, comm.impl_, &impl_);
		if(status != MPI_SUCCESS) throw std::runtime_error("cannot init send");
	}
};
#endif

template<class ContRequestIterator, class Size>
void wait_all_n(ContRequestIterator it, Size n){
	MPI_Waitall(n, &detail::data(it)->impl_, MPI_STATUSES_IGNORE);
}

template<class ContRequestIterator>
void wait_all(ContRequestIterator it1, ContRequestIterator it2){
	wait_all_n(it1, std::distance(it1, it2));
}

//auto wait(request& r){return r.wait();}

//MPI_Request move_impl(request&& r){
//	MPI_Request ret = r.impl_;
//	r.impl_ = MPI_REQUEST_NULL;
//	return ret;
//}

template<class... Args>
void wait(Args&&... args){
	auto move_impl = [](request&& r)->MPI_Request{	MPI_Request ret = r.impl_;
		r.impl_ = MPI_REQUEST_NULL;
		return ret;
	};
	std::vector<MPI_Request> v{move_impl(std::move(args))...};
	MPI_Waitall(v.size(), v.data(), MPI_STATUSES_IGNORE);
}

template<class ContiguousIterator, class Size>
ContiguousIterator wait_any_n(ContiguousIterator it, Size n){
	int index = -1;
	int s = MPI_Waitany(n, &detail::data(it)->impl_, &index, MPI_STATUS_IGNORE);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot wait any");
	return it + index;
}

template<class ContiguousIterator>
ContiguousIterator wait_any(ContiguousIterator first, ContiguousIterator last){
	return wait_any_n(first, std::distance(first, last));
}

template<class ContiguousIterator, class Size>
std::vector<int> wait_some_n(ContiguousIterator it, Size n){
	int outcount = -1;
	std::vector<int> indices(n);
	int s = MPI_Waitsome(n, &detail::data(it)->impl_, &outcount, indices.data(), MPI_STATUSES_IGNORE);
	if(s != MPI_SUCCESS) throw std::runtime_error("cannot wait some");
	indices.resize(outcount);
	return indices;
}

template<class ContiguousIterator>
std::vector<int> wait_some(ContiguousIterator first, ContiguousIterator last){
	return wait_some_n(first, std::distance(first, last));
}

#if 0
template<class ContiguousIterator, class Size, typename value_type, class datatype>
request communicator::ireceive_n(ContiguousIterator I, Size count, int source, int tag){
	request ireceive_n;
	MPI_Irecv(mpi3::detail::data(I), count, detail::datatype<value_type>::value, source, tag, impl_, &ireceive_n.impl_);
	return ireceive_n;
}

template<class ContiguousIterator, typename value_type, class datatype>
request communicator::ireceive(ContiguousIterator first, ContiguousIterator last, int source, int tag){
	return ireceive_n(first, std::distance(first, last), source, tag);
}

template<class ContIt>
send_request communicator::send_init(ContIt first, ContIt last, int dest, int tag){
	return send_request(first, last, dest, tag, *this);
}
template<class ContIt, class Size>
send_request communicator::send_init_n(ContIt first, Size n, int dest, int tag){
	return send_init(first, first + n, dest, tag);
}
template<class ContIt>
receive_request communicator::receive_init(ContIt first, ContIt last, int source, int tag){
	return receive_request(first, last, source, tag, *this);
}
template<class ContIt, class Size>
receive_request communicator::receive_init_n(ContIt first, Size n, int source, int tag){
	return receive_init(first, first + n, source, tag);
}

template<
	class ContiguousIterator, class Size, 
	class value_type, class datatype
>
request communicator::isend_n(ContiguousIterator I, Size count, int dest, int tag){
	request ret;
	MPI_Isend(mpi3::detail::data(I), count, datatype{}, dest, tag, impl_, &ret.impl_);
	return ret;
}
template<class ContiguousIterator>
request communicator::isend(ContiguousIterator first, ContiguousIterator last, int dest, int tag){
	return isend_n(first, std::distance(first, last), dest, tag);
}
#endif

}}

#ifdef _TEST_BOOST_MPI3_REQUEST

#include "../mpi3/environment.hpp"
#include<iostream>

using std::cout;
namespace mpi3 = boost::mpi3;

int main(int argc, char* argv[]){
	mpi3::environment env(argc, argv);
	std::vector<int> buf(10);

#if 0
//	mpi3::send_
	mpi3::request r = env.world().send_init_n(buf.begin(), buf.size(), 0);

	std::vector<int> rbuf(10);
	if(env.world().rank() == 0){
		std::vector<mpi3::request> rr;//(env.world().size());
		for(int i = 0; i != env.world().size(); ++i)
			rr.emplace_back(env.world().ireceive(rbuf.begin(), rbuf.end(), i));
		r.start();
		r.wait();
		wait_all(rr.begin(), rr.end());
	}else{
		r.start();
		r.wait();
	}
#endif

#if 0
	if(env.world().rank() == 0){
	//	mpi3::receive_
		mpi3::request r = env.world().receive_init(rbuf.begin(), rbuf.end());
		mpi3::request sr = env.world().isend(buf.begin(), buf.end(), 0);
		for(int i = 0; i != env.world().size(); ++i){
			r.start();
			r.wait();
		}
		sr.wait();
	}else{
		env.world().send(buf.begin(), buf.end(), 0);
	}
#endif

	return 0;
}
#endif
#endif

