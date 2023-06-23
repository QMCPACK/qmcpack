// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_PROCESS_HPP
#define BOOST_MPI3_PROCESS_HPP

#include "../mpi3/communicator.hpp"

#include <optional>

#include "config/NODISCARD.hpp"

namespace boost {
namespace mpi3 {

using std::optional;

class process {
	communicator& comm_;
	int rank_;
	friend boost::mpi3::communicator;

	process(communicator& comm, int rank) : comm_{comm}, rank_{rank} {}

 public:
	communicator& comm() const {return comm_;}

	int rank() const{return rank_;}
	template<class T>
	optional<T> operator+=(T const& t) && {
		T val = comm_.reduce_value(t, std::plus<>{}, rank_);
		if(rank_ != comm_.rank()) {return {};}
		return optional<T>(val);
	}

	template<class T>
	process& operator>>(T& t) & {
		comm_.receive_n(&t, 1, rank_);
	//	comm_.receive_value(t, rank_);
		return *this;  // NOLINT(hicpp-move-const-arg,performance-move-const-arg) TODO(correaa)
	}
	template<class T>
	process&& operator&(T& t) && {
		comm_.broadcast_value(t, rank_);
		return std::move(*this);  // NOLINT(hicpp-move-const-arg,performance-move-const-arg) TODO(correaa)
	}
};

template<class T>
auto operator<<(process&& p, const T& value) -> decltype(std::move(p << value)) {
	return std::move(p << value);
}

template<class T>
auto operator>>(process&& p, T&& value) -> decltype(std::declval<process&>() >> value) {
	return p >> value;
}

template<class T> 
process& operator<<(process& self, T const& t) {
	self.comm().send_value(t, self.rank());
	return self;
}

inline auto communicator::operator[](int rank) -> process {return {*this, rank};}

//inline auto communicator::iterator::operator*() const -> process {return {*commP_, rank_};}

template<class T>
auto operator&(communicator& comm, T&& t)
->decltype(comm.all_to_all(begin(std::forward<T>(t))), std::forward<T>(t)) {
	assert( t.size() == comm.size() );
//	using std::begin;
	auto e = comm.all_to_all(begin(std::forward<T>(t)));
	using std::end;
	assert( e == end(t) );
	return std::forward<T>(t);
}

template<class T> 
//NODISCARD("do not ignore result when second argument is const")
auto operator&(communicator& comm, T const& t)
->decltype(comm.all_to_all(t.begin(), std::declval<T>().begin()), T(comm.size())) {
	assert(t.size() == comm.size());
	T ret(comm.size()); 
	comm.all_to_all(t.begin(), ret.begin());
	return ret;
}

template<class T>
auto operator||(process&& self, T& t) {self.comm().broadcast_value(t, self.rank());}

template<class T>
communicator& operator>>(communicator& comm, T& t) {
	comm.receive_n(&t, 1);
//	comm.receive_value(t);
	return comm;
}
template<class T>
std::vector<T> operator|=(communicator& comm, T const& t) {
	return comm.all_gather_value(t);
}

template<class T>
std::vector<T> operator|=(process&& self, T const& t) {
	return self.comm().gather_value(t, self.rank());
}

template<class T>
std::pair<T, process> communicator::max_location(T const& t) {
	auto const ml = max_loc(t);
	return std::pair<T, process>{
		ml.value,
		process{*this, ml.location}
	};
}

}  // end namespace mpi3
}  // end namespace boost
#endif
