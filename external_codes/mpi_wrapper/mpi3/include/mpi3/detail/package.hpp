// Copyright 2018-2023 Alfredo A. Correa

#ifndef BMPI3_DETAIL_PACKAGE_HPP
#define BMPI3_DETAIL_PACKAGE_HPP

// TODO(correaa) move from detail to root

#include "../../mpi3/vector.hpp"

#include "../../mpi3/detail/basic_communicator.hpp"
#include "../../mpi3/detail/buffer.hpp"
#include "../../mpi3/detail/iterator.hpp"
#include "../../mpi3/detail/iterator_traits.hpp"
#include "../../mpi3/detail/value_traits.hpp"

namespace boost {
namespace mpi3 {

class communicator;

namespace detail {

struct package : buffer {
 private:
	basic_communicator& bcomm_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members) TODO(correaa) reevaluate if a reference here is the right thing

 public:
	explicit package(communicator& comm, buffer::size_type n = 0)
	: buffer(n), bcomm_{reinterpret_cast<basic_communicator&>(comm)} {}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) TODO(correaa) break cyclic dependency of classes
	package(package const&) = delete;
	package(package     &&) = delete;

	package& operator=(package const&) = delete;
	package& operator=(package     &&) = delete;

	~package() noexcept = default;

	template<class It, typename Size>
	void pack_n(It first, Size count){
		bcomm_.pack_n(first, count, static_cast<buffer&>(*this));
	}
	template<class It>
	auto pack(It first, It last){
		bcomm_.pack(first, last, static_cast<buffer&>(*this));
	}
	template<class It, class Size>
	void unpack_n(It first, Size count){
		bcomm_.unpack_n(static_cast<buffer&>(*this), first, count);
	}
	template<class It>
	void unpack(It first, It last){
		bcomm_.unpack(static_cast<buffer&>(*this), first, last);
	}
	explicit operator bool() const {return pos < static_cast<std::ptrdiff_t>(size());}

	template<class T>
	package& operator>>(T& t){
		unpack_n(std::addressof(t), 1);
		return *this;
	}
	auto send(int dest, int tag = 0) {
		return bcomm_.send(static_cast<buffer&>(*this), dest, tag);
	}
	#if not defined(EXAMPI)
	auto receive(int dest, int tag = 0) {
		return bcomm_.receive(static_cast<buffer&>(*this), dest, tag);
	}
	#endif
};

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost
#endif
