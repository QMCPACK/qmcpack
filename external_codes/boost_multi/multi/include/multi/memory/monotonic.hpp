#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef BOOST_MULTI_MEMORY_MONOTONIC_HPP
#define BOOST_MULTI_MEMORY_MONOTONIC_HPP

#include "../memory/block.hpp"
#include "../memory/allocator.hpp"

#include<cstddef> // max_align_t
#include<stdexcept>

namespace boost {
namespace multi {

// template<class T> auto raw_pointer_cast(T* p) -> T* {return p;}  // with the meaning of std::to_address

namespace memory {

template<class T>
T* align_up(T* ptr, std::size_t bytes = alignof(std::max_align_t)) {
//	return
//		reinterpret_cast<T*>(
//			  (reinterpret_cast<std::uintptr_t>(ptr) + (bytes-1))
//			& ~(bytes-1)
//		)
//	;
	using uintptr_t = std::uint64_t;
	static_assert( sizeof(uintptr_t) == sizeof(T*), "this function works in 64 bit systems" );
	return reinterpret_cast<T*>( bytes * ((reinterpret_cast<uintptr_t&>(ptr) + (bytes - 1)) / bytes) );
}

template<class Ptr>  // TODO test with actual fancy ptr
constexpr
Ptr align_up(Ptr ptr, std::size_t bytes = alignof(std::max_align_t)) {
	using multi::to_address;
	auto p_(to_address(ptr));
////  using multi::raw_pointer_cast;
////  auto p_{raw_pointer_cast(p)};

	static_assert( sizeof(*p_)==1 , "!"); // crash
//	auto q_ = reinterpret_cast<decltype(p_)>(
//		(reinterpret_cast<std::uintptr_t>(p_) + (align-1))
//		& ~(align-1)
//	);
	auto q_ = align_up(p_, bytes);
	return ptr + std::distance(p_, q_);
//	return reinterpret_cast<Ptr&>( bytes * ((reinterpret_cast<std::uintptr_t&>(ptr) + (bytes - 1)) / bytes) );  // maybe using uint64_t and static_assert sizeof(void*) == uint64_t
}

template<typename Ptr = byte*, std::size_t Align = alignof(std::max_align_t)>
class null_t {
 public:
	using pointer = Ptr;
	using size_type = typename std::pointer_traits<pointer>::difference_type;

 private:
	using void_pointer = typename std::pointer_traits<pointer>::template rebind<void>;

 public:
	template<std::size_t AA = Align>
	void_pointer allocate(size_type required_bytes, size_type = AA) {
		if(required_bytes > 0) {throw std::bad_alloc{};}
		return nullptr;
	}
	void deallocate(void_pointer p, size_type /*discarded_bytes*/) {
		if(p != nullptr) {throw std::bad_alloc{};}
	}
};

template<typename Ptr = byte*, std::size_t Align = alignof(std::max_align_t)>
class monotonic : protected block<Ptr> {
 protected:
	using void_pointer = typename std::pointer_traits<typename monotonic::pointer>::template rebind<void>;

 public:
	using block<Ptr>::start_;
	using block<Ptr>::size_;
	using block<Ptr>::block;
	typename monotonic::pointer position_ = start_;
	void reset() {position_ = this->start_;}
	template<std::size_t AA = Align>
	typename monotonic::void_pointer allocate(
		typename monotonic::size_type required_bytes,
		typename monotonic::size_type align = AA//alignof(std::max_align_t)
	) {
		auto ret = align_up(this->position_, align);
		auto new_position_ = ret + required_bytes;
		using std::distance;
		if(not this->contains(new_position_-1)) {
			throw overflow(required_bytes, this->size_ - distance(start_, position_));
		}
		this->position_ = new_position_;
		return ret;
	}
	bool owns(typename monotonic::void_pointer p) const {
		return this->contains(static_cast<typename monotonic::pointer>(p));
	}
	void deallocate(
		typename monotonic::void_pointer p, 
		typename monotonic::size_type /*discarded_bytes*/
	) {
		if(not owns(p)) {throw std::bad_alloc{};}
	}
	struct overflow : public std::bad_alloc{
		using size_type = typename monotonic::size_type;
		size_type required;
		size_type available;
	//	constexpr auto to_string = [](auto a){return std::to_string(a);};
		std::string msg;
		overflow(size_type required, size_type available)
		: required{required}, available{available},
			msg{"required "+std::to_string(required)+" while only "+std::to_string(available)+" bytes available"} {}
		virtual const char* what() const throw() override {return msg.c_str();}// + std::to_string(required)).c_str();}
	};
};

template<class T = void> 
using monotonic_allocator = multi::memory::allocator<T, monotonic<char*>>;

}  // end namespace memory
}  // end namespace multi
}  // end namespace boost

#if not __INCLUDE_LEVEL__ // _TEST_BOOST_MULTI_MEMORY_MONOTONIC

#include "../../multi/array.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
using std::cout;

int main(){
{
	multi::memory::null_t<char*> mr;
	try {
		mr.allocate(1*sizeof(double), alignof(double));
	} catch(...) {}
}
{
	alignas(double) std::array<char, 256*sizeof(double)> buffer;  // char buffer[256*sizeof(double)];
	multi::memory::monotonic<char*> m(buffer.data(), buffer.size());
	auto p1 = m.allocate(1*sizeof(double), alignof(double));
	auto p2 = m.allocate(255*sizeof(double), alignof(double));
	m.deallocate(p2, 255*sizeof(double));
	m.deallocate(p1, 1*sizeof(double));
	try {
		m.deallocate((char*)p1 + 10000, 1*sizeof(double));
	} catch(...){}
}
{
	alignas(double) std::array<char, 300*sizeof(double)> buffer;  // char buffer[300*sizeof(double)];
	multi::memory::monotonic<char*> m(buffer.data(), buffer.size());
	multi::memory::monotonic_allocator<double> alloc(&m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> A({10, 10}, &m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> B({10, 10}, &m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> C({10, 10}, &m);
}
}
#endif
#endif
