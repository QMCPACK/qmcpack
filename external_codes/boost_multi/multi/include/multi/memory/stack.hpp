#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef BOOST_MULTI_MEMORY_STACK_HPP
#define BOOST_MULTI_MEMORY_STACK_HPP

#include "../memory/monotonic.hpp"
#include "../memory/allocator.hpp"

#include<cstddef> // max_align_t
#include<stack>
#include<cassert>

namespace boost{
namespace multi{
namespace memory{

template<typename Ptr = byte*, std::size_t Align = alignof(std::max_align_t)>
class stack : protected monotonic<Ptr, Align>{
private:
	std::stack<typename stack::void_pointer> positions_ = {};
	typename stack::size_type total_requested_ = 0;
	typename stack::size_type total_discarded_ = 0;
	typename stack::size_type max_needed_ = 0;
	long hits_ = 0;
public:
	long hits() const{return hits_;}
	void reset(){monotonic<Ptr, Align>::reset(); positions_.clear();}
	typename stack::size_type max_needed() const{return max_needed_;}
	using monotonic<Ptr, Align>::monotonic;
	template<std::size_t AA = Align>
	typename stack::void_pointer allocate(
		typename stack::size_type required_bytes,
		typename stack::size_type align = AA
	){
		total_requested_ += required_bytes;
		max_needed_ = std::max(max_needed_, total_requested_ - total_discarded_);
		positions_.push(monotonic<Ptr, Align>::template allocate<AA>(required_bytes, align));
		++hits_;
		return positions_.top();
	}
	void deallocate(
		typename stack::void_pointer p, 
		typename stack::size_type discarded_bytes
	){
		total_discarded_ += discarded_bytes;
		monotonic<Ptr, Align>::deallocate(p, discarded_bytes);//positions_.top(), discarded_bytes);
		assert( p == positions_.top() && "stack violation" );
		this->position_ -= discarded_bytes;
		positions_.pop();
	}
};

template<class T> 
using stack_allocator = multi::memory::allocator<T, stack<char*>>;//, alignof(T)>>;

}}}

#if not __INCLUDE_LEVEL__
#define BOOST_TEST_MODULE "C++ Unit Tests for Multi stack memory resource"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../multi/array.hpp"
#include<vector>

namespace multi = boost::multi;
namespace memory = multi::memory;

BOOST_AUTO_TEST_CASE(multi_memory_allocator) {
	alignas(double) std::array<char, 256*sizeof(double)> buffer;
	memory::stack<char*> stck(buffer.data(), buffer.size());
	auto p1 = stck.allocate(1*sizeof(double), alignof(double));
	BOOST_REQUIRE( stck.max_needed() == 1*sizeof(double) );

	auto p2 = stck.allocate(255*sizeof(double), alignof(double));
	BOOST_REQUIRE( stck.max_needed() == 256*sizeof(double) );

	stck.deallocate(p2, 255*sizeof(double));
	BOOST_REQUIRE( stck.max_needed() == 256*sizeof(double) );

	stck.deallocate(p1, 1*sizeof(double));
	BOOST_REQUIRE( stck.max_needed() == 256*sizeof(double) );

	auto p3 = stck.allocate(100*sizeof(double)); (void)p3;
}
#endif
#endif
