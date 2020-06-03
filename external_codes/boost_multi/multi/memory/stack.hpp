#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&& clang++ -std=c++17 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MULTI_MEMORY_STACK $0.cpp -o $0x && $0x && rm $0x $0.cpp; exit
#endif
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
		typename stack::size_type align = AA//alignof(std::max_align_t)
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

#if _TEST_BOOST_MULTI_MEMORY_STACK

#include "../../multi/array.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
namespace memory = multi::memory;

using std::cout;

int main(){
{
	alignas(double) char buffer[256*sizeof(double)];
	memory::stack<char*> stck(buffer);
	auto p1 = stck.allocate(1*sizeof(double), alignof(double));	  assert( stck.max_needed() == 1*sizeof(double) );
	auto p2 = stck.allocate(255*sizeof(double), alignof(double)); assert( stck.max_needed() == 256*sizeof(double) );
	stck.deallocate(p2, 255*sizeof(double));	                  assert( stck.max_needed() == 256*sizeof(double) );
	stck.deallocate(p1, 1*sizeof(double));
	assert( stck.max_needed() == 256*sizeof(double) );
	auto p3 = stck.allocate(100*sizeof(double)); (void)p3;
}
{
	alignas(double) char buffer[256*sizeof(double)];
	multi::memory::stack<char*> sm(buffer);
	{
		std::vector<double, multi::memory::stack_allocator<double>> v(10, &sm);
		std::vector<double, multi::memory::stack_allocator<double>> w(10, &sm);
	}
	std::vector<double, multi::memory::stack_allocator<double>> w(5, &sm);
	assert( sm.max_needed()/sizeof(double) == 20 );
}
}
#endif
#endif

