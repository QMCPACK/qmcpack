// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR_HPP
#define BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR_HPP

#include <algorithm>   // for max
#include <cassert>
#include <cstddef>     // for max_align_t
#include <memory>
#include <stack>

#include "generic_allocator.hpp"

namespace boost {
namespace multi {

template<typename Ptr = char*>
struct block : std::pointer_traits<Ptr>{
	template<std::size_t N> block(char(&t)[N]) : start_{t}, lenght_{N} {}
	typename block::pointer start_;
	typename block::size_type lenght_;
	bool contains(typename block::pointer p) const{
		using std::distance;
		return distance(start_, p) < static_cast<difference_type>(lenght_);
	}
};

template<class Ptr = char*, std::size_t A = alignof(std::max_align_t)>
class monotonic_buffer : block<Ptr> {
	size_type hits_ = 0;
	size_type misses_ = 0;
	size_type allocated_bytes_ = 0;
	size_type deallocated_bytes_ = 0;
//  size_type position_ = 0;
	static std::size_t max_alignment = A;

	static size_type align_up(size_type n) noexcept {
		return (n + (max_alignment-1)) & ~(max_alignment-1);
	}

 public:
	size_type size() const {return size_;}
	size_type hits() const {return hits_;}
	size_type misses() const {return misses_;}

	size_type allocated_bytes() const {return allocated_bytes_;}
	size_type deallocated_bytes() const {return deallocated_bytes_;}

	using block<Ptr>::block;
	monotonic_buffer& operator=(block<Ptr> const& b) {
		position_ = 0;
		block<Ptr>::operator=(b); return *this;
	}

	monotonic_buffer(monotonic_buffer const&) = delete;
	monotonic_buffer& operator=(monotonic_buffer const&) = delete;

	~monotonic_buffer() {
#ifndef _BOOST_MULTI_RELAX_MONOTONIC_LEAK_CONDITION
		assert(allocated_bytes() == deallocated_bytes());
#endif
	}
	using void_pointer = monotonic_buffer;

	template<size_type RequiredAlignment = sizeof(std::max_align_t)>
	void_pointer allocate(size_type req_bytes, size_type al = RequiredAlignment) {
		static_assert(RequiredAlignment <= max_alignment, "!");  // requested alignment is too large for this MR

		assert( al <= max_alignment );
		auto bytes = this->align_up(req_bytes);
		if(position_ + bytes < size_) {
			auto old_position_ = position_;
			position_ += bytes;
			++hits_;
			allocated_bytes_ += bytes;
			return buffer_ + old_position_;
		}
		++misses_;
		auto p = alloc_.allocate(bytes/sizeof(std::max_align_t));
		if(p) allocated_bytes_ += bytes;
		return p;
	}
	void deallocate(void_pointer p, size_type req_bytes){
		auto bytes = align_up(req_bytes);
		deallocated_bytes_ += bytes;
		if(not in_buffer(static_cast<char_pointer>(p))) {
			alloc_.deallocate(static_cast<pointer>(p), bytes);
		}
	}
};

template<class T = void>
using monotonic_allocator = multi::generic_allocator<T, multi::monotonic_buffer<std::allocator<char>>>;

}  // end namespace multi
}  // end namespace boost

#if _TEST_BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR

#include "../../multi/array.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
using std::cout;

int main() {
	std::cout << sizeof(std::max_align_t) <<" "<< alignof(std::max_align_t) << std::endl;

	 {
		multi::monotonic_buffer<> buf(250*sizeof(double));
		 {
			multi::array<double, 2, multi::monotonic_allocator<> > A({10, 10}, &buf);
			multi::array<double, 2, multi::monotonic_allocator<> > B({10, 10}, &buf);
			multi::array<double, 2, multi::monotonic_allocator<> > C({10, 10}, &buf);
		}
		assert( buf.hits() == 2 );
		assert( buf.misses() == 1 );
		cout
			<<"size: "<< buf.size()
			<<"\nsaved: "<< buf.hits()
			<<"\nmisses "<< buf.misses()
			<<"\nallocated(bytes) "<< buf.allocated_bytes()
			<<"\nreleased(bytes) "<< buf.deallocated_bytes()
			<< std::endl;
	}
	cout<<"----------"<<std::endl;
	 {
		std::size_t guess = 0;
		multi::monotonic_buffer<std::allocator<char>> buf(guess*sizeof(double));

		for(int i = 0; i != 3; ++i) {
			cout<<"pass "<< i << std::endl;
			 {
				multi::array<double, 2, multi::monotonic_allocator<> > A({10, 10}, &buf);
				multi::array<double, 2, multi::monotonic_allocator<> > B({10, 10}, &buf);
				multi::array<double, 2, multi::monotonic_allocator<> > C({10, 10}, &buf);
			}
			cout
				<<"  size: "<< buf.size()
				<<"\n  save: "<< buf.hits()
				<<"\n  misses "<< buf.misses()
				<<"\n  allocated(bytes) "<< buf.allocated_bytes()
				<<"\n  released(bytes) "<< buf.deallocated_bytes()
				<< std::endl;
		//  guess = std::max(guess, buf.allocated_bytes());
			buf.reset(buf.allocated_bytes());
		}
	}
	cout<<"----------monotonic"<<std::endl;
	 {
		std::size_t guess_bytes = 120;

		for(int i = 0; i != 3; ++i) {
			cout<<"pass "<< i << std::endl;
			multi::monotonic_buffer<std::allocator<char>> buf(guess_bytes);

			 {
				multi::array<double, 2, multi::monotonic_allocator<> > A({10, 10}, &buf);

				for(int i = 0; i != 3; ++i) {
					multi::array<double, 2, multi::monotonic_allocator<> > B({10, 10}, &buf);
					std::vector<int, multi::monotonic_allocator<int>> v(3, &buf);
					v.push_back(33); v.push_back(33);
				}

				multi::array<double, 2, multi::monotonic_allocator<> > C({10, 10}, &buf);
			}

			cout
				<<"  size: "<< buf.size()
				<<"\n  hits: "<< buf.hits()
				<<"\n  misses "<< buf.misses()
				<<"\n  allocated(bytes) "<< buf.allocated_bytes()
				<<"\n  released(bytes) "<< buf.deallocated_bytes()
				<< std::endl;
			guess_bytes = std::max(guess_bytes, buf.allocated_bytes());
		}
	}
}
#endif
#endif

