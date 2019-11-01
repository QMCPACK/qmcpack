#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && c++ -std=c++14 -Wall -Wextra -Wfatal-errors -D_TEST_BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR $0x.cpp -o $0x.x && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
#ifndef BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR_HPP
#define BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR_HPP

#include<cstddef> // max_align_t
#include<memory>
#include<stack>
#include<cassert>
#include<iostream>

#include "generic_allocator.hpp"

namespace boost{
namespace multi{

template<class Alloc = std::allocator<char>, typename std::allocator_traits<Alloc>::size_type MaxAlignemnt = sizeof(std::max_align_t)>
struct monotonic_buffer : std::allocator_traits<Alloc>::template rebind_alloc<std::max_align_t>{
	using allocator_type = typename std::allocator_traits<Alloc>::template rebind_alloc<std::max_align_t>;
private:
	using allocator_traits = std::allocator_traits<allocator_type>;
public:
	using size_type       = typename allocator_traits::size_type;
	using pointer         = typename allocator_traits::pointer;
	using void_pointer    = typename std::pointer_traits<pointer>::template rebind<void>;
	using char_pointer    = typename std::pointer_traits<pointer>::template rebind<char>;
	using difference_type = typename allocator_traits::difference_type;
	constexpr static const size_type max_alignment = MaxAlignemnt;
protected:
	allocator_type alloc_;
	char_pointer buffer_;
	size_type size_;
	size_type position_ = 0;
	static size_type align_up(size_type n) noexcept{return (n + (max_alignment-1)) & ~(max_alignment-1);}
	long hits_ = 0;
	long misses_ = 0;
	size_type allocated_bytes_ = 0;
	size_type deallocated_bytes_ = 0;
	bool in_buffer(char_pointer p) const{
		return std::distance(buffer_, static_cast<char_pointer>(p)) <= static_cast<difference_type>(size_);
	}
public:
	allocator_type& allocator(){return alloc_;}
	size_type size() const{return size_;}
	long hits() const{return hits_;}
	long misses() const{return misses_;}
	size_type allocated_bytes() const{return allocated_bytes_;}
	size_type deallocated_bytes() const{return deallocated_bytes_;}
	monotonic_buffer(size_type bytes = 0, allocator_type alloc = {}) : 
		allocator_type{alloc}, buffer_{static_cast<char_pointer>(static_cast<void_pointer>(allocator_type::allocate(align_up(bytes)/sizeof(std::max_align_t))))}, size_{align_up(bytes)}{}
	void reset(size_type bytes = 0){
		if(bytes > size_){
			allocator_type::deallocate(static_cast<pointer>(static_cast<void_pointer>(buffer_)), size_/sizeof(std::max_align_t));
			buffer_ = static_cast<char_pointer>(static_cast<void_pointer>(allocator_type::allocate(align_up(bytes)/sizeof(std::max_align_t))));
			size_ = align_up(bytes);
		}
	}
	monotonic_buffer(monotonic_buffer const&) = delete;
	monotonic_buffer& operator=(monotonic_buffer const&) = delete;
	~monotonic_buffer(){
		#ifndef _BOOST_MULTI_RELAX_MONOTONIC_LEAK_CONDITION
		assert(allocated_bytes() == deallocated_bytes());
		#endif
		allocator_type::deallocate(static_cast<pointer>(static_cast<void_pointer>(buffer_)), size_/sizeof(std::max_align_t));
	}
	template<size_type RequiredAlignment = sizeof(std::max_align_t)>
	void_pointer allocate(size_type req_bytes, size_type al = RequiredAlignment){
		static_assert( RequiredAlignment <= max_alignment, "!"); // requested alignment is too large for this MR
		assert( al <= max_alignment );
		auto bytes = this->align_up(req_bytes);
		if(position_ + bytes < size_){
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
		if(not in_buffer(static_cast<char_pointer>(p))) 
			alloc_.deallocate(static_cast<pointer>(p), bytes);
	}
};

template<class T = void>
using monotonic_allocator = multi::generic_allocator<T, multi::monotonic_buffer<std::allocator<char>>>;

}}

#if _TEST_BOOST_MULTI_DETAIL_MONOTONIC_ALLOCATOR

#include "../../multi/array.hpp"

#include<iostream>
#include<vector>
#include<cmath>

namespace multi = boost::multi;
using std::cout;

int main(){
	std::cout << sizeof(std::max_align_t) << " " << alignof(std::max_align_t) << std::endl;
//	return 0;
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
			<< std::endl
		;
	}
	cout<<"----------"<<std::endl;
	{
		std::size_t guess = 0;
		multi::monotonic_buffer<std::allocator<char>> buf(guess*sizeof(double));
		for(int i = 0; i != 3; ++i){
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
				<< std::endl
			;
		//	guess = std::max(guess, buf.allocated_bytes());
			buf.reset(buf.allocated_bytes());
		}
	}
	cout<<"----------monotonic"<<std::endl;
	{
		std::size_t guess_bytes = 120;
		for(int i = 0; i != 3; ++i){
			cout<<"pass "<< i << std::endl;
			multi::monotonic_buffer<std::allocator<char>> buf(guess_bytes);
			{
				multi::array<double, 2, multi::monotonic_allocator<> > A({10, 10}, &buf);
				for(int i = 0; i != 3; ++i){
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
				<< std::endl
			;
			guess_bytes = std::max(guess_bytes, buf.allocated_bytes());
		}
	}
}
#endif
#endif

