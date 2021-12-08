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

namespace boost{
namespace multi{
namespace memory{

template<class Ptr> // TODO test with actual fancy ptr
Ptr align_up(Ptr p, std::size_t align = alignof(std::max_align_t)){
	using multi::to_address;
	auto p_(to_address(p));
	static_assert( sizeof(*p_)==1 , "!"); // crash
	auto q_ = reinterpret_cast<decltype(p_)>(
		(reinterpret_cast<std::uintptr_t>(p_) + (align-1))
		& ~(align-1)
	);
	return p+std::distance(p_,q_);
}

template<class T>
T* align_up(T* p, std::size_t align = alignof(std::max_align_t)){
	return 
		reinterpret_cast<T*>(
			  (reinterpret_cast<std::uintptr_t>(p) + (align-1)) 
			& ~(align-1)
		)
	;
}

template<typename Ptr = byte*, std::size_t Align = alignof(std::max_align_t)>
class null_t{
public:
	using pointer = Ptr;
	using size_type = typename std::pointer_traits<pointer>::difference_type;
private:
	using void_pointer = typename std::pointer_traits<pointer>::template rebind<void>;
public:
	template<std::size_t AA = Align>
	void_pointer allocate(size_type required_bytes, size_type = AA){
		if(required_bytes > 0) throw std::bad_alloc{};		
		return nullptr;
	}
	void deallocate(void_pointer p, size_type /*discarded_bytes*/){
		if(p != nullptr) throw std::bad_alloc{};
	}
};

template<typename Ptr = byte*, std::size_t Align = alignof(std::max_align_t)>
class monotonic : protected block<Ptr>{
protected:
	using void_pointer = typename std::pointer_traits<typename monotonic::pointer>::template rebind<void>;
public:
	using block<Ptr>::start_;
	using block<Ptr>::size_;
	using block<Ptr>::block;
	typename monotonic::pointer position_ = start_;
	void reset(){position_ = this->start_;}
	template<std::size_t AA = Align>
	typename monotonic::void_pointer allocate(
		typename monotonic::size_type required_bytes,
		typename monotonic::size_type align = AA//alignof(std::max_align_t)
	){
		auto ret = align_up(this->position_, align);
		auto new_position_ = ret + required_bytes;
		using std::distance;
		if(not this->contains(new_position_-1))
			throw overflow(required_bytes, this->size_ - distance(start_, position_));
		this->position_ = new_position_;
		return ret;
	}
	bool owns(typename monotonic::void_pointer p) const{
		return this->contains(static_cast<typename monotonic::pointer>(p));
	}
	void deallocate(
		typename monotonic::void_pointer p, 
		typename monotonic::size_type /*discarded_bytes*/
	){
		if(not owns(p)) throw std::bad_alloc{};
	}
	struct overflow : public std::bad_alloc{
		using size_type = typename monotonic::size_type;
		size_type required;
		size_type available;
	//	constexpr auto to_string = [](auto a){return std::to_string(a);};
		std::string msg;
		overflow(size_type required, size_type available) 
		: required{required}, available{available}, 
			msg{"required "+std::to_string(required)+" while only "+std::to_string(available)+" bytes available"}
		{}
		virtual const char* what() const throw() override {return msg.c_str();}// + std::to_string(required)).c_str();}
	};
};

template<class T = void> 
using monotonic_allocator = multi::memory::allocator<T, monotonic<char*>>;

}}}

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
	try{
		mr.allocate(1*sizeof(double), alignof(double));
	}catch(...){}
}
{
	alignas(double) char buffer[256*sizeof(double)];
	multi::memory::monotonic<char*> m(buffer);
	auto p1 = m.allocate(1*sizeof(double), alignof(double));
	auto p2 = m.allocate(255*sizeof(double), alignof(double));
	m.deallocate(p2, 255*sizeof(double));
	m.deallocate(p1, 1*sizeof(double));
	try{
		m.deallocate((char*)p1 + 10000, 1*sizeof(double));
	}catch(...){}
}
{
	alignas(double) char buffer[300*sizeof(double)];
	multi::memory::monotonic<char*> m(&buffer[0], 300*sizeof(double));
	multi::memory::monotonic_allocator<double> alloc(&m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> A({10, 10}, &m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> B({10, 10}, &m);
	multi::array<double, 2, multi::memory::monotonic_allocator<double>> C({10, 10}, &m);
}
//	m.allocate(156*sizeof(double), alignof(double));



//	m.allocate(256*sizeof(double));
//	typename monotonic::pointer

#if 0
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
#endif
}
#endif
#endif

