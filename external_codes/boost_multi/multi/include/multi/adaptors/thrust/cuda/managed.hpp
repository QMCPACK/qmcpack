#pragma once


#include "../../cuda/runtime/error.hpp"

#include <thrust/system/cuda/pointer.h>

#include<new> // bad_alloc
#include<cassert>

namespace boost{
namespace multi{

namespace thrust{
namespace cuda{
namespace managed{

template<class> class pointer;

template<class T>
class reference : public ::thrust::cuda::reference<T>{
	using base_type = ::thrust::cuda::reference<T>;
public:
	constexpr explicit reference(::thrust::cuda::reference<T> const& other) : base_type{other}{}
	constexpr explicit reference(T& other) : base_type{&other}{}
	constexpr operator T&()&&{return raw_reference_cast(static_cast<base_type&>(*this));}
	constexpr pointer<T> operator&(){return pointer<T>{base_type::operator&()};}
	using ::thrust::cuda::reference<T>::operator=;
};

template<class T>
class pointer {
	::thrust::cuda::pointer<T> impl_;

 public:
	constexpr explicit pointer(::thrust::cuda::pointer<T> const& other) : impl_{other} {}
	constexpr explicit pointer(T* other) : impl_(other) {}

	using difference_type   = typename ::thrust::iterator_traits<::thrust::cuda::pointer<T>>::difference_type;
	using value_type        = typename ::thrust::iterator_traits<::thrust::cuda::pointer<T>>::value_type;
#pragma push
#pragma nv_diag_suppress = class_and_member_name_conflict  // for nvcc warning: declaration of a member with the same name as its class TODO(correaa) switch to new stype pragma for nvcc and nvc++
	using pointer           = pointer<T>;
#pragma pop
	using reference         = managed::reference<T>;
	using iterator_category = typename ::thrust::iterator_traits<::thrust::cuda::pointer<T>>::iterator_category;

	using element_type = T;

	constexpr operator T*() const{return raw_pointer_cast(impl_);}
	constexpr operator ::thrust::cuda::pointer<T>() const{return impl_;}

	constexpr pointer& operator++(){impl_.operator++(); return *this;}
	constexpr pointer& operator--(){impl_.operator--(); return *this;}
	constexpr auto operator++(int i){return pointer{impl_.operator++(i)};}
	constexpr auto operator--(int i){return pointer{impl_.operator--(i)};}
	constexpr pointer& operator+=(difference_type n){impl_.operator+=(n); return *this;}
	constexpr pointer& operator-=(difference_type n){impl_.operator-=(n); return *this;}
	constexpr pointer operator+(difference_type n) const{return pointer{impl_ + n};}
	constexpr pointer operator-(difference_type n) const{return pointer{impl_ - n};}

	constexpr reference operator*() const{return reference{impl_.operator*()};}
	constexpr reference operator[](difference_type n){return *((*this)+n);}

	friend auto raw_pointer_cast(pointer const& p){return raw_pointer_cast(p.impl_);}
};

struct bad_alloc : std::bad_alloc{};

template<class T = void>
class allocator{// : cuda::allocator<T>{
	static_assert( std::is_same<T, std::decay_t<T>>{}, "!" );
public:
	using value_type = T;
	using pointer = managed::pointer<T>;
	using size_type = ::size_t; // as specified by CudaMalloc
	pointer allocate(typename allocator::size_type n, const void* = 0){
		if(n == 0) return pointer{nullptr};
		T* p = nullptr;
		namespace cudart = boost::multi::cuda::runtime;
		auto e = static_cast<cudart::error>(cudaMallocManaged(&p, n*sizeof(T)));
		switch(e){
			case cudart::success : break;
			case cudart::memory_allocation: throw bad_alloc{};
			default: throw std::system_error{e, "cannot allocate "+std::to_string(n*sizeof(T))+" bytes in '"+__PRETTY_FUNCTION__+"'"};
		}
		auto ret = static_cast<pointer>(p);
		if(!ret) throw bad_alloc{};
		return ret;
	}
	void deallocate(pointer p, size_type){
		namespace cudart = boost::multi::cuda::runtime;
		auto e = static_cast<cudart::error>(cudaFree(raw_pointer_cast(p)));
		if(e!=cudart::success){
			throw std::system_error{e, std::string{"cannot "}+ __PRETTY_FUNCTION__};
		}
	}
	template<class P, class... Args>
	void construct(P p, Args&&... args){ // remove?
		::new(p.rp_) T(std::forward<Args>(args)...);
	}
	template<class P, class... Args>
	void construct(P* p, Args&&... args){ // remove?
		::new(p) T(std::forward<Args>(args)...);
	}
	template<class P> void destroy(P p){p.rp_->~T();} // remove?
	template<class P> void destroy(P* p){p->~T();}    // remove?
	constexpr bool operator==(allocator<T> const&) const{return true;}
	constexpr bool operator!=(allocator<T> const&) const{return false;}

	template<class InputIt, class ForwardIt>
	constexpr ForwardIt alloc_uninitialized_copy(InputIt first, InputIt last, ForwardIt d_first) const{
		return ForwardIt{adl_uninitialized_copy(first, last, d_first)};
	}
	template<class InputIt, class Size, class ForwardIt>
	constexpr ForwardIt alloc_uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first) const{
		return ForwardIt{adl_uninitialized_copy_n(first, count, d_first)};
	}
	template<class ForwardIt, class Size>
	constexpr ForwardIt alloc_uninitialized_default_construct_n(ForwardIt first, Size n) const{
		return ForwardIt{adl_uninitialized_default_construct_n(first, n)};
	}
	template<class ForwardIt, class Size>
	constexpr ForwardIt alloc_destroy_n(ForwardIt first, Size n) const{return ForwardIt{destroy_n(first, n)};}
};

}}}

}}

#if not __INCLUDE_LEVEL__

#include<memory>
#include<iostream>
#include "../../../../array.hpp"

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

int main(){

	multi::array<double, 1, multi::memory::cuda::managed::allocator<double> > A(32);
	A[17] = 3.;
	assert( A[17] == 3. );

}
#endif
