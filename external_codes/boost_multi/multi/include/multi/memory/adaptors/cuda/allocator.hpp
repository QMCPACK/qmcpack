#ifdef COMPILATION// -*- indent-tabs-mode:t;c-basic-offset:4;tab-width:4; -*-
$CXXX $CXXFLAGS $0 -o $0x -lcudart -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#ifndef MULTI_MEMORY_ADAPTORS_CUDA_ALLOCATOR_HPP
#define MULTI_MEMORY_ADAPTORS_CUDA_ALLOCATOR_HPP

#include<cuda_runtime.h> // cudaMalloc

#include "../../adaptors/cuda/ptr.hpp"
#include "../../adaptors/cuda/algorithm.hpp"

#include "../../adaptors/cuda/clib.hpp"    // cuda::malloc
#include "../../adaptors/cuda/cstring.hpp" // cuda::memcpy
#include "../../adaptors/cuda/malloc.hpp"

#include<new>      // bad_alloc
#include<cassert>
#include<iostream> // debug

#include<complex>

#include<cstddef>

namespace boost{namespace multi{
namespace memory{namespace cuda{

struct bad_alloc : std::bad_alloc {};

//struct allocation_counter {
//	static long n_allocations;
//	static long n_deallocations;
//	static long bytes_allocated;
//	static long bytes_deallocated;
//};

//long allocation_counter::n_allocations = 0;
//long allocation_counter::n_deallocations = 0;
//long allocation_counter::bytes_allocated = 0;
//long allocation_counter::bytes_deallocated = 0;

template<class T = void>
class allocator {//: protected allocation_counter {
	static_assert(std::is_same<T, std::decay_t<T>>{},
		"allocated type should be a value type, not a reference or decorated type");

 public:
	using value_type = T;
	using pointer = ptr<T>;
	using const_pointer = ptr<T const>;
	using void_pointer = ptr<void>;
	using const_void_pointer = ptr<void const>;
	using difference_type = typename pointer::difference_type;
	template<class TT> using rebind = allocator<TT>;
	using size_type = ::size_t; // as specified by CudaMalloc

	allocator() = default;
	template<class U>
	allocator(allocator<U> const& /*other*/) noexcept {}

	pointer allocate(size_type n, const_void_pointer = 0) {//const void* = 0) {
		if(n == 0) return pointer{nullptr};
		auto ret = static_cast<pointer>(cuda::malloc(n*sizeof(T)));
		if(not ret) throw bad_alloc{};
	//  ++n_allocations; bytes_allocated+=sizeof(T)*n;
		return ret;
	}
	void deallocate(pointer p, size_type n) {
		cuda::free(p);
	//  ++n_deallocations; bytes_deallocated+=sizeof(T)*n;
	}

	std::true_type  operator==(allocator const&) const {return {};}
	std::false_type operator!=(allocator const&) const {return {};}

	template<class P, class... Args>
	[[deprecated("cuda slow")]]
	void construct(P p, Args&&... args) = delete;/*{
		if(sizeof...(Args) == 0 and std::is_trivially_default_constructible<T>{})
			cuda::memset(p, 0, sizeof(T));
		else{
			char buff[sizeof(T)];
			::new(buff) T(std::forward<Args>(args)...);
			cuda::memcpy(p, buff, sizeof(T));
		}
	}*/
	template<class P> 
	[[deprecated("cuda slow")]]
	void destroy(P p) {
		if(not std::is_trivially_destructible<T>{}) {
			std::array<char, sizeof(T)> buff;  // char buff[sizeof(T)];
			cuda::memcpy(buff.data(), p, buff.size());
			((T*)buff)->~T();
		}
	}

#if 0
	template<class InputIt, class Size, class ForwardIt>//, typename T1 = typename std::iterator_traits<ForwardIt>::value_type>
	auto alloc_uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first)
	DECLRETURN(adl_uninitialized_copy_n(first, count, d_first))

	template<class InputIt, class Size, class ForwardIt, typename T1 = typename std::iterator_traits<ForwardIt>::value_type>
	auto alloc_uninitialized_move_n(InputIt first, Size count, ForwardIt d_first)
	DECLRETURN(adl_uninitialized_move_n(first, count, d_first))

	template<class InputIt, class ForwardIt, typename T1 = typename std::iterator_traits<ForwardIt>::value_type>
	auto alloc_uninitialized_copy(InputIt first, InputIt last, ForwardIt d_first)
	DECLRETURN(adl_uninitialized_copy(first, last, d_first))

//	DECLRETURN(adl_uninitialized_copy(first, count, d_first))
	template<class InputIt, class Size, class ForwardIt, typename T1 = typename std::iterator_traits<ForwardIt>::value_type>
	auto alloc_uninitialized_copy(InputIt first, Size count, ForwardIt d_first) 
	DECLRETURN(uninitialized_copy(first, count, d_first))
	template<class Ptr, class Size, class V = typename Ptr::element_type>//, std::enable_if_t<typename Ptr::element_type> >
	auto alloc_uninitialized_value_construct_n(Ptr p, Size n)
	DECLRETURN(uninitialized_value_construct_n(p, n))

	template<class Ptr, typename Size>
	auto alloc_uninitialized_default_construct_n(Ptr p, Size n)
	DECLRETURN(uninitialized_default_construct_n(p, n))

	template<class Ptr, class Size, class V, std::enable_if_t<std::is_trivially_copy_constructible<V>{}, int> =0>// = typename Ptr::element_type>
	Ptr alloc_uninitialized_fill_n(Ptr p, Size n, V const& v){
		return uninitialized_fill_n(p, n, v);}
	template<class TT> 
	static std::true_type  is_complex_(std::complex<TT>);
	static std::false_type is_complex_(...);
	template<class TT> struct is_complex : decltype(is_complex_(TT{})){};
	template<
		class Ptr, class Size, class V = typename Ptr::element_type, 
		std::enable_if_t<std::is_trivially_default_constructible<V>{} or is_complex<V>{}, int> = 0
	>
	Ptr alloc_uninitialized_default_construct_n(Ptr const& p, Size n) const{return p + n;}
	template<class Ptr, class Size, class V = typename Ptr::element_type>
	Ptr alloc_destroy_n(Ptr p, Size n){
		if(std::is_trivially_destructible<V>{}) {
		} else {assert(0);}
		return p + n;
	}
#endif
};

template<>
class allocator<std::max_align_t> {//: allocation_counter{
 public:
	using T = std::max_align_t;
	using value_type = T;
	using pointer = ptr<T>;

//	using void_pointer = ptr<void>;
	using const_void_pointer = ptr<void const>;
//	using difference_type = typename pointer::difference_type;

	using size_type = ::size_t; // as specified by CudaMalloc
	auto allocate(size_type n, const_void_pointer = 0) {  // const void* = 0){
		if(n == 0) return pointer{nullptr};
		auto ret = static_cast<pointer>(cuda::malloc(n*sizeof(T)));
		if(not ret) throw bad_alloc{};
	//  ++n_allocations; bytes_allocated+=sizeof(T)*n;
		return ret;
	}
	void deallocate(pointer p, size_type n) {
		cuda::free(p);
	//  ++n_deallocations; bytes_deallocated+=sizeof(T)*n;
	}
	std::true_type  operator==(allocator<std::max_align_t> const&) const {return {};} // template explicit for nvcc
	std::false_type operator!=(allocator<std::max_align_t> const&) const {return {};}

	template<class P, class... Args>
	void construct(/*[[maybe_unused]]*/ P p, Args&&...) {(void)p; assert(0);} // TODO investigate who is calling this
	template<class P>
	void destroy(P) {} // TODO(correaa) investigate who is calling this
};

}}}}

namespace std {

#if __NVCC__ // this solves this error with nvcc error: ‘template<class _Tp> using __pointer = typename _Tp::pointer’ is protected within this context
template<class T>
class allocator_traits<boost::multi::memory::cuda::allocator<T>> {
	using Alloc = boost::multi::memory::cuda::allocator<T>;

 public:
	using allocator_type = Alloc;
	using value_type = typename Alloc::value_type;
	using pointer = typename Alloc::pointer;
	using const_pointer = typename Alloc::const_pointer;
	using void_pointer = typename Alloc::void_pointer;
	using const_void_pointer = typename Alloc::const_void_pointer;
	using difference_type = typename Alloc::difference_type;
	using size_type = typename Alloc::size_type;
	using propagate_on_container_copy_assignment = std::false_type;
	using propagate_on_container_move_assignment = std::false_type;
	using propagate_on_container_swap 	         = std::false_type;
	template<class T2>
	using rebind_alloc = typename Alloc::template rebind<T2>;

	static constexpr Alloc select_on_container_copy_construction(Alloc const& a) {return a;}

	template<class...As> static auto deallocate(allocator_type& a, As&&... as) {return a.deallocate(std::forward<As>(as)...);}
	template<class...As> static auto   allocate(allocator_type& a, As&&... as) {return a.  allocate(std::forward<As>(as)...);}
};
#endif

}  // end namespace std

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__
//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi memory allocator"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include<memory>
//#include<iostream>
//#include<complex>

//#include "../../../array.hpp"
//#include "../cuda/algorithm.hpp"

//namespace multi = boost::multi;
//namespace cuda  = multi::memory::cuda;

//void add_one(double& d){d += 1.;}
//template<class T> void add_one(T&& t){std::forward<T>(t) += 1.;}

//template<class T> void what(T&&) = delete;
//using std::cout;

//BOOST_AUTO_TEST_CASE(multi_memory_allocator){
//	{
//		multi::static_array<double, 1> A(32, double{}); A[17] = 3.;
//		multi::static_array<double, 1, cuda::allocator<double>> A_gpu = A;
//		BOOST_REQUIRE( A_gpu[17] == 3 );
//	}
//}
//#endif
#endif

