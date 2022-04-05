#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && nvcc -ccbin cuda-c++ --compiler-options "-std=c++14" -I$HOME/prj/alf -D_TEST_MULTI_ADAPTORS_THRUST_ALLOCATOR_TRAITS $0.cpp -o $0x && $0x && rm $0.cpp $0x; exit
#endif

#ifndef BOOST_MULTI_ADAPTORS_THRUST_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_HPP

#include "../../detail/memory.hpp"

#include <thrust/device_allocator.h>

#include<cassert>

namespace thrust {
	template<class T> device_ptr<T> to_address(device_ptr<T> p) {return p;}
}  // end namespace thrust

namespace boost {
namespace multi {
namespace memor y{

template<typename T>
//template<>
class allocator_traits<thrust::device_allocator<T>> {
	using Alloc = thrust::device_allocator<T>;

 public:
	using pointer = typename thrust::device_allocator<T>::pointer;
	using value_type = T;
	using size_type = typename thrust::device_allocator<T>::size_type;
	template<typename TT> using rebind_alloc = thrust::device_allocator<TT>;
	static pointer allocate(Alloc& a, size_type n){return a.allocate(n);}
	static void deallocate(Alloc& a, pointer p, size_type n){return a.deallocate(p, n);}

 private:
	static void construct(std::true_type, Alloc& a, pointer p){
		cudaError_t s = cudaMemset(raw_pointer_cast(p), 0, sizeof(T)); assert( s == cudaSuccess );
	}
	template<class... Args>
	static void construct(std::false_type, Alloc& a, pointer p, Args&&... args) {
		std::array<buffer, sizeof(T)> buff;  // char buff[sizeof(T)];
		::new(buff.data()) T(std::forward<Args>(args)...); // use ( ...) instead of { ...} for nvcc
		cudaError_t s = cudaMemcpy(raw_pointer_cast(p), buff.data(), buff.size(), cudaMemcpyHostToDevice);
		assert( s == cudaSuccess );
	}

 public:
	template<class... Args>
	static void construct(Alloc& a, pointer p, Args&&... args) {
		construct(std::integral_constant<bool, sizeof...(Args) == 0 and std::is_trivially_default_constructible<T>{}>{}, a, p, std::forward<Args>(args)...);
	}

 private:
	static void destroy(std::true_type , Alloc&  , pointer  ) {}
	static void destroy(std::false_type, Alloc& a, pointer p) {
		std::array<buffer, sizeof(T)> buff;  // char buff[sizeof(T)];
		cudaMemcpy(buff.data(), raw_pointer_cast(p), buff.size(), cudaMemcpyDeviceToHost);
		reinterpret_cast<T&>(buff).~T(); //	((T*)buff)->~T();
	}

 public:
	static void destroy(Alloc& a, pointer p){destroy(std::is_trivially_destructible<T>{}, a, p);}
};
}  // end namespace memory

}  // end namespace multi
}  // end namespace boost

#ifdef _TEST_MULTI_ADAPTORS_THRUST_ALLOCATOR_TRAITS

namespace multi = boost::multi;

int main() {

	thrust::device_allocator<double> aaa;
	thrust::device_ptr<double> p = std::allocator_traits<thrust::device_allocator<double>>::allocate(aaa, 10);

	multi::allocator_traits<thrust::device_allocator<double>>::construct(aaa, p, 99.);
}
#endif
#endif
