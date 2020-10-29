#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cpp) && nvcc -ccbin cuda-c++ --compiler-options "-std=c++14" -I$HOME/prj/alf -D_TEST_MULTI_ADAPTORS_THRUST_ALLOCATOR_TRAITS $0.cpp -o $0x && $0x && rm $0.cpp $0x; exit
#endif

#ifndef BOOST_MULTI_ADAPTORS_THRUST_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_HPP

//#include <thrust/host_vector.h>
//#include <thrust/device_vector.h>
#include "../../detail/memory.hpp"

#include <thrust/device_allocator.h>

#include<cassert>


namespace thrust{
	template<class T> device_ptr<T> to_address(device_ptr<T> p){return p;}
}

namespace boost{
namespace multi{
namespace memory{

template<typename T>
//template<>
class allocator_traits<thrust::device_allocator<T>>{
//	using T = double;
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
	static void construct(std::false_type, Alloc& a, pointer p, Args&&... args){
		char buff[sizeof(T)]; 
		::new(buff) T(std::forward<Args>(args)...); // use ( ...) instead of { ...} for nvcc
		cudaError_t s = cudaMemcpy(raw_pointer_cast(p), buff, sizeof(T), cudaMemcpyHostToDevice); 
		assert( s == cudaSuccess );
	}
public:
	template<class... Args>
	static void construct(Alloc& a, pointer p, Args&&... args){
		construct(std::integral_constant<bool, sizeof...(Args) == 0 and std::is_trivially_default_constructible<T>{}>{}, a, p, std::forward<Args>(args)...);
	}
private:
	static void destroy(std::true_type , Alloc&  , pointer  ){}
	static void destroy(std::false_type, Alloc& a, pointer p){
		char buff[sizeof(T)];
		cudaMemcpy(buff, raw_pointer_cast(p), sizeof(T), cudaMemcpyDeviceToHost);
		reinterpret_cast<T&>(buff).~T(); //	((T*)buff)->~T();
	}
public:
	static void destroy(Alloc& a, pointer p){destroy(std::is_trivially_destructible<T>{}, a, p);}
};
}

}}

#ifdef _TEST_MULTI_ADAPTORS_THRUST_ALLOCATOR_TRAITS

namespace multi = boost::multi;

int main(){

	thrust::device_allocator<double> aaa;
	thrust::device_ptr<double> p = std::allocator_traits<thrust::device_allocator<double>>::allocate(aaa, 10);
//	aaa.deallocate(p, 10);

	multi::allocator_traits<thrust::device_allocator<double>>::construct(aaa, p, 99.);
//	assert( p[0] == 99. );

#if 0
	namespace multi = boost::multi;
	multi::array<double, 1, thrust::device_allocator<double>> A({100}, 0.);
	A[20] = 44.;

	multi::array<double, 1> A_host({100}, 99.);
	{
	//	multi::array<double, 1, thrust::device_allocator<double>> Adev({10}, 0.); std::iota(begin(Adev), end(Adev), 0.);
		multi::array<double, 2, thrust::device_allocator<double>> Adev({3, 3}, 0.); std::iota(begin(Adev[2]), end(Adev[2]), 5.);

		std::cout <<"iota? "<< Adev[0][0] <<" "<< Adev[0][1] <<" "<< Adev[0][2] <<" " << std::endl;
		std::cout <<"iota? "<< Adev[1][0] <<" "<< Adev[1][1] <<" "<< Adev[1][2] <<" " << std::endl;
		std::cout <<"iota? "<< Adev[2][0] <<" "<< Adev[2][1] <<" "<< Adev[2][2] <<" " << std::endl;
		std::cout <<"----"<< std::endl;

		multi::array<double, 2> Ahos({3, 3}, 0.);
	//	assert( Ahos[1].size() == Adev[2].size());
		Ahos.rotated()[1] = Adev[2];
	//	Ahos.rotated()[0] = Adev[2];
		std::cout <<"iota? "<< Ahos[0][0] <<" "<< Ahos[0][1] <<" "<< Ahos[0][2] <<" "<< std::endl;
		std::cout <<"iota? "<< Ahos[1][0] <<" "<< Ahos[1][1] <<" "<< Ahos[1][2] <<" "<< std::endl;
		std::cout <<"iota? "<< Ahos[2][0] <<" "<< Ahos[2][1] <<" "<< Ahos[2][2] <<" "<< std::endl;
	}
	return 0;
	{
		namespace multi = boost::multi;
		multi::array<double, 2, thrust::device_allocator<double>> A({100, 100}, 5.);
		multi::array<double, 2> B({100, 100}, 3.);
	//	B = A;
	//	assert( B[3][2] == 5. );
	//	multi::array<double, 2, thrust::device_allocator<double>> A({100, 100}, 0.);
	}

	assert( A_host[25] == 99. );
	copy(begin(A), end(A), begin(A_host));
	A_host = A; 
	assert( size(A_host) == 100 );
	std::cout << A_host[20] << std::endl;
	assert(A_host[20] == 44. );

    // H has storage for 4 integers
    thrust::host_vector<int> H(4);

    // initialize individual elements
    H[0] = 14;
    H[1] = 20;
    H[2] = 38;
    H[3] = 46;
    
    // H.size() returns the size of vector H
    std::cout << "H has size " << H.size() << std::endl;

    // print contents of H
    for(int i = 0; i < H.size(); i++)
        std::cout << "H[" << i << "] = " << H[i] << std::endl;

    // resize H
    H.resize(2);
    
    std::cout << "H now has size " << H.size() << std::endl;

    // Copy host_vector H to device_vector D
    thrust::device_vector<int> D = H;
    
    // elements of D can be modified
    D[0] = 99;
    D[1] = 88;
    
    // print contents of D
    for(int i = 0; i < D.size(); i++)
        std::cout << "D[" << i << "] = " << D[i] << std::endl;

    // H and D are automatically deleted when the function returns
    return 0;
#endif
}
#endif
#endif

