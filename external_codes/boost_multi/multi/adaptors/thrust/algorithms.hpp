#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0.cu) && nvcc -I$HOME/prj/alf -D_TEST_MULTI_THRUST_ALGORITHMS $0.cu -o $0x && $0x && rm $0.cu $0x; exit
#endif

#ifndef BOOST_MULTI_ADAPTORS_THRUST_ALGORITHMS_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_ALGORITHMS_HPP

#include "../thrust/algorithms/copy.hpp"

//#include "../../array_ref.hpp"
//#include<thrust/device_ptr.h>
//#include <iostream>
#if 0
namespace boost{
namespace multi{
	template<typename From, typename To, typename = std::enable_if_t<std::is_trivially_assignable<To&, From>{}> >
	array_iterator<To, 1, To*> copy(
		array_iterator<From, 1, thrust::device_ptr<To>> f, 
		array_iterator<From, 1, thrust::device_ptr<To>> l, 
		array_iterator<To, 1, To*> d
	){
		assert(f.stride() == l.stride()); static_assert(sizeof(From) == sizeof(To), "!");
		auto n = std::distance(f, l);
		if(f.stride()==1 and d.stride()==1){
			auto s = cudaMemcpy(d.data(), raw_pointer_cast(f.data()), n*sizeof(To), cudaMemcpyDeviceToHost); assert( s == cudaSuccess );
		}else{
			auto s = cudaMemcpy2D(d.data(), d.stride()*sizeof(To), raw_pointer_cast(f.data()), f.stride()*sizeof(To), sizeof(To), n, cudaMemcpyDeviceToHost);
			assert( s == cudaSuccess );
		}
		return d + n;
	}

#if 1
	template<typename From, typename To, typename = std::enable_if_t<std::is_trivially_assignable<To&, From>{}> >
	void fill(
		array_iterator<To, 1, thrust::device_ptr<To>> f, 
		array_iterator<To, 1, thrust::device_ptr<To>> l, 
		From const& value
	){
		assert( f.stride() == l.stride() ); static_assert(sizeof(From) == sizeof(To), "!");
		auto n = std::distance(f, l);
		if(f.stride()==1){
		//	auto s = cudaMemcpy(f.data(), thrust::raw_pointer_cast(&value), n*sizeof(To), cudaMemcpyHostToDevice); assert(s == cudaSuccess);
		}else{
		//	auto s = cudaMemcpy2D(f.data(), f.stride()*sizeof(To), raw_pointer_cast(&value), 0*sizeof(To), sizeof(To), n, cudaMemcpyHostToDevice);
		//	assert( s == cudaSuccess );
		}
		return;
	}
#endif

}}
#endif

#ifdef _TEST_MULTI_THRUST_ALGORITHMS

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/device_allocator.h>

#include "../../array.hpp"

#include "../../adaptors/thrust/allocator_traits.hpp"

#include<numeric> // iota

int main(){
#if 0
	thrust::device_allocator<double> aaa;
	auto p = aaa.allocate(10);
	aaa.deallocate(p, 10);

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
  
	multi::array<int, 1, thrust::device_allocator<int>> H2(4, 99); assert(size(H2) == 4);
	assert( H2[2] == 99 );
	copy( begin(H), end(H), begin(H2) );
	assert( H2[1] == 20 );

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

