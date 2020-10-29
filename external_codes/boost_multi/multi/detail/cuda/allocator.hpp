#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ -std=c++17 -I/usr/include/cuda -Wfatal-errors -D_TEST_BOOST_MULTI_DETAIL_MEMORY_CUDA_ALLOCATOR $0x.cpp -o $0x.x -lcudart && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif

#include<cuda_runtime.h> // cudaMalloc

#include<new> // bad_alloc
#include<cassert>
#include "../cuda/ptr.hpp"

namespace boost{
namespace multi::detail::memory{

namespace cuda{

struct bad_alloc : std::bad_alloc{};

struct allocation_counter{
	static long n_allocations;
	static long n_deallocations;
	static long bytes_allocated;
	static long bytes_deallocated;
};

long allocation_counter::n_allocations = 0;
long allocation_counter::n_deallocations = 0;
long allocation_counter::bytes_allocated = 0;
long allocation_counter::bytes_deallocated = 0;

template<class T=void> class allocator : allocation_counter{
public:
	using value_type = T;
	using pointer = ptr<T>;
	using size_type = ::size_t; // as specified by CudaMalloc
	pointer allocate(size_type n, const void* = 0){
		T* p;
		cudaError_t s = cudaMalloc(&p, n*sizeof(T));
		if(s != cudaSuccess or p==nullptr) throw bad_alloc{};
		++n_allocations;
		bytes_allocated+=sizeof(T)*n;
		return {p};
	}
	void deallocate(pointer p, size_type n){
		cudaError_t s = cudaFree(p.impl_); assert(s == cudaSuccess);
		++n_deallocations;
		bytes_deallocated+=sizeof(T)*n;
	}
	std::true_type operator==(allocator const&) const{return {};}
	std::false_type operator!=(allocator const&) const{return {};}
	template<class P, class... Args>
	void construct([[maybe_unused]] P p, Args&&... args){
	//	if constexpr(sizeof...(Args) == 0 and std::is_trivially_default_constructible<T>{}){
		//	cudaError_t s = cudaMemset(p.impl_, 0, sizeof(T)); assert( s == cudaSuccess );
	//	}else{
	//	if constexpr( std::is_same<T, max_align_t>{} ) return;
		if(sizeof...(Args) == 0 and std::is_trivially_default_constructible<T>{}){
			cudaError_t s = cudaMemset(p.impl_, 0, sizeof(T)); assert( s == cudaSuccess );
		}else{
			char buff[sizeof(T)];
			::new(buff) T{std::forward<Args>(args)...};
			cudaError_t s = cudaMemcpy(p.impl_, buff, sizeof(T), cudaMemcpyHostToDevice); assert( s == cudaSuccess );
		}
	}
	template<class P>
	void destroy(P p){
		if(not std::is_trivially_destructible<T>{}){
			char buff[sizeof(T)];
			cudaMemcpy(buff, p.impl_, sizeof(T), cudaMemcpyDeviceToHost);
			((T*)buff)->~T();
		}
	}
};

}}}

#ifdef _TEST_BOOST_MULTI_DETAIL_MEMORY_CUDA_ALLOCATOR

#include<memory>
#include<iostream>
#include "../../../multi/array.hpp"

namespace boost{
namespace multi::cuda{
	template<class T, dimensionality_type D>
	using array = multi::array<T, D, multi::detail::memory::cuda::allocator<T>>;
}
}

namespace multi = boost::multi;
namespace cuda = multi::detail::memory::cuda;

void add_one(double& d){d += 1.;}
template<class T>
void add_one(T&& t){std::forward<T>(t) += 1.;}

using std::cout;

int main(){
//	std::allocator_traits<cuda::allocator<double>>::difference_type d = 0;
//	std::allocator_traits<cuda::allocator<double>>::pointer dd;
	static_assert(
		std::is_same<
			std::allocator_traits<cuda::allocator<int>>::rebind_alloc<double>,
			cuda::allocator<double>
		>{}
	);
	{
	multi::array<double, 1> arr(multi::array<double, 1>::extensions_type{100}, 0.); assert(size(arr) == 100);
	cuda::allocator<double> calloc;
	assert(calloc == calloc);
	cuda::ptr<double> p = calloc.allocate(100);
	p[33] = 123.;
	p[99] = 321.;
//	p[33] += 1;
	add_one(p[33]);
	double p33 = p[33];
	assert( p33 == 124. );
	assert( p[33] == 124. );
	assert( p[33] == p[33] );
	swap(p[33], p[99]);
	assert( p[99] == 124. );
	assert( p[33] == 321. );
	std::cout << p[33] << std::endl;
	calloc.deallocate(p, 100);

	multi::array<double, 1, cuda::allocator<double>> arr2(multi::array<double, 1>::extensions_type{100l}, 999.);
		
	assert(size(arr2) == 100);
	}
	cout<<"n_alloc/dealloc "<< cuda::allocation_counter::n_allocations <<"/"<< cuda::allocation_counter::n_deallocations <<"\n"
		<<"bytes_alloc/dealloc "<< cuda::allocation_counter::bytes_allocated <<"/"<< cuda::allocation_counter::bytes_deallocated <<"\n";

}
#endif

