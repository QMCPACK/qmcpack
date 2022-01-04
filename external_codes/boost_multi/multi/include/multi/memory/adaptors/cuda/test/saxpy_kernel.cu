#ifdef COMPILATION
/usr/local/cuda/bin/nvcc -g -pg -O3 --expt-relaxed-constexpr --extended-lambda $0 -o $0x&&time $0x&&rm $0x;exit
#/usr/local/cuda/bin/nvcc -g -pg -O3 --expt-relaxed-constexpr --extended-lambda --expt-extended-lambda --Werror=cross-execution-space-call --compiler-options -Ofast,-g,-pg,-Wall,-Wfatal-errors -g $0 - $0x &&$0x&&rm $0x; exit
#endif

#include "../../cuda/managed/allocator.hpp"
#include "../../../../array.hpp"

#include <cuda.h>

#include <stdio.h>

namespace boost::multi{

template<class Kernel>
__global__ void cuda_run_kernel_1(unsigned const n, Kernel k){
	auto i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i < n) k(i);
}

template<class Kernel>
void run(std::size_t n, Kernel k){
	constexpr int CUDA_BLOCK_SIZE = 256;
	unsigned nblock = (n + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;
	cuda_run_kernel_1<<<nblock, CUDA_BLOCK_SIZE>>>(n, k);
	auto code = cudaGetLastError();
	if(code != cudaSuccess){
		fprintf(stderr,"GPU assert: %s %s %d\n", cudaGetErrorString(code), __FILE__, __LINE__);
		exit(-1);
	}
	cudaDeviceSynchronize();
}


template <class kernel_type>
__global__ void cuda_run_kernel_2(unsigned sizex, unsigned sizey, unsigned dim2, kernel_type kernel){
	auto i1 = blockIdx.x*blockDim.x + threadIdx.x;
	auto i2 = blockIdx.y*blockDim.y + threadIdx.y;
	auto i3 = blockIdx.z*blockDim.z + threadIdx.z;
	
	auto ix = i1;
	auto iy = i2 + dim2*i3;
	if(ix < sizex && iy < sizey) kernel(ix, iy);
}

//finds fact1, fact2 < thres such that fact1*fact2 >= val
inline void factorize(const std::size_t val, const std::size_t thres, std::size_t & fact1, std::size_t & fact2){
	fact1 = val;
	fact2 = 1;
	while (fact1 > thres){
		fact1 = (fact1 + 1)/2;
		fact2 *= 2;
	}

	assert(fact1*fact2 >= val);
}

template <class kernel_type>
void run(size_t sizex, size_t sizey, kernel_type kernel){
	constexpr int CUDA_BLOCK_SIZE = 256;
	unsigned nblock = (sizex + CUDA_BLOCK_SIZE - 1)/CUDA_BLOCK_SIZE;

	constexpr int CUDA_MAX_DIM23 = 65535;

	size_t dim2, dim3;
	factorize(sizey, CUDA_MAX_DIM23, dim2, dim3);
	
	struct dim3 dg{nblock, unsigned(dim2), unsigned(dim3)};
	struct dim3 db{CUDA_BLOCK_SIZE, 1, 1};
	cuda_run_kernel_2<<<dg, db>>>(sizex, sizey, dim2, kernel);
	
	assert(cudaGetLastError() == cudaError_t(CUDA_SUCCESS));
		
	cudaDeviceSynchronize();
}

}

namespace multi = boost::multi;

#include "../../../../adaptors/cuda.hpp"
#include "../../../../array.hpp"
#include<numeric>
#include<thrust/complex.h>

template<class T> constexpr void what(T&&) = delete;

int main(){
	using T = thrust::complex<float>;

	{
		int N = 1<<20; // 1<<20 == 1048576

		multi::array<T, 1> x(N); std::iota(begin(x), end(x), 1.f);
		multi::array<T, 1> y(N); std::iota(begin(x), end(x), 2.f);

		multi::cuda::managed::array<T, 1> Dx = x;
		multi::cuda::managed::array<T, 1> Dy = y;

		assert( Dx.size() == Dy.size() );

		float a = 2.f;

	//	for(int i = 0; i < Dx.size(); ++i)
	//		Dy[i] = a*Dx[i] + Dy[i];

		multi::run(
			Dx.size(), 
			[a, Dy=begin(Dy), Dx=begin(Dx)] __device__ (int i){
				Dy[i] = a*Dx[i] + Dy[i];
			}
		);
/*
		multi::run(
			Dx.size(), 
			[a, DyP=&Dy(), DxP=&Dx()] __device__ (int i){
				(*DyP)[i] = a*(*DxP)[i] + (*DyP)[i];
			}
		);
*/
		y = Dy;

		assert( y[0] == 4.0f );
	}
	
	{
		int N = 1000;
		multi::array<T, 2> a({N, N}); 
		std::iota(a.elements().begin(), a.elements().end(), 1.f);

		multi::array<T, 2> b({N, N}); 
		std::iota(b.elements().begin(), b.elements().end(), 2.f);

		multi::cuda::managed::array<T, 2> const Da = a;
		multi::cuda::managed::array<T, 2> Db = b;

		multi::run(
			N, N,
			[Da=begin(Da), Db=begin(Db)] __device__ (int i, int j){
				Db[i][j] = Da[i][j] + Db[i][j];
			}
		);
		
		multi::array<T, 2> b_cpy = Db;
		assert( b_cpy[11][23] == a[11][23] + b[11][23] );
	}

}

