#ifdef COMPILATION
/usr/local/cuda/bin/nvcc -DBOOST_ALL_NO_LIB -DBOOST_SERIALIZATION_DYN_LINK -DUSING_Libxc -I/home/correaa/prj/inq/external_libs/libxc/src -I/home/correaa/prj/inq/build/cuda11/external_libs/libxc -I/home/correaa/prj/inq/build/cuda11/external_libs/libxc/gen_funcidx -I/home/correaa/prj/inq/build/cuda11 -I/home/correaa/prj/inq/build/cuda11/external_libs/pseudopod -I/home/correaa/prj/inq/src/. -I/home/correaa/prj/inq/src/../external_libs  -DNDEBUG -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include  -g -pg -D_DISABLE_CUDA_SLOW -O3 --expt-relaxed-constexpr --expt-extended-lambda --Werror=cross-execution-space-call --compiler-options -Ofast,-g,-pg,-std=c++14,-Wall,-Wfatal-errors -g   -x cu $0 -o $0x&&$0x&&rm $0x;exit
#/usr/local/cuda/bin/nvcc -g -pg -O3 --expt-relaxed-constexpr --expt-extended-lambda --Werror=cross-execution-space-call --compiler-options -Ofast,-g,-pg,-Wall,-Wfatal-errors -g $0 - $0x &&$0x&&rm $0x; exit
#endif

#include "../../cuda/managed/allocator.hpp"
#include "../../../../array.hpp"

#include <cuda.h>

template<class Kernel>
__global__ void cuda_run_kernel_1(unsigned n, Kernel k){
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

#include<thrust/complex.h>

namespace multi = boost::multi;
namespace cuda = multi::memory::cuda;

using complex = thrust::complex<double>;

int main(){


	int set_size = 32;
	int basis_size = 100;
	
	multi::array<complex, 2, cuda::managed::allocator<complex>> phi1({set_size, basis_size}, 1.);
	multi::array<complex, 2, cuda::managed::allocator<complex>> phi2({set_size, basis_size}, 1.);
	
	multi::array<complex, 1, cuda::managed::allocator<complex>> overlap(set_size);

#if 0 // using iterators	
	auto phi1_p = begin(phi1);
	auto phi2_p = begin(phi2);
	
	auto overlap_p = begin(overlap);

	run(set_size, 
		[phi1_p, phi2_p, overlap_p, basis_size] __device__ (int ist){
			complex a = 0.;
			for(int ip = 0; ip < basis_size; ip++){
				auto p1 = phi1_p[ist][ip];
				auto p2 = phi2_p[ist][ip];
				a += conj(p1)*p2;
			}
			overlap_p[ist] = a;
		}
	);
#endif
#if 0
	run(set_size, 
		[phi1_p = begin(phi1), phi2_p = begin(phi2), overlap_p = begin(overlap), basis_size] __device__ (int ist){
			complex a = 0.;
			for(int ip = 0; ip < basis_size; ip++){
				auto p1 = phi1_p[ist][ip];
				auto p2 = phi2_p[ist][ip];
				a += conj(p1)*p2;
			}
			overlap_p[ist] = a;
		}
	);
#endif

	run(set_size, 
		[phi1_p = &phi1(), phi2_p = &phi2(), overlap_p = &overlap()] __device__ (int ist){
			(*overlap_p)[ist] = 0;
			for(int ip = 0; ip < overlap_p->size(); ip++)
				(*overlap_p)[ist] += conj((*phi1_p)[ist][ip])*(*phi2_p)[ist][ip];
		}
	);

	assert( overlap[21] == double(basis_size) );
	
#if 0	
	{
		auto npoints = phi1.basis().part().local_size();
		auto vol_element = phi1.basis().volume_element();
		auto phi1p = begin(phi1.matrix());
		auto phi2p = begin(phi2.matrix());
		auto overlap = begin(overlap_vector);
			
		//OPTIMIZATION: here we should parallelize over points as well 
		gpu::run(phi1.set_size(),
						 [phi1p, phi2p, overlap, npoints,vol_element] __host__ __device__ (int ist){
							 type aa = 0.0;
							 for(int ip = 0; ip < npoints; ip++){
								phi1p[ip];
							//	 auto p1 = phi1p[ip][ist];
							//	 auto p2 = phi2p[ip][ist];
							//	 aa += conj(p1)*p2;
									 
							 }
								 
							 overlap[ist] = vol_element*aa;
						 });
	}
#endif
	return 0;
}

