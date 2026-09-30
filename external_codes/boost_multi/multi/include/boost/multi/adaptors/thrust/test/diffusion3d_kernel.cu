#define ENABLE_GPU 1

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/adaptors/thrust/reduce_by_index.hpp>
#include <boost/multi/array.hpp>

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/reduce.h>
#include <thrust/universal_allocator.h>

#include <boost/core/lightweight_test.hpp>

#include <chrono>

template<class Tp>
inline
#if defined(_MSC_VER)
	__forceinline
#else
	__attribute__((always_inline))
#endif
	void
	DoNotOptimize(Tp const& value) {  // NOLINT(readability-identifier-naming)
#if defined(_MSC_VER)
	_ReadWriteBarrier();
	(void)value;
#else
	asm volatile("" : : "r,m"(value) : "memory");  // NOLINT(hicpp-no-assembler)
#endif
}

// A simple macro for checking CUDA API calls
#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
void check(cudaError_t result, char const* const func, char const* const file, int const line) {
	if(result) {
		fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line, static_cast<unsigned int>(result), cudaGetErrorName(result), func);
		exit(99);
	}
}

// The 3D diffusion kernel to be tested.
__global__ void diffusion3d_gpu(int nx, int ny, int nz, float cc, float ce, float cw, float cn, float cs, float ct, float cb, float const* f, float* fn) {
	int const i = blockIdx.x * blockDim.x + threadIdx.x;
	int const j = blockIdx.y * blockDim.y + threadIdx.y;
	int const k = blockIdx.z * blockDim.z + threadIdx.z;

	if(i < nx && j < ny && k < nz) {
		int const ix = nx * ny * k + nx * j + i;
		int const ip = i == nx - 1 ? ix : ix + 1;
		int const im = i == 0 ? ix : ix - 1;
		int const jp = j == ny - 1 ? ix : ix + nx;
		int const jm = j == 0 ? ix : ix - nx;
		int const kp = k == nz - 1 ? ix : ix + nx * ny;
		int const km = k == 0 ? ix : ix - nx * ny;

		fn[ix] = cc * f[ix] + ce * f[ip] + cw * f[im] + cn * f[jp] + cs * f[jm] + ct * f[kp] + cb * f[km];
	}
}

// A CPU version of the function for verification purposes.
void diffusion3d_cpu(int nx, int ny, int nz, float cc, float ce, float cw, float cn, float cs, float ct, float cb, float const* f, float* fn_cpu) {
	for(int k = 0; k < nz; k++) {
		for(int j = 0; j < ny; j++) {
			for(int i = 0; i < nx; i++) {
				int const ix = nx * ny * k + nx * j + i;
				int const ip = i == nx - 1 ? ix : ix + 1;
				int const im = i == 0 ? ix : ix - 1;
				int const jp = j == ny - 1 ? ix : ix + nx;
				int const jm = j == 0 ? ix : ix - nx;
				int const kp = k == nz - 1 ? ix : ix + nx * ny;
				int const km = k == 0 ? ix : ix - nx * ny;

				fn_cpu[ix] = cc * f[ix] + ce * f[ip] + cw * f[im] + cn * f[jp] + cs * f[jm] + ct * f[kp] + cb * f[km];
			}
		}
	}
}

auto main() -> int {
    // Problem dimensions
    const int nx = 128;
    const int ny = 128;
    const int nz = 128;
    {
        const int arraySize = nx * ny * nz;
        const int memSize = arraySize * sizeof(float);

        // Host pointers
        float *f_h, *fn_h, *fn_cpu;

        // Device pointers
        float *f_d, *fn_d;

        // Allocate host memory
        f_h = (float*)malloc(memSize);
        fn_h = (float*)malloc(memSize);
        fn_cpu = (float*)malloc(memSize);

        // Initialize host input data
        for (int i = 0; i < arraySize; i++) {
            f_h[i] = (float)i;
        }

        // Allocate device memory and copy data from host to device
        checkCudaErrors(cudaMalloc((void**)&f_d, memSize));
        checkCudaErrors(cudaMalloc((void**)&fn_d, memSize));
        checkCudaErrors(cudaMemcpy(f_d, f_h, memSize, cudaMemcpyHostToDevice));

        // Diffusion coefficients
        float cc = 0.5f, ce = 0.1f, cw = 0.1f, cn = 0.1f, cs = 0.1f, ct = 0.1f, cb = 0.1f;

        // Define CUDA grid and block dimensions
        const int threadsPerBlock_x = 16;
        const int threadsPerBlock_y = 16;
        const int threadsPerBlock_z = 4;
        dim3 threadsPerBlock(threadsPerBlock_x, threadsPerBlock_y, threadsPerBlock_z);
        dim3 numBlocks((nx + threadsPerBlock_x - 1) / threadsPerBlock_x,
                    (ny + threadsPerBlock_y - 1) / threadsPerBlock_y,
                    (nz + threadsPerBlock_z - 1) / threadsPerBlock_z);

        // Launch the CUDA kernel
        diffusion3d_gpu<<<numBlocks, threadsPerBlock>>>(nx, ny, nz, cc, ce, cw, cn, cs, ct, cb, f_d, fn_d);
        checkCudaErrors(cudaGetLastError());

        // Copy the result back from device to host
        checkCudaErrors(cudaMemcpy(fn_h, fn_d, memSize, cudaMemcpyDeviceToHost));

        // Perform CPU calculation for verification
        diffusion3d_cpu(nx, ny, nz, cc, ce, cw, cn, cs, ct, cb, f_h, fn_cpu);

        // Compare results
        int errors = 0;
        for (int i = 0; i < arraySize; i++) {
            if (abs((fn_h[i] - fn_cpu[i])/fn_cpu[i]) > 1e-6) {
                if (errors < 10) {
                    printf("Error at index %d: GPU result %f, CPU result %f\n", i, fn_h[i], fn_cpu[i]);
                }
                errors++;
            }
        }

        if (errors == 0) {
            printf("Success! GPU and CPU results match.\n");
        } else {
            printf("Verification failed: Found %d errors.\n", errors);
            return 1;
        }

        // Free all allocated memory
        free(f_h);
        free(fn_h);
        free(fn_cpu);
        checkCudaErrors(cudaFree(f_d));
        checkCudaErrors(cudaFree(fn_d));
    }
    {
        namespace multi = boost::multi;
        multi::thrust::device_array<float, 3> f_d({nx, ny, nz});
        // f_d.elements() = [] __host__ __device__ (multi::thrust::device_array<float, 3>::index i) { return i; } ^ multi::extents_t<1>(f_d.num_elements());
    }
    return 0;
}
