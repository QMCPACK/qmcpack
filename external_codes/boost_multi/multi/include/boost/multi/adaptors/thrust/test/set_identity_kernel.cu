// Copyright 2020-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/adaptors/thrust.hpp>

#include <thrust/complex.h>
#include <thrust/system/cuda/memory.h>

#include <chrono>

namespace multi = boost::multi;

#define CUDA_CHECKED(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line, bool /*abort*/=true) {
   if (code != cudaSuccess) {
        std::cerr<<"error: "<< cudaGetErrorString(code) <<" "<< file <<":"<< line <<std::endl;
        assert(0);
    }
}

#define REQUIRE(ans) { require((ans), __FILE__, __LINE__); }
inline void require(bool code, const char *file, int line, bool /*abort*/=true) {
   if (!code) {
        std::cerr<<"error: "<< file <<":"<< line <<std::endl;
        exit(666);
    }
}

template<typename Array2DCursor>
__global__ void kernel_setIdentity(Array2DCursor home, int m, int n) {
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = threadIdx.y + blockDim.y * blockIdx.y;
    if ((i < m) && (j < n)) {
        if (i == j) {
            home[i][j] = 1.0;
        } else {
            home[i][j] = 0.0;
        }
    }
}

template<class Array2D>
auto set_identity(Array2D&& arr) -> Array2D&&{
    int xblock_dim = 16;
    auto [m, n] = arr.sizes();
    auto xgrid_dim  = static_cast<int>((m + xblock_dim - 1) / xblock_dim);
    auto ygrid_dim  = static_cast<int>((n + xblock_dim - 1) / xblock_dim);
    dim3 block_dim(xblock_dim, xblock_dim);
    dim3 grid_dim(xgrid_dim, ygrid_dim);
    kernel_setIdentity<<<grid_dim, block_dim>>>(arr.home(), static_cast<int>(m), static_cast<int>(n));
    CUDA_CHECKED(cudaGetLastError());
    // CUDA_CHECKED(cudaDeviceSynchronize());
    return std::forward<Array2D>(arr);
}

int main() {
    using T = thrust::complex<double>;

    {
        multi::array<double, 2, thrust::cuda::allocator<double>> A({10000, 10000});
        auto const size = A.num_elements()*sizeof(T)/1e9;
        std::cout<<"size is "<< size << "GB\n";

        auto start_time = std::chrono::high_resolution_clock::now();

        thrust::fill(A.elements().begin(), A.elements().end(), 0.0);
        thrust::fill(A.diagonal().begin(), A.diagonal().end(), 1.0);

        REQUIRE( A[0][0] == 1.0 );
        REQUIRE( A[1][1] == 1.0 );
        REQUIRE( A[2][1] == 0.0 );

        auto time = std::chrono::high_resolution_clock::now() - start_time;
        auto rate = size/time.count();
        std::cout<<"algorithm rate = "<< rate <<" GB/s (ratio = 1)\n";
    }

    {
        multi::array<double, 2, thrust::cuda::allocator<double>> A({10000, 10000});
        auto const size = A.num_elements()*sizeof(T)/1e9;
        std::cout<<"size is "<< size << "GB\n";

        auto start_time = std::chrono::high_resolution_clock::now();
        set_identity(A);

        REQUIRE( A[0][0] == 1.0 );
        REQUIRE( A[1][1] == 1.0 );
        REQUIRE( A[2][1] == 0.0 );

        auto time = std::chrono::high_resolution_clock::now() - start_time;
        auto rate = size/time.count();
        std::cout<<"kernel rate = "<< rate <<" GB/s (ratio = 1)\n";
    }
}
