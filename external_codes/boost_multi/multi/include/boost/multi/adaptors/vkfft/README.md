<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.]MultiAdaptors.VkFFT

_© Alfredo A. Correa, 2026_

`Multi` is a modern C++ library that provides access and manipulation of data in multidimensional arrays.
Algorithms on multidimensional array data structures are fundamental to several branches of computing.
Multiple libraries implement these algorithms, and some are specially tuned to specific systems and hardware.

## Contents
[[_TOC_]]

## VkFFT

[VkFFT](https://github.com/DTolm/VkFFT) is a GPU library for computing the discrete Fourier transform (DFT) in one or more dimensions.
Unlike cuFFT, it is not tied to a single vendor: the same source targets Vulkan, CUDA, HIP, OpenCL, Level Zero and Metal, selected at compile time.
VkFFT is header-only and generates and compiles its compute kernels at runtime, tailored to the exact transform requested.

This adaptor drives VkFFT through the same interface as the `cufft` and [`fftw`](../fftw/README.md) adaptors.
It currently covers the **CUDA** and **HIP** backends and interleaved `complex<double>` (C2C) transforms; `execute` accepts the memory pointer as-is (`thrust::cuda::pointer`, HIP pointer, or a raw device pointer) and unwraps it, the same way the `thrust` adaptor does.
There is no CPU backend in VkFFT; for host transforms use the `fftw` adaptor.

### Plans

Plans sample the input and output *layouts* (dimensionality `D`) and are compiled on construction.

```cpp
auto p = multi::vkfft::plan<D>({which...}, in_layout, out_layout);
```

Input and output layouts must have the same sizes.
`{which...}` is a set of (at most `D`) boolean values that determine which dimensions are transformed: `{true, true, ...}` transforms every axis, `{false, true, false, ...}` transforms the second axis only, and the all-`false` mask is rejected (it would be an element-wise copy, which VkFFT does not do).
Non-transformed axes become batch dimensions; VkFFT expresses these natively (`omitDimension`), including a non-transformed axis *between* two transformed ones, which the cuFFT adaptor cannot do in a single plan.

The only constraint on the layout is that the fastest-varying axis must be contiguous (unit stride) — it need not be one of the transformed axes.
A transposed or `rotated` view where the unit-stride axis is not last is handled by sorting the axes; a view whose innermost axis is strided is not.

Plans are executed (many times if needed) with the direction as a third argument:

```cpp
p.execute(in_base, out_base, multi::vkfft::forward);   // or ::backward
```

`multi::vkfft::forward` has exponent sign −1 and `multi::vkfft::backward` sign +1, the same convention as FFTW and cuFFT.
The inverse transform is unnormalized (divide by the number of transformed elements to recover the input).
Passing the same device pointer as `in_base` and `out_base` performs the transform in place.

### Functions

There are convenience functions that build (or reuse) a plan and execute it:

```cpp
template<class In, class Out>
auto&& multi::vkfft::dft({which...}, In const& in, Out&& out, multi::vkfft::sign dir) {
	multi::vkfft::cached_plan<D>({which...}, in.layout(), out.layout()).execute(in.base(), out.base(), dir);
	return std::forward<Out>(out);
}

multi::vkfft::dft_forward ({which...}, in, out);   // dir = forward
multi::vkfft::dft_backward({which...}, in, out);   // dir = backward
auto out = multi::vkfft::dft({which...}, in, dir); // allocates the result
```

Because a plan compiles VkFFT kernels at construction (via NVRTC on CUDA, HIPRTC on HIP), the first use of a given transform is slow.
`multi::vkfft::cached_plan<D>` keeps a thread-local cache keyed on `{which, in-layout, out-layout}`, so repeated calls with the same shape pay the compilation cost only once.
The plan also allocates GPU scratch (twiddle LUTs, and a temp buffer for large transforms) at construction; it is released when the plan is destroyed, so a long-lived `cached_plan` keeps that memory for every distinct shape it has seen.

### Building

The adaptor needs the VkFFT headers and, for the CUDA backend, the CUDA **driver API** and **NVRTC** (HIP needs **HIPRTC**).
The CMake target `multi::vkfft` bundles all of this:

```cmake
find_package(VkFFT)   # FindVkFFT.cmake, or set -DVkFFT_ROOT=<prefix>; falls back to fetching VkFFT
target_link_libraries(my_target PRIVATE multi::vkfft)
```

`VKFFT_BACKEND` is selected automatically from `__CUDACC__` / `__HIPCC__` when the translation unit is compiled as CUDA or HIP; `multi::vkfft` also sets it explicitly (`1` for CUDA, `2` for HIP).
`VKFFT_MAX_FFT_DIMENSIONS` (default 12 in the adaptor header) bounds the plan dimensionality; raise it with `-DVKFFT_MAX_FFT_DIMENSIONS=N` if you need transforms of more dimensions.
