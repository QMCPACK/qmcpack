<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.]MultiAdaptors.FFTW

> **Disclosure: This is not an official or accepted Boost library and is unrelated to the std::mspan proposal.**

_© Alfredo A. Correa, 2018-2023_

`Multi` is a modern C++ library that provides access and manipulation of data in multidimensional arrays.
Algorithms on multidimensional array data structures are fundamental to several branches of computing.
Multiple libraries implement these algorithms, and some are specially tuned to specific systems and hardware.

Linear algebra and Fourier transforms are some examples of operations with algorithms on regularly contiguous (strided) multidimensional array data structures.
Although not generic, these libraries are the best options in specific systems for specific element types.

## Contents
[[_TOC_]]

## FFTW

FFTW is a C library for computing the discrete Fourier transform (DFT) in one or more dimensions for real and complex data.
It is the defacto interface for many implementations, including Intel's MKL FFT.

The FFTW adaptor provides two ways to use the library: one is through plan objects, and the other is through functions.

Plans are runtime-optimized algorithms that tune the DFT operation for specific array sizes and layouts known in advance.
Plans reserve resources and precalculate parameters utilized during the execution.

Plans are created from array layouts with dimensionality `D` that sample the input and output.

```cpp
auto p = multi::fftw::plan::[forward|backward]({which...}, in_layout, out_layout);
```

Input and output layout must have the same associated sizes.
`{which...}` is a set of (at most D) boolean value that determined which dimensions are transformed; for example `{true, true, ...}` performs the FFT on all directions, ``{false, true, false, ...}` for the second dimension only and `{false, false, ...}` doesn't perform any Fourier transform, effectively performing a element-wise copy or transposition.

The plans can be later executed (many times if necessary) as:
```cpp
p.execute(in_base, out_base);
```

Executions of the same plan (or over the same data) are not thread-safe (the plan has internal buffers that are modified), so the `.execute` function is not marked `const`.

The use pattern of the FFTW adaptor (and the original FFTW) interface is somewhat entangled.
The plan construction takes the arrays, and the execute takes the internal pointers to array data.

There is a convenience function that generates and executes the plan consistently:
```cpp
template<class In, class Out>
auto&& multi::fftw::dft::forward({which, ...}, In const& in, Out&& out) {
	multi::fftw::plan::forward(in.layout(), out.layout()).execute(in.base(), out.base());
	return std::forward<Out>(out);
}
```

<!-- Finally, FFTW offers plans creation that are fast at the expense of complexity.

```cpp
auto p = multi::fftw::measured_plan::[forward|backward]({which...}, in_buffer, in_layout, out_buffer, out_layout);
```

These plans require passing memory pointers that is possibly overwritten and are more complex to use.

The buffer must be (at least) the "hull" size of the respective layouts.
This is typically satisfied by the (sub)arrays base pointers themselves, however, output will be overwritten.

```cpp
auto p = multi::fftw::measured_plan::[forward|backward]({which...}, in.base(), in.layout(), out.base(), out.layout());
```

This type of plans are harder to use.
 -->
