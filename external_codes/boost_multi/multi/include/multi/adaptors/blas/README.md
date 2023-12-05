<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.]Multi BLAS Adaptor

(not an official Boost library)

_© Alfredo A. Correa, 2018-2021_

The BLAS Adaptor provides an interface for BLAS-like libraries.

## Contents
[[_TOC_]]

## Numeric Arrays, Conjugation Real and Imaginary parts

This functions produce views (not copies) related to conjugation, real and imaginary parts.

```cpp
	using complex = std::complex<double>; 
	complex const I{0, 1};
	multi::array<complex, 2> B = {
		{1. - 3.*I, 6. + 2.*I},
		{8. + 2.*I, 2. + 4.*I},
		{2. - 1.*I, 1. + 1.*I}
	};

	namespace blas = multi::blas;
	multi::array<complex, 2> conjB = blas::conj(B);

	assert( blas::conj(B)[2][1] == std::conj(B[2][1]) );

	assert( blas::transposed(B)[1][2] == B[2][1] );
	assert( blas::transposed(B) == ~B );

	assert( blas::hermitized(B)[2][1] == blas::conj(B)[1][2] );
	assert( blas::hermitized(B)       == blas::conj(blas::transposed(B)) );

	assert( blas::real(B)[2][1] == std::real(B[2][1]) );
	assert( blas::imag(B)[2][1] == std::imag(B[2][1]) );

	multi::array<double, 2> B_real_doubled = {
		{ 1., -3., 6., 2.},
		{ 8.,  2., 2., 4.},
		{ 2., -1., 1., 1.}
	};
	assert( blas::real_doubled(B) == B_real_doubled );
```

Usage:
```cpp
	multi::array<double, 2> const a_real = {
		{ 1., 3., 1.},
		{ 9., 7., 1.},
	};

	multi::array<complex, 2> const b = {
		{ 11.+1.*I, 12.+1.*I, 4.+1.*I, 8.-2.*I},
		{  7.+8.*I, 19.-2.*I, 2.+1.*I, 7.+1.*I},
		{  5.+1.*I,  3.-1.*I, 3.+8.*I, 1.+1.*I}
	};

	multi::array<complex, 2> c({2, 4});

	blas::real_doubled(c) = blas::gemm(1., a_real, blas::real_doubled(b)); // c = a_real*b
```

## Installation and Tests

...

