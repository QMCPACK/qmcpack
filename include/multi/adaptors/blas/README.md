<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# Multi BLAS Adaptor

_© Alfredo A. Correa, 2018-2023_

The BLAS Adaptor provides an interface for BLAS and BLAS-like libraries (namely cuBLAS).

## Contents
[[_TOC_]]

## Numeric Arrays, Conjugation Real and Imaginary parts

These functions produce views (not copies) related to conjugation, real and imaginary parts.

```cpp
	using complex = std::complex<double>; 
	complex const I{0.0, 1.0};
	multi::array<complex, 2> B = {
		{1.0 - 3.0*I, 6.0 + 2.0*I},
		{8.0 + 2.0*I, 2.0 + 4.0*I},
		{2.0 - 1.0*I, 1.0 + 1.0*I}
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
		{ 1.0, -3.0, 6.0, 2.0},
		{ 8.0,  2.0, 2.0, 4.0},
		{ 2.0, -1.0, 1.0, 1.0}
	};
	assert( blas::real_doubled(B) == B_real_doubled );
```

```cpp
	multi::array<double, 2> const a_real = {
		{ 1.0, 3.0, 1.0},
		{ 9.0, 7.0, 1.0},
	};

	multi::array<complex, 2> const b = {
		{ 11.0 + 1.0*I, 12.0 + 1.0*I, 4.0 + 1.0*I, 8.0 - 2.0*I},
		{  7.0 + 8.0*I, 19.0 - 2.0*I, 2.0 + 1.0*I, 7.0 + 1.0*I},
		{  5.0 + 1.0*I,  3.0 - 1.0*I, 3.0 + 8.0*I, 1.0 + 1.0*I}
	};

	multi::array<complex, 2> c({2, 4});

	blas::real_doubled(c) = blas::gemm(1., a_real, blas::real_doubled(b)); // c = a_real*b
```

## Table of features

All these operations are now supported for CPU and GPU memory, real and complex.

scalars: `aa` ($`\alpha`$), `bb` ($`\beta`$) \
vectors: `x`, `y` \
matrices: `A`, `B`, `C`

vector operations: `C` (`*`) conjugation (element-wise) \
matrix operations: `J` (`*`) conjugation (element-wise) (use `C` for vectors), `T` transpose, `H` transpose conjugate (also `C`, discouraged), `U`/`L` upper or lower triangular part (logical zeroing other side)


| BLAS   | mutable form           | effect                        | operator form [³]        | functional form | thrust/STL [¹] |
|---     |---                     | ---                           | ---                  | ---             | --- |
| SWAP   |`blas::swap(x, y)`      | $`x_i \leftrightarrow y_i`$ | `(x^y)` |    | `swap_ranges(begin(x), end(x), begin(y))` |
| COPY   |`blas::copy(x, y)`      | $`y_i \leftrightarrow x_i`$ | `y << x` |  `y = blas::copy(x)` | `copy(begin(x), end(x), begin(y))` |
| ASUM   |`blas::asum(x, res)`    | $`r \leftarrow \sum_i \|\Re x_i\| + \|\Im x_i\|`$ | `x==0`/`x!=0` `isinf(x)` `isnan(x)`[²] | `res = blas::asum(x)` | `transform_reduce(begin(x), end(x), 0.0, plus<>{}, [](auto const& e){return abs(e.real()) + abs(e.imag());})` |
| NRM2   |`blas::nrm2(x, res)`     | $`r \leftarrow \sqrt{\sum_i \|x_i\|^2}`$ | `abs(x)` | `res = blas::nrm2(x);` | `sqrt(trasnform_reduce(begin(x), end(x), 0.0, plus<>{}, [](auto const& e){return norm(e);}));` |
| SCAL   |`blas::scal(aa, x);`    | $`x_i \leftarrow \alpha x_i`$ | `x*=aa;`           |                 | `for_each(begin(x), end(x), [aa](auto& e){return e*=aa;})` |
| AXPY   |`blas::axpy(aa, x, y)`  | $`y_i \leftarrow \alpha x_i + y_i`$ | `y+=x` `y-=x` `y+=aa*x` `y-=aa*x` |     | `transform(x.begin(), x.end(), y.begin(), y.begin(), [aa](auto ex, auto ey) {return aa*ex + ey;}` |
| DOT    | `blas::dot(x, y, res)` | $`r = \sum_i x_i y_i`$        | `res = (x, y);`       | `res = blas::dot(x, y)`                | `inner_product(begin(x), end(x), begin(y), T{});` |
|        | `blas::dot(blas::C(x), y, res)` |  $`r = \sum_i \bar x_i y_i`$ | `res = (*x, y);`  | `res = blas::dot(blas::C(x), y)`             | `inner_product(begin(x), end(x), begin(y), T{}, plus<>{}, [](T const& t1, T const& t2) {return conj(t1)*t2;});` |
|        | `blas::dot(x, blas::C(y), res)` |  $`r = \sum_i x_i \bar y_i`$ | `res = (x, *y);`  | `res = blas::dot(x, blas::C(y));`             | `inner_product(x.begin(), x.end(), y.begin(), T{}, plus<>{}, [](T const& t1, T const& t2) {return t1*conj(t2);});` |
|        | ~~`blas::dot(blas::C(x), blas::C(y), res)`~~ | $`r = \sum_i \bar x_i \bar y_i`$ not implemented in BLAS, conjugate result |  | | `auto res = conj(inner_product(x.begin(), x.end(), y.begin(), T{});` |
| GEMV   | `blas::gemv(aa, A, x, bb, y)` | $`y_i \leftarrow \alpha\sum_j A_{ij}x_j + \beta y_i`$ | `y=A%x` `y=aa*A%x` `y+=A%x` `y+=aa*A%x`[¤] | `y=blas::gemv(aa, A, x)` `y+=blas::gemv(aa, A, x)`  | `transform(begin(A), end(A), begin(y), [&x, aa] (auto const& Ac) {return aa*blas::dot(Ac, x);})` |
|        | `blas::gemv(aa, blas::T(A), x, bb, y)` | $`y_i \leftarrow \alpha\sum_j A_{ji}x_j + \beta y_i`$ | `y= ~A % x` `y=aa*(~A)%x` `y+=(~A)%x` `y+=aa*(~A)%x` | `y=blas::gemv(aa, blas::T(A), x)` `y+=blas::gemv(aa, blas::T(A), x)`  | `transform(begin(transposed(A)), end(transposed(A)), begin(y), [&x, aa] (auto const& Ac) {return aa*blas::dot(Ac, x);})` |
|        | `blas::gemv(aa, blas::J(A), x, bb, y)` | $`y_i \leftarrow \alpha\sum_j A_{ij}^*x_j + \beta y_i`$ | `y= *A % x` `y=aa*(*A)%x` `y+=(*A)%x` `y+=aa*(*A)%x` | `y=blas::gemv(aa, blas::J(A), x)` `y+=blas::gemv(aa, blas::J(A), x)`  | `transform(begin(A), end(A), begin(y), [&x, aa] (auto const& Ac) {return aa*blas::dot(*Ac, x);})` |
|        | ~~`blas::gemv(aa, blas::H(A), x, bb, y)`~~ | $`y_i \leftarrow \alpha\sum_j A_{ji}^*x_j + \beta y_i`$ (not BLAS-implemented)|    |   | `transform(begin(transposed(A)), end(transposed(A)), begin(y), [&x, aa] (auto const& Ac) {return aa*blas::dot(*Ac, x);})` |
| GEMM | `blas::gemm(aa, A, B, bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ik} B_{kj} + \beta C_{ij}`$ | `C = aa*(A*B)` | `C = blas::gemm(aa, A, B)` `C += blas::gemm(aa, A, B)` | `transform(begin(A), end(A), begin(C), begin(C), [&B, aa, bb] (auto const& Ar, auto&& Cr) {return blas::gemv(aa, blas::T(B), Ar, bb, move(Cr));})` |
|      | `blas::gemm(aa, A, blas::T(B), bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ik} B_{jk} + \beta C_{ij}`$ | `C = aa*(A* ~B)` | `C = blas::gemm(aa, A, blas::T(B))` `C += blas::gemm(aa, A, blas::T(B))` | `transform(begin(A), end(A), begin(C), begin(C), [&B, aa, bb] (auto const& Ar, auto&& Cr) {return blas::gemv(aa, B, Ar, bb, move(Cr));})` |
|        | `blas::gemm(aa, blas::T(A), B, bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ki} B_{kj} + \beta C_{ij}`$ | `C =~A * B` `C = aa*(~A * B)` `C+=~A * B` `C+=aa*(~A * B)` | `C = blas::gemm(aa, blas::T(A), B, bb, C)` (or `+=`) | `transform(begin(transposed(A)), end(transposed(A)), begin(C), begin(C), [&B, aa, bb] (auto const& Ar, auto&& Cr) {return blas::gemv(aa, blas::T(B), Ar, bb, std::move(Cr));})` |
|        | `blas::gemm(aa, blas::T(A), blas::T(B), bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ki} B_{jk} + \beta C_{ij}`$ | `C =~A * ~B` `C = aa*(~A * ~B)` `C+=~A * ~B` `C+=aa*(~A * ~B)` | `C = blas::gemm(aa, blas::T(A), blas::T(B), bb, C)` (or `+=`) | `transform(begin(transposed(A)), end(transposed(A)), begin(C), begin(C), [&B, aa, bb] (auto const& Ar, auto&& Cr) {return blas::gemv(aa, B, Ar, bb, std::move(Cr));})` |
|      | <s>`blas::gemm(aa, A, blas::J(B), bb, C)`</s> (use `blas::gemm(..., blas::T(B), blas::H(A), ..., HC)` and conjtranspose result) | $`C_{ij} \leftarrow \alpha \sum_k A_{ik} B_{kj}^* + \beta C_{ij}`$ (not BLAS-implemented) |    |   | `transform(begin(A), end(A), begin(C), begin(C), [BT=transposed(B)](auto const& Ar, auto&& Cr) {transform(begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto&& c) {return aa*blas::dot(Ar, blas::C(Bc)) + bb*c;}); return std::move(Cr);});` |
|      | ~~`blas::gemm(aa, blas::J(A), B, bb, C)`~~ | $`C_{ij} \leftarrow \alpha \sum_k A_{ik}^* B_{kj} + \beta C_{ij}`$ (not BLAS-implemented) |    |   | `transform(begin(A), end(A), begin(C), begin(C), [BT=transposed(B)](auto const& Ar, auto&& Cr) {transform(begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto&& c) {return aa*blas::dot(blas::C(Ar), Bc) + bb*c;}); return std::move(Cr);});` |
|      | <s>`blas::gemm(aa, blas::J(A), blas::J(B), bb, C)`</s> | $`C_{ij} \leftarrow \alpha \sum_k \bar{A_{ik}} \bar{B_{kj}} + \beta C_{ij}`$ (not BLAS-implemented) |    |   | `transform(begin(A), end(A), begin(C), begin(C), [BT=transposed(B)](auto const& Ar, auto&& Cr) {transform(begin(BT), end(BT), begin(Cr), begin(Cr), [&Ar](auto const& Bc, auto&& c) {return aa*blas::dot(blas::C(Ar), blas::C(Bc)) + bb*c;}); return std::move(Cr);});` |
|      | `blas::gemm(aa, A, blas::H(B), bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ik} \bar B_{jk} + \beta C_{ij}`$ | `C = aa*(A* ~*B)` (or `+=`) | `C = blas::gemm(aa, A, blas::H(B))` `C += blas::gemm(aa, A, blas::H(B))` | `transform(begin(A), end(A), begin(CC), begin(CC), [&](auto const& Ar, auto&& Cr){return blas::gemv(aa, blas::J(B), Ar, bb, move(Cr));})` |
|      | `blas::gemm(aa, blas::H(A), B, bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k \bar A_{ki} B_{kj} + \beta C_{ij}`$ | `CC=~*A *B` | `C=blas::gemm(aa, blas::H(A), B)` | `transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [BT=transposed(B)](auto const& Ac, auto&& Cr) {transform(begin(BT), end(BT), begin(Cr), begin(Cr), [&Ac](auto const& Bc, auto&& c){return aa*blas::dot(blas::C(Ac), Bc) + bb*c;}); return move(Cr);})` |
|      | `blas::gemm(aa, blas::H(A), blas::H(B), bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k \bar A_{ki} \bar B_{jk} + \beta C_{ij}`$ | `CC=~*A * ~*B` | `C=blas::gemm(aa, blas::H(A), blas::H(B))` | `transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [&B](auto const& Ac, auto&& Cr) {transform(begin(B), end(B), begin(Cr), begin(Cr), [&Ac](auto const& Bc, auto&& c) {return conj(std::transform_reduce(begin(Ac), end(Ac), begin(Bc), 0.0*c, std::plus<>{}, [](auto const& a, auto const& b) {return a*b;}));}); return move(Cr);})` |
|      | `blas::gemm(aa, blas::T(A), blas::H(B), bb, C)` | $`C_{ij} \leftarrow \alpha \sum_k A_{ki} \bar B_{jk} + \beta C_{ij}`$ | `CC=~A * ~*B` | `C=blas::gemm(aa, blas::T(A), blas::H(B))` | `transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [&B](auto const& Ac, auto&& Cr) {transform(begin(B), end(B), begin(Cr), begin(Cr), [&Ac](auto const& Bc, auto&& c) {return std::transform_reduce(begin(Ac), end(Ac), begin(Bc), 0.0*c, std::plus<>{}, [](auto const& a, auto const& b) {return a*conj(b);});}); return move(Cr);})` |
|      | ~~`blas::gemm(aa, blas::H(A), blas::T(B), bb, C)`~~ | $`C_{ij} \leftarrow \alpha \sum_k \bar A_{ki} B_{jk} + \beta C_{ij}`$ (not BLAS-implemented) |  |  | `transform(begin(transposed(A)), end(transposed(A)), begin(CC), begin(CC), [&B](auto const& Ac, auto&& Cr) {transform(begin(B), end(B), begin(Cr), begin(Cr), [&Ac](auto const& Bc, auto&& c) {return std::transform_reduce(begin(Ac), end(Ac), begin(Bc), 0.0*c, std::plus<>{}, [](auto const& a, auto const& b) {return conj(a)*b;});}); return move(Cr);})` |
|      | ~~`blas::gemm(aa, blas::J(A), blas::H(B), bb, C)`~~ | $`C_{ij} \leftarrow \alpha \sum_k \bar A_{ik} \bar B_{jk} + \beta C_{ij}`$ (not BLAS-implemented) |  |  |  |
|      | ~~`blas::gemm(aa, blas::H(A), blas::J(B), bb, C)`~~ | $`C_{ij} \leftarrow \alpha \sum_k \bar A_{ki} \bar B_{kj} + \beta C_{ij}`$ (not BLAS-implemented) |  |  |  |
| TRSM | `blas::trsm(blas::side::right, aa, blas::U(A), B)` | $`B\leftarrow B.U^{-1}`$  | `B /= U(A)` | | TODO |  
|      | `blas::trsm(blas::side::right, aa, blas::L(A), B)` | $`B\leftarrow B.L^{-1}`$  | `B /= L(A)` | | TODO |
|      | `blas::trsm(blas::side::left, aa, blas::U(A), B)` | $`B\leftarrow U^{-1}.B`$  | `B \|= U(A)` | | TODO |  
|      | `blas::trsm(blas::side::left, aa, blas::L(A), B)` | $`B\leftarrow L^{-1}.B`$  | `B \|= L(A)` | | TODO |
|      | ~~`blas::trsm(blas::side::right, aa, blas::U(A), blas::J(B))`~~ | $`B*\leftarrow B*.U^{-1}`$ $`B\leftarrow B.U*^{-1}`$ | | | TODO |
|      | ~~`blas::trsm(blas::side::right, aa, blas::L(A), blas::J(B))`~~ | $`B*\leftarrow B*.L^{-1}`$ $`B\leftarrow B.L*^{-1}`$ | | | TODO |
|      | `blas::trsm(blas::side::right, aa, blas::U(A), blas::H(B))` | $`B^\dagger\leftarrow B^\dagger.U^{-1}`$ $`B\leftarrow U^\dagger^{-1}.B`$ | | | TODO |
|      | `blas::trsm(blas::side::right, aa, blas::L(A), blas::H(B))` | $`B^\dagger\leftarrow B^\dagger.L^{-1}`$ $`B\leftarrow L^\dagger^{-1}.B`$ | | | TODO |

[¹]: for reference, not optimal. \
[²]: `asum` is interpreted as a mechanism to detect null vectors or vectors containing NaN or infinities. \
[³]: needs explicit invocation `using namespace multi::operators` namespace or of specific symbols `using multi::operator*`/`operator/=`/etc. \
[¤]: `y *=bb +=aa*A%x` (`gemv(aa, A, x, bb, y)`) would also be possible.
