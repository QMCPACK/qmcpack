<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# Multi BLAS Adaptor

_© Alfredo A. Correa, 2018-2024_

(documentation in progress)

The BLAS Adaptor provides an interface for the BLAS and BLAS-like linear algebra libraries (cuBLAS and hipBLAS).
Although BLAS is not strictly a multidimensional array library, as it only works on 1D (vectors) and 2D arrays (matrices), it is an extremely popular numeric library.

The adaptor library has a two-fold purpose:

First, it allows the abstracting of the stride information and the conjugation/transposition in the BLAS calls, simplifying the interface enormously and making it consistent with their GPU counterparts, such as cuBLAS.
"View manipulators" automatically handle cases related to conjugation, transposition, and real and imaginary parts.

Second, it provides a functional interface to the BLAS calls, which is easier to use than the C-style interface and plays well with STL algorithms that assume a functional programming style.

This functions in this adaptor library strictly uses BLAS operations, the data is not processed outside the BLAS calls.

## Contents
[[_TOC_]]

## Interfaces

In the _functional_ interface, most functions return special "views" rather than direct results.
The results are computed when converted to value types or assigned to other views.
Value types can be 2D (`multi::array<T, 2>`), 1D (`multi::array<T, 1>`) or 0D (`multi::array<T, 0>` or scalars `T`).
Views can be assigned to subarrays (e.g. `multi::subarray<T, 2>` or `multi::subarray<T, 1>`).

In this interface, functions like `gemv` generates views that can be assigned to values or from which constructors can be called in the library without unnecessary copies or allocations.
Expressions such as `multi::blas::gemv(alpha, A, x)` produce a range object that can be used in larger expressions, such as construction and assignemnt.

- Construction:
```cpp
multi::array<double, 1> const y = multi::blas::gemv(alpha, A, x);  // same effect as multi::array<double, 1> y({A.size()}); multi::blas::gemv(alpha, A, x, 0.0, y);
```
Among other advantages, the functional style gives the possibility of creating constant variables for results.

Like other parts of the library, when using `auto` and the unary `operator+` help generating concrete values.

```cpp
auto const y = +multi::blas::gemv(alpha, A, x);  // y variable is deduced as multi::array<double, 1> and latter constructed from the gemv operation
```

- Assignment:
```cpp
 multi::array<double, 1> y;  // empty vector
 y = multi::blas::gemv(alpha, A, x);  // same as  multi::blas::gemv(alpha, A, x, 0.0, y), y is resized if necessary
```

- Assignment (to subarray):
```cpp
multi::array<double, 2> Y;  // empty vector
Y[0] = multi::blas::gemv(alpha, A, x);  // same as  multi::blas::gemv(alpha, A, x, 0.0, Y[0]), Y[0] can't be resized because it is a subarray must have the correct size, 
```

- Compound-assign:
```cpp
multi::array<double, 1> y(A.size());
y += multi::blas::gemv(alpha, A, x);  // same as multi::blas::gemv(alpha, A, x, 1.0, y)
```

This interface plays well with the style of the STL algorithms.
For example, suppose we have a container of vectors, all of which need to be multiplied by a given array.

```cpp
std::list<multi::array<double, 1> > vs = ...;  // using std::list to avoid confusion
std::list<multi::array<double, 1> > ws = ...;
multi::array<double, 2> const A = ...;

std::transform(vs.begin(), vs.end(), ws.begin(), [&A](auto const& v) {return multi::blas::gemv(1.0, A, v);})
```

Although it shares some of the goals, this interface is independent of the [C++26 Linear Algebra Proposal](https://en.cppreference.com/w/cpp/numeric/linalg).
The main difference with other BLAS adaptors is that this library aims to offer a functional interface.

## Numeric Arrays, Conjugation Real and Imaginary parts

Just as with BLAS, the library supports element of real (`double` and `float`) and complex (`std::complex<double>` and `std::complex<float>`) types.
Other types that are semantically equivalent and binary-compatible (such as `thrust::complex`) also work directly.

## View manipulators

These functions produce views (not copies) related to conjugation, and transposition.
These typically replace the 'T', 'C' and 'N' characted arguments of the BLAS calls in the C or Fortran interfaces.

### `auto multi::blas::C(`_complex/real vector/matrix_`) -> `_complex/real vector/matrix view_

The conjugation operation is a unary operation that conjugates each element of the array, producing a view of the array that preserves the shape of the original array.

### `multi::blas::T(`_complex/real vector/matrix_`) -> `_complex/real vector/matrix view_

The transposition operation is a unary operation that transposes an array, producing a view of the array that transposed the elements (and the shape) of the original array.

### `multi::blas::N(`_complex/real vector/matrix_`) -> `_complex/real vector/matrix view_

This view returns the same array, implies no operations on the array; it is provided for completeness.

### `multi::blas::H(`_complex/real vector/matrix_`) -> `_complex/real vector/matrix view_

```cpp
using complex = std::complex<double>; 
complex const I{0.0, 1.0};

multi::array<complex, 2> B = {
    {1.0 - 3.0*I, 6.0 + 2.0*I},
    {8.0 + 2.0*I, 2.0 + 4.0*I},
    {2.0 - 1.0*I, 1.0 + 1.0*I}
};

namespace blas = multi::blas;
multi::array<complex, 2> conjB = blas::C(B);

assert( blas::C(B)[1][2] == std::conj(B[1][2]) );
assert( blas::T(B)[1][2] ==           B[2][1]  );
assert( blas::N(B)[1][2] ==           N[1][2]  );
assert( blas::H(B)[1][2] == std::conj(B[2][1]) );
```

Note that views do not play well with self-assignment.
```cpp
multi::array<double, 2> A({10, 10});
A = multi::blas::T(A);   // undefined behavior, this is not the right way to transpose a matrix in-place
```
The main purpose of these functions is to manipulate arguments to BLAS interface functions.

## BLAS level 1

(https://godbolt.org/z/Kjfa48d4P)

The functions in this level operate on one-dimensional arrays (vectors).
Here, we use `multi::array<T, 1>` as representative of a vector, but a one-dimensional subarray, such as a row or a column of a 2D array, can also be used as a vector.

### `auto multi::blas::copy(`_complex/real vector_`) -> `_convertible to complex/real vector_

Copies the values of a vector to another.

This is similar to assigment, except that it used the underlying BLAS function (including parallelization if offered by the BLAS implementation) and has marginal utility.
However, this case serves as illustration of the _functional_ interface, used in the rest of the library:
`multi::blas::copy(v)` doesn't copy or allocates anything, it creates a "view" that can serve different purposes, illustrated in 3 different cases:
1) The view can be used to construct a new vector (needing allocation),
Once again `operator+` helps with automatic type deduction.

```cpp
multi::array<double, 1> const v = {1.0, 2.0, 3.0};
multi::array<double, 1> const v2 = multi::blas::copy(v);  // case 1: allocates space for 3 elements and copies (using BLAS)
// auto const v2 = +v_copy;  // same effect as line above
```

(Note that `auto const v2 = v_copy;` would not create a value or perform a copy, it will simply hold a variable with the "copy range".
This is not recommended as it can be confused and create a dangling range.)

2) to assign to an existing vector (and resize it if is needed and possible)

```cpp
multi::array<double, 1> v3;  // mutable vector
v3 =  multi::blas::copy(v);  // case 2: resizes v3 (allocates space for 3 elements) and copies
```

```cpp
multi::array<double, 1> v4({3}, 0.0);  // allocates space for 3 elements
v4 =  multi::blas::copy(v);  // case 2: assigns copies (no allocation necessary)
```

3) to assign to a 1D subarray vector that is _not_ resizable.
The importance of this case is that it guarantees that no allocations are performed.

```cpp
 multi::array<double, 1> cA({3}, 0.0);  // allocates space for 3 elements
 v4({0, 2}) =  multi::blas::copy(v);  // case 3: LHS is not resizable, assigns copies to subrange (resizing is not possible or necessary, no allocations)
```

### `auto multi::blas::swap(`_complex/real vector_`, `_complex/real vector_`) -> void`

Swaps the values of two vectors.
Vector extensions must match.

Note that the utility of `multi::blas::copy` and `multi::blas::swap` is redundant with native features of the library (such as plain assignment, copy construction and swap), the only difference is that these operations will be performed using the BLAS operations elementwise one dimensional arrays (vectors) on real or complex data only.

These function do not work on 2D arrays (matrices) as the BLAS functions do not support this.
Copying or swapping 2D arrays with arbitrary layouts using BLAS could be done row-by-row: `std::transform(A2.begin(), A2.end(), B2.begin(), [](auto const& arow) {return multi::blas::copy(arow);})` (or, for efficiency column-by-column depending on the layout).

### `auto multi::blas::nrm2(`_complex/real vector_`) -> `_convertible to real scalar_

Unary operation that computes the norm-2 (Euclidean norm) of a vector.
The result is convertible to a real scalar

```cpp
multi::array<double, 1> const v = {1.0, 2.0, 3.0};
double const n = multi::blas::nrm2(v);
// auto const n = +multi::blas::nrm2(v);
```

```cpp
multi::array<std::complex<double>, 2> const v = { {1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0} };
double const n = multi::blas::nrm2(v[0]);  // acting on a row view
// auto const n = +multi::blas::nrm2(v[0]);
```

### `auto multi::blas::asum(`_complex/real vector_`) -> `_convertible to real scalar_

Returns the sum of the absolute values of the elements of a vector (norm-1).

### `auto multi::blas::iamax(`_complex/real vector_`) -> `_index_type_

Index of the element with the largest absolute value (zero-based)

### `auto multi::blas::dot(`_complex/real vector_, _complex/real vector_`) -> `_convertible to complex/real scalar_

Returns the dot product of two vectors with complex or real elements (`T`).

```cpp
multi::array<double, 1> const v = {1.0, 2.0, 3.0};
multi::array<double, 1> const w = {4.0, 5.0, 6.0};
double const d = multi::blas::dot(v, w); 
// auto const d = +multi::blas::dot(v, w);
```

Conjugation can be applied to either vector argument,

```cpp
    using multi::blas::dot;
    using multi::blas::C;

    auto const d1 = +dot(  v ,   w );
    auto const d2 = +dot(C(v),   w );
    auto const d3 = +dot(  v , C(w));
    auto const d4 = +dot(C(v), C(w));
```

It is important to note that the left hand side of the assignment can be a scalar that is part of a heap allocation.
In this case, the result is going to directly put at this location.

```cpp
multi::array<double, 1> z = {0.0, 0.0, 0.0};
z[1] = multi::blas::dot(v, w);
```

This feature regarding scalar results is essential when operating on GPU memory since the whole operation can be performed on the device.

> In CPUs BLAS `dot` functions has known bugs in different implementations.
> For example BLAS Apple Accelerate has bug in `sdot` while BLAS 32bit has a bug in `cdot`.
> In addition, some implementations of BLAS functions return the complex result in the stack and other write into a pointer.
> For this reason and for consistency, the library uses BLAS's `gemv` functions in place of `dot` in these cases.
> This should not affect the results.

### `auto multi::blas::scal(`_complex/real scalar`, `_complex/real vector_`)`

Scales a vector.

### `auto multi::blas::axpy(`_complex/real scalar`, `_complex/real vector_`) -> `_convertible to complex/real_

Vector addition.

```cpp
multi::array<double, 1> const x  = ...;
multi::array<double, 1> y = ...;
y += blas::axpy(2.0, x);  // same as blas:::axpy(+2.0, x, y)
y -= blas::axpy(2.0, x);  // same as blas:::axpy(-2.0, x, y)
```

## BLAS level 2

These functions operate on vectors and arrays.
Again, we use `multi::array<T, 1>` as representative of a vector, but a one-dimensional subarray, such as a row or a column of a 2D array, can also be used as a vector.
`multi::array<T, 2>` as representative of a matrices, but a two-dimensional subarray or larger of higher dimensional arrays can be used as long as one of the two interternal strides in 1.
This is limitation of BLAS, that only acts on certain layouts of 2D arrays.

### `auto multi::blas::gemv(`_complex/real scalar_ `,` _complex/real matrix_`) -> `_convertible to complex/real vector_

```cpp
multi::array<double, 2> const A({4, 3});
multi::array<double, 1> const x = {1.0, 2.0, 3.0};
multi::array<double, 1> const x = {1.0, 2.0, 3.0, 4.0};

y = blas::gemv(5.0, A, x);  // y <-  5.0 A * x
```

The gemv expression can be used for addition and subtraction,

```
y += blas::gemv(1.0, A, x);  // y <-  + A * x + y
y -= blas::gemv(1.0, A, x);  // y <-  - A * x + y
```

### GEMM

```cpp
#include<multi/array.hpp>
#include<multi/adaptors/blas.hpp>

namespace multi = boost::multi;

int main() {
 multi::array<double, 2> const A({2, 2});
 multi::array<double, 2> const B({2, 2});

 multi::array<double, 2> const C1 = multi::blas::gemm(1.0, A, B);
    auto const C2 = + multi::blas::gemm(1.0, A, B);
}
```
https://godbolt.org/z/d1E7donWM

(need linking to BLAS to work, e.g. `-lblas` or `-lopenblas` or `-lmkl`)

## Table of features

All these operations are now supported for CPU and GPU memory, real and complex.

scalars: `aa` ($`\alpha`$), `bb` ($`\beta`$) \
vectors: `x`, `y` \
matrices: `A`, `B`, `C`

vector operations: `C` (`*`) conjugation (element-wise) \
matrix operations: `J` (`*`) conjugation (element-wise) (use `C` for vectors), `T` transpose, `H` transpose conjugate (also `C`, discouraged), `U`/`L` upper or lower triangular part (logical zeroing other side)


| BLAS   | mutable form           | effect                        | operator form [³]        | functional form | thrust/STL [¹] |
|---     |---                     | ---                           | ---                  | ---             | --- |
| SWAP   |`blas::swap(x, y)`      | $x_i \leftrightarrow y_i$ | `(x^y)` |    | `swap_ranges(begin(x), end(x), begin(y))` |
| COPY   |`blas::copy(x, y)`      | $`y_i \leftarrow x_i`$ | `y << x` |  `y = blas::copy(x)` | `copy(begin(x), end(x), begin(y))` |
| ASUM   |`blas::asum(x, res)`    | $`r \leftarrow \sum_i \|\Re x_i\| + \|\Im x_i\|`$ | `x==0`/`x!=0` `isinf(x)` `is an(x)`[²] | `res = blas::asum(x)` | `transform_reduce(begin(x), end(x), 0.0, plus<>{}, [](auto const& e){return abs(e.real()) + abs(e.imag());})` |
| NRM2   |`blas::nrm2(x, res)`     | $`r \leftarrow \sqrt{\sum_i \|x_i\|^2}`$ | `abs(x)` | `res = blas::nrm2(x);` | `sqrt(trasnform_reduce(begin(x), end(x), 0.0, plus<>{}, [](auto const& e){return norm(e);}));` |
| SCAL   |`blas::scal(aa, x);`    | $`x_i \leftarrow \alpha x_i`$ | `x*=aa;`           |                 | `for_each(begin(x), end(x), [aa](auto& e){return e*=aa;})` |
| AXPY   |`blas::axpy(aa, x, y)`  | $`y_i \leftarrow \alpha x_i + y_i`$ | `y+=x` `y-=x` `y+=aa*x` `y-=aa*x` |     | `transform(x.begin(), x.end(), y.begin(), y.begin(), [aa](auto ex, auto ey) {return aa*ex + ey;}` |
| DOT    | `blas::dot(x, y, res)` | $`r = \sum_i x_i y_i`$        | `res = (x, y);`       | `res = blas::dot(x, y)`                | `inner_product(begin(x), end(x), begin(y), T{});` |
|        | `blas::dot(blas::C(x), y, res)` | $`r = \sum_i \bar x_i y_i`$ | `res = (*x, y);`  | `res = blas::dot(blas::C(x), y)`             | `inner_product(begin(x), end(x), begin(y), T{}, plus<>{}, [](T const& t1, T const& t2) {return conj(t1)*t2;});` |
|        | `blas::dot(x, blas::C(y), res)` | $`r = \sum_i x_i \bar y_i`$ | `res = (x, *y);`  | `res = blas::dot(x, blas::C(y));`             | `inner_product(x.begin(), x.end(), y.begin(), T{}, plus<>{}, [](T const& t1, T const& t2) {return t1*conj(t2);});` |
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
