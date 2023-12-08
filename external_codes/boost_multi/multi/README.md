<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.]Multi

> **Disclosure: This is not an official or accepted Boost library and is unrelated to the std::mdspan proposal.**

_© Alfredo A. Correa, 2018-2023_

_Multi_ is a modern C++ library that provides access and manipulation of data in multidimensional arrays, for both CPU and GPU memory.

Multidimensional array data structures are fundamental to several branches of computing, such as data analysis, image processing, and scientific simulations, and in combination with GPUs to Artificial Intelligence and Machine Learning.

This library offers array containers and views in arbitrary dimensions with well-behaved value semantics, offering total compatibility with the Standard Algorithms (STL), special memory (including GPUs), and following modern C++ design principles.
It requires at least C++17.

Some features of this library:

* Value semantics of multi-dimensional array containers
* Well-defined referential semantics of subarray (view) types
* Interoperability with other libraries, STL, ranges, thrust (CUDA and AMD GPUs), Boost, and C-libraries
* Fast access to elements and subarrays (views) types
* Arbitrary pointer types (fancy pointers, memory spaces)
* Simplified implementation (~4000 lines)

Do not confuse this library with [Boost.MultiArray](https://www.boost.org/doc/libs/1_69_0/libs/multi_array/doc/index.html), or with the standard MDSpan proposal `std::mdspan`.
`Multi` shares some of their goals but at a different level of generality.
The code is completely independent and with important differences in the implementation and semantics.

## Contents
[[_TOC_]]

## Using the library, installation and tests

You can try the library [online](https://godbolt.org/z/dvacqK8jE) before using it.

_Multi_ doesn't require installation, a single header `#include <multi/array.hpp>` is enough to use the full core library.
_Multi_ has no dependencies (except for the standard C++ library) and can be used immediately after downloading it.

```bash
git clone https://gitlab.com/correaa/boost-multi.git
```

Although installation is not necessary, the library can still be installed with CMake.
The header (and cmake) files will typically end up in `/usr/local/include/multi` and `/usr/local/share/multi`.

```bash
cd boost-multi
mkdir -p build && cd build
cmake ..  # --install-prefix=$HOME/.local
cmake --install .  # or sudo ...
```

_Testing_ the library requires Boost.Test library, installed for example via `sudo apt install cmake git g++ libboost-test-dev make` or `sudo dnf install boost-devel cmake gcc-c++ git`.
A CMake build system is provided to compile and run basic tests.

```bash
cmake --build .
ctest
```

Once installed, other CMake projects (targets) can depend on Multi by adding a simple `add_subdirectory(my_multi_path)` or by `find_package`:

```cmake
find_package(multi)  # see https://gitlab.com/correaa/boost-multi#using-the-library-installation-and-tests
```

Alternatively to `find_package` the library can be fetched on demand:
```cmake
include(FetchContent)
FetchContent_Declare(multi GIT_REPOSITORY https://gitlab.com/correaa/boost-multi.git)
FetchContent_MakeAvailable(multi)
...
target_link_libraries(my_target PUBLIC multi)
```

The code requires compilers with standard C++17 support, for reference any of:
LLVM's       `clang` [(5.0+)](https://godbolt.org/z/51E1hjfnn) (`libc++` and `libstdc++`), 
GNU's        `g++` [(7.1+)](https://godbolt.org/z/1nGEbKc5a), 
Nvidia's    [`nvcc`](https://godbolt.org/z/abdT73PqM) (11.4+) and `nvc++` (22.7+), 
Intel's      `icpc` (2021.2.0+) and `icpx` (2022.0.0+), 
Baxter's    [`circle`](https://www.circle-lang.org/) (build 187+),
and 
Microsoft's [MSVC](https://visualstudio.microsoft.com/vs/features/cplusplus/) (+19.14 in [conformant mode](https://godbolt.org/z/vrfh1fxWK)).

Optional "adaptor" sublibraries (included in `multi/adaptors/`) have specific dependencies, Boost.Serialization, fftw, blas, lapack, thurst, CUDA
(which can be installed with `sudo apt install libboost-serializa
tion-dev libfftw3-dev libblas64-dev liblapack64-dev libthrust-dev libcudart11.0` or `sudo dnf install blas-devel fftw-devel`.)
HIP support is experimental.

## Types

* `multi::array<T, D, A = std::allocator<T>>`: 
Array of integer positive dimension `D`, it has value semantics if element type `T` has value semantics. 
Memory is requested by an allocator of type `A` (supports stateful and polymorphic allocators).
* `multi::array_ref<T, D, P = T*>`: 
Array interpretation of a random access range, usually a contiguous memory block. 
It has reference semantics.
`P` does not need to be a language-pointer, it can be anything that behaves like a pointer (derreference and random access arithmetics).
`array<T, D, A>` is implicitly an `array_ref<T, D, A::pointer>`, the reverse conversion is only explicit.
* Other derived "unspecified types" fulfill a `MultiSubarray` concept, for example by taking partial indices or rotations (transpositions).
These derived referencial types can be named by life-time extensions `auto&&` or `auto const&`,
and they can decay to value types (copied elements).
`array`s automatically fulfill the concept being a `MultiSubarray`.
* `MultiSubarray::(const_)iterator`:
Iterator to subarrays of lower dimension.
For `D == 1` this is an iterator to an element (or scalar, with zero dimension).
This types are generated by `begin` and `end` (member) functions.
* `MultiSubarray::(const_)reference`:
Reference to subarrays of lower dimension.
For `D > 1`, these references are not true language-references, but types that emulate them (with reference semantics).
For `D == 1` this is a language reference to an element type (`T&`).
These types are generated by dereferencing iterators, (e.g. `*begin(MA)`), or by indexing (`MA[0]`).

## Basic Usage

The following code declares an array by specifying the element type and the dimension;
individual elements can be initialized from a nested rectangular list.
```cpp
multi::array<double, 2> A = {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0}
};

assert( A.size() == 2 );
assert( A.num_elements() == 6 );

assert( std::get<0>(A.sizes()) == 2 );
assert( std::get<1>(A.sizes()) == 3 );
```

The value of an array can be copied, (moved,) and compared;
copies are equal but independent.

```cpp
std::array<double, 2> B = A;
assert( extensions(B) == extensions(A) );
assert(  B       ==  A                 );
assert(  B[0][1] ==  A[0][1]           );
assert( &B[0][1] != &A[0][1]           );
```

Individual elements can be accessed by the multidimensional indices, either with square brackets (one index at a time) or with parenthesis (comma separated).

```cpp
assert( &A(1, 2) ==  &A[1][2] );
```

An arrays can be initialized from its sizes alone, in which case the element values are defaulted (or uninitialized):

```cpp
multi::array<double, 3> C({3, 4, 5});
assert( num_elements(C) == 3*4*5 );   // 60 elements with unspecified values
```

Arrays can be passed by value or by reference.
Most of the time, arguments should be passed through generic parameters to allow functions to also work with parts (subblocks, slices, etc.) of an array.
Usually, the most useful functions work on the _concept_ of an array rather than on a concrete type, for example:

```cpp
template<class ArrayDouble2D>  // instead of the overspecific argument std::array<double, 2>
auto element_1_1(ArrayDouble2D const& m) -> double const& {return m[1][1];}
...
assert( &element_1_1(A) == &A[1][1] );
```

The function expects any array or subarray of dimension 2 and return an element with type `double`. 

The generic function template arguments that are not intended to be modified are passed by `const&`; otherwise, they are passed by forward-reference `&&`.
In this way, the functions can be applied on subblocks of larger matrices.

```cpp
assert( &element_1_1(C3D[0]) == &C3D[0][1][1] );
```

## Advanced Usage

We can create a static C-array of `double`s, and refer to it via a bidimensional array `multi::array_ref<double, 2>`.

```cpp
#include "multi/array.hpp"

#include<algorithm>  // for sort
#include<iostream>  // for print

namespace multi = boost::multi;

int main() {
	double d_data[20] = {
		150.0, 16.0, 17.0, 18.0, 19.0,
		 30.0,  1.0,  2.0,  3.0,  4.0,
		100.0, 11.0, 12.0, 13.0, 14.0,
		 50.0,  6.0,  7.0,  8.0,  9.0
	};  // block of 20 elements ...
	multi::array_ref<double, 2> d2D_ref{&d_data[0], {4, 5}};  // interpreted as a 4 by 5 array
	...
```

Note that the syntax of creating a reference array involves passing the pointer to a memory block (20 elements here) and the logical dimensions of that memory block (4 by 5 here).

Next we print the elements in a way that corresponds to the logical arrangement:

```cpp
	...
	auto [is, js] = d2D_ref.extensions();
	for(auto i : is) {
		using std::cout;
		for(auto j : js) {
			cout<< d2D_ref[i][j] <<' ';
		}
		cout <<'\n';
	}
	...
```

This will output:

> ```
> 150 16 17 18 19  
> 30 1 2 3 4  
> 100 11 12 13 14  
> 50 6 7 8 9
> ```

The arrays provide iterator-based access, which allows it to interface with algorithms and implement new ones.

It is sometimes said (by Sean Parent) that the whole of STL algorithms can be seen as intermediate pieces to implement `std::stable_sort`. 
Pressumably, if one can sort over a range, one can perform any other standard algorithm.

```cpp
		...
		std::stable_sort( begin(d2D_ref), end(d2D_ref) );
		...
```

If we print the result, we will get:

> ```
> 30 1 2 3 4  
> 50 6 7 8 9  
> 100 11 12 13 14  
> 150 16 17 18 19
> ```

The array has been changed to be in row-based lexicographical order.
Since the sorted array is a reference to the original data, the original C-array has changed.

(Note that `std::sort` cannot be applied directly to a multidimensional C-array or to Boost.MultiArray types, among other libraries.
The arrays implemented by this library are, to the best of my knowledge, the only ones that support all STL algorithms directly.)

If we want to order the matrix in a per-column basis we need to "view" the matrix as range of columns.
This is done in the bidimensional case, by accessing the matrix as a range of columns:

```cpp
		...
		std::stable_sort( rotated(d2D_ref).begin(), rotated(d2D_ref).end() );
	}
```

Which will transform the matrix into:

> ```
> 1 2 3 4 30  
> 6 7 8 9 50  
> 11 12 13 14 100  
> 16 17 18 19 150 
> ```

In other words, by combining index rotations and transpositions, an array of dimension `D` can be viewed simultaneously as `D!` (D-factorial) different ranges of different "transpositions" (rotation/permutation of indices.)

## Initialization

`array_ref` is initialized from a preexisting contiguous range, the index extensions should compatible with the total number of elements.

```cpp
double* dp = new double[12];
multi::array_ref<double, 2> A({3, 4}, dp);
multi::array_ref<double, 2> B({2, 6}, dp);
...
delete[] dp;
```
Array references do not own memory and, just as language references, can not be rebinded (resized or "reseated") to refer to a different location.
Since `array_ref` is an array reference, it can "dangle" if the the original memory is deallocated.

Array objects (`multi::array`), in constrast, own the elements it contains and can be resized;
`array` is initialized by specifying the index extensions (and optionally a default value).

```cpp
multi::array<double, 1> A1({3}      , 11.0);  // {11.0, 11.0, 11.0}

multi::array<double, 2> A2({2, 3}   , 22.0);  // { {22.0, 22.0, 22.}, {22.0, 22.0, 22.0} }

multi::array<double, 3> A3({3, 2, 2}, 33.0);  // { { { 33., ...}, { ... }, ... } }
```
... or alternatively from a rectangular list.

```cpp
multi::array<double, 1> A1 = {1.0, 2.0, 3.0};
assert( num_elements(A1)==3 );

multi::array<double, 2> A2 {
	 { 1.0, 2.0, 3.0},
	 { 4.0, 5.0, 6.0}
};

assert( num_elements(A2) == 2*3);

multi::array<double, 3> const A3 = {
    {{ 1.2,  0.0}, { 2.4, 1.0}},
    {{11.2,  3.0}, {34.4, 4.0}},
    {{15.2, 99.0}, {32.4, 2.0}}
};

assert( A3.num_elements() == 3 * 2 * 2 );
```

In all cases, constness (`const` declaration) is honored in the expected way.

## Copy and assigment

The library offer value semantics for the `multi::array<T, D>` family of classes.
Constructing or assigning from an existing array generates a copy of the original object, that is, and object that is independent but equal in value.

```cpp
auto B2 = A2;  // same as multi::array<double, 2> B2 = A2;

assert(  B2       ==  A2       );  // copies have the same value (and also the same shape)
assert(  B2[0][0] ==  A2[0][0] )
assert( &B2[0][0] != &A2[0][0] );  // but they are independent
```

A (mutable) array can be assigned at any moment, independently of the previous state or shape (extensions).
The dimensionalities must match.
```cpp
B2 = A2;
```

(The operation can fail if there is no enough memory to hold a copy.)

Sometimes it is necessary to generate copies from views or subblocks.
```cpp
multi::array<double, 3> C2 = A2( {0, 2}, {0, 2} );
```
or equivalently,:
```cpp
auto C2 = + A2( {0, 2}, {0, 2} );
```
Note the use of the prefix `+` as an indicator that a copy must be created (it has no arithmetic implications).
Due to limitations of the language, omiting the `+` will create effectively another reference  non-indepdent view of the left-hand-side, which is generally undesired.

Subviews can also assigned but only if the shape of the left-hand side (LHS) and right-hand side (RHS) match.
Otherwise the behavior is undefined (in debug mode the program will fail an `assert`).

```cpp
C2( {0, 2}, {0, 2} ) = A2( {0, 2}, {0, 2} );  // both are 2x2 views of arrays, *elements* are copied
```

Introducing the same or overlapping arrays in the RHS and LHS produces undefined behavior in general (and the library doesn't check);
Notably, this instruction does not transpose the array, but produces an undefined result:

```cpp
A2 = A2.transposed();
```

While this instead does produce a transposition, at the cost of making a copy (`+`) of the tranposed array first and assigning (or moving) it back to the original array.

```cpp
A2 = + A2.transposed();
```

In-place transposition is an active subject of research; _optimal_ in speed and memory transpositions might require specially designed libraries.

Finally, arrays can be efficiently moved by transferring ownership of the internal data.

```cpp
auto B2 = std::move(A2);  // A2 is empty after this
```

Subarrays do not own the data therefore they cannot be moved in the same sense.
However, indivial elements of a view can be moved, this is particularly useful if the elements are expensive to copy.
A "moved" subview is simply another kind view of the elements.

```cpp
multi::array<std::vector<double>, 2> A({10, 10});
multi::array<std::vector<double>, 2> B({10, 10});
...
B[1] = A[2].element_moved();  // 10 *elements* of the third row of A is moved into the second row of B.
```

## Change sizes (extents)

Arrays can change their size while _preserving elements_ with the `reextents` method.

```cpp
multi::array<double, 2> A {
	 {1.0, 2.0, 3.0},
	 {4.0, 5.0, 6.0}
};

A.reextents({4, 4});

assert( A[0][0] == 1.0 );
```

Arrays can be emptied (zero-size) and memory is freed with `.clear()` (equivalent to `.reextents({0, ...})`).

The main purpose of `reextents` is element preservation.
Allocations are not amortized; 
except for trivial cases, all calls to reextend allocates and deallocates memory.
If element preservation is not desired, a simple assignment (move) from a new array expresses the intention better and it is more efficient since it doesn't need to copy preexisiting elements.

```cpp
A = multi::array<double, 2>({4, 4});  // like A.reextents({4, 4}) but elements are not preserved.
```

An alternative syntax, `.reextents({...}, value)` sets _new_ (not preexisting) elements to a specific value.

Subarrays or views cannot change their size or be emptied (e.g. `A[1].reextents({4})` or `A[1].clear()` will not compile).
For the same reason, subarrays cannot be assigned from an array or another subarray of a different size.

Changing the size of arrays by `reextents`, `clear`, or assignment generally invalidates existing iterators and ranges/views.
In contrast, `static_array<T, D>`, which can be used in many cases as a replacement of `array<T, D>`, doesn't have operations that invalidate iterators as it cannot be resized or assigned from arrays of different size.
(Preventing iterator invalidation can be a helpful technique in multithreaded programs.)

## Iteration

Accessing arrays by iterators (`begin`/`end`) enables the use of many iterator-based algorithms (see the sort example above).
`begin(A)/end(A)` (or equivalently `A.begin()/A.end()`) gives iterators that are linear and random access in the leading dimension.

Other non-leading dimensions can be obtained by "rotating" indices first.
`A.rotated().begin()/.end()` gives access to a range of subarrays in second dimension number (first dimension is put at the end).

`cbegin/cend` give constant (read-only access).
For example in a three dimensional array,

```cpp
	(cbegin(A)+1)->operator[](1).begin()[0] = 342.4;  // error, read-only
	( begin(A)+1)->operator[](1).begin()[0] = 342.4;  // assigns to A[1][1][0]
	assert( ( begin(A)+1)->operator[](1).begin()[0] == 342.4 );
```

As an example, this function allows printing arrays of arbitrary dimension into a linear comma-separated form.

```cpp
void flat_print(double const& d) { cout<<d; };  // terminating overload

template<class MultiArray>
void flat_print(MultiArray const& ma) {
	cout << "{";
	if(not ma.empty()) {
		flat_print(*cbegin(ma));  // first element
		std::for_each(cbegin(ma)+1, cend(ma), [](auto&& e) { cout<<", "; flat_print(e);});  // rest
	}
	cout << "}";
}
...
print(A);
```
> ```
> {{{1.2, 1.1}, {2.4, 1}}, {{11.2, 3}, {34.4, 4}}, {{15.2, 99}, {32.4, 2}}}
> ```

Except for those corresponding to the one-dimensional case, derreferencing iterators generally produce "proxy"-references (i.e. objects that behave in a large degree like language references).
These references can be given a name; using `auto` can be misleading since the resulting variable does not have value semantics.

```cpp
auto row = *begin(A);  // accepted by the language but misleading, row is not an independent value
```

In my experience, however, the following usage pattern produces a more consistent idiom for generating references (still without copying elements):

```cpp
auto&&       row0 = * begin(A);  // same as decltype(A)::      reference  row0 = * begin(A);
auto const& crow0 = *cbegin(A);  // same as decltype(A)::const_reference crow0 = *cbegin(A);

auto&&       row1 =               A [1];  // same as decltype(A)::      reference  row1 =               A [1];
auto const& crow1 = std::as_const(A)[1];  // same as decltype(A)::const_reference crow0 = std::as_const(A)[1];
```

If a new value is desired, these (equivalent) options express the intention more explicitly:

```cpp
decltype(A)::value_type row =   *begin(A);  // there is a real copy of the row
                   auto row = + *begin(A);  // there is another copy, note the use of '+' (unary plus)
```

### "Pointer" to subarray

The library strongly relies on value-sematics, and it doesn't entertain the concept of "shallow" copy; however, it supports refenece- and pointer-sematics.

Subarrays (e.g., rows in a 2D array) are reference-like objects with a concrete address-like value that identifies them uniquely.
These addresses, which behave like pointers, can be helpful to "mark" subviews; these markers can be copied and stored in arrays.

```cpp
auto A = multi::array<double, 2>({4, 4});

auto row2_ptr = &A[2];  // A[2] is a row of A (not an element)
assert( row2_ptr == &*(A.begin() + 2) );
```

The expression `A[2]` above is technically a C++ temporary object, and therefore it doesn't have a C++ address (taking `std::addressof` gives a compilation error). 
However, in the library's abstraction, `A[2]` references an existing part of the original array, i.e. it is a "library reference", whose "library address" can be obtained operator `&`. 
The case is an illustration that, in the library, operator `&` is, for subarrays, different than the `std::addressof` operator; the latter may not be defined and even not compile for some expressions.

Comparing these markers/pointers with different provenance, i.e., originating from different arrays, is generally undefined.

## Indexing

Arrays provide random access to elements or subviews.
Many algorithms on arrays are oriented to linear algebra,
which are ubiquitously implemented in terms of multidimensional index access.

Iterator access and index access are two alternatives for accessing elements.
For example `*(begin(A) + n)` and `A[n]` are equivalent
and the range defined by the pair `begin(A), end(A)` is equivalent to `A(extension(A))` and, in turn, to `A()` (even for a multidimensional array, `D > 1`).
The syntax can be combined in arbitrary ways, for example `*begin(A[n])` is equivalent to `A[n][0]`.

### Element access and partial access

Index access mimics that of C-fixed sizes arrays. 
For example, a 2-dimensional array will access to an element by specifying two indices `A[1][2]`,
which can be used for direct write and read operations; 
while _partial_ index arguments `A[1]` generate a view 1-dimensional object (a reference).

```cpp
A        // is a 2D value array
A[0]     // is a 1D "reference"/"view" array
A[0][0]  // is a an element reference, zero-D
```

Transpositions are also multidimensional arrays _views_ in which the index are *logically* rearranged, for example `rotated(m)[2][3][1] == m[1][2][3]`.
(`rotated`/`unrotated` refers to the fact that the logical _indices_ are rotated to the left/right.)

As an illustration of an algorithm based on index access (as opposed to iterators),
this example code implements Gauss Jordan Elimination without pivoting:

```cpp
template<class Matrix, class Vector>
auto gj_solve(Matrix&& A, Vector&& y) -> decltype(y[0]/=A[0][0], y) {
	std::ptrdiff_t Asize = size(A);
	for(std::ptrdiff_t r = 0; r != Asize; ++r) {
		auto&& Ar = A[r];
		auto&& Arr = Ar[r];
		for(std::ptrdiff_t c = r + 1; c != Asize; ++c) {Ar[c] /= Arr;}
		auto const yr = (y[r] /= Arr);
		for(std::ptrdiff_t r2 = r + 1; r2 != Asize; ++r2) {
			auto&& Ar2 = A[r2];
			auto const& Ar2r = Ar2[r];  // auto&& Ar = A[r];
			for(std::ptrdiff_t c = r + 1; c != Asize; ++c) {Ar2[c] -= Ar2r*Ar[c];}
			y[r2] -= Ar2r*yr;
		}
	}
	for(std::ptrdiff_t r = Asize - 1; r > 0; --r) {
		auto const& yr = y[r];
		for(std::ptrdiff_t r2 = r-1; r2 >=0; --r2) {y[r2] -= yr*A[r2][r];}
	}
	return y;
}
```

This function can be applied to a `multi::array` container:

```cpp
multi::array<double, 2> A = {{-3.0, 2.0, -4.0},{0.0, 1.0, 2.0},{2.0, 4.0, 5.0}};
multi::array<double, 1> y = {12.0, 5.0, 2.0};  // (M); assert(y.size() == M); iota(y.begin(), y.end(), 3.1);
gj_solve(A, y);
```

and also to a combination of `MultiArrayView`-type objects (including standard vectors):

```cpp
multi::array<double, 2> A({6000, 7000}); std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 0.1);
std::vector<double> y(3000); std::iota(y.begin(), y.end(), 0.2);  // could be also a multi::array<double, 1> y({3000});
gj_solve(A({1000, 4000}, {0, 3000}), y);
```

### Slices and strides

Given an array, a slice in the first dimension can be taken with the `sliced` function. 
`sliced` takes two arguments, the first index of the slice and the last index (not included) of the slice. For example,

```cpp
multi::array<double, 2> A({4, 5});  // A is a value
assert( std::get<0>(A.sizes()) == 4 );
assert( std::get<1>(A.sizes()) == 5 );

auto&& A_sliced = A.sliced(1, 3); // {{d2D[1], d2D[2]}}
assert( std::get<0>(A_sliced.sizes()) == 2 );
assert( std::get<1>(A_sliced.sizes()) == 5 );
```

The number of rows in the sliced matrix is 2 because we took only two rows, row 1 and row 2 (row 3 is excluded).

In the same way a strided view of the original array can be taken with the `strided` function.

```cpp
auto&& d2D_strided = d2D.strided(2); // {{ d2D[0], d2D[1] }};
assert( d2D_strided.size(0) == 2 and d2D_strided.size(1) == 5 );
```

In this case the number of rows is 2 because, out of the 4 original rows we took one every two.

Operations can be combined in a single line:

```cpp
auto&& d2D_slicedstrided = d2D.sliced(1, 3).strided(2); // {{ d2D[1] }};
assert( std::get<0>(d2D_slicedstrided.sizes()) == 1 and std::get<1>(d2D_slicedstrided.sizes()) == 5 );
```

For convenience, `A.sliced(a, b, c)` is the same as `A.sliced(a, b).strided(c)`.

By combining `rotated`, `sliced` and `strided` one can take sub arrays at any dimension index.
For example in a two dimensional array one can take a subset of columns by defining.

```cpp
auto&& subA = A.rotated().sliced(1, 3).strided(2).unrotated();
```

Other notations are available, for example this is equivalent to `A(multi::all, {1, 3, /*every*/2})` or `~(~A)({1, 3, 2})`.
The `rotated/strided/sliced/rotated` and combinations of them provides the most control over the subview operations.

Blocks (slices) in multidimensions can be obtained by pure index notation using parentheses `()` (`.operator()`):

```cpp
auto        A = multi::array<double, 2>({6, 7});  // 2D value array

auto&&      A_block1 = A({1, 4}, {2, 4});  // 2D subarray reference (modifiable)
auto const& A_block2 = A({1, 4}, {2, 4});  // 2D subarray reference (non-modifiable)

auto        A_block3 = A({1, 4}, {2, 4});  // works but it can be confusing, use `auto&&` instead
```

Sometimes copies are necessary, specifically from a subarray block, this can be done by constructing a new array. 
The value array can be deduced by using `auto` and the `decay` member, which in turn is equivalent to the prefix `+` operator.

```cpp
multi::array<double, 2> block_value_1 =   A({1, 4}, {2, 4})        ;
auto                    block_value_2 =   A({1, 4}, {2, 4}).decay();
auto                    block_value_3 = + A({1, 4}, {2, 4})        ;
```

Any parenthesis argument can be either a range (with or without stride) or an index. 
Range argument can be substituted by `multi::all` to obtain the whole range.

## Conversions

Conversion between arrays of distinct types is possible if the underlying elements allow it.
The result is as if elements are converted one by one; array sizes (extensions) are preserved.
Allowed conversions can be implicit or explicit and reflect the behavior of the element types.

```cpp
// implicit conversions from real to complex is allowed ...
double                  d = 5.0;     std::complex<double>                  z = d;
// ... therefore it is also allowed from array of reals to arrays of complex
multi::array<double, 2> D({10, 10}); multi::array<std::complex<double>, 2> Z = D;
// (implicit or explicit) conversions from real to complex are disallowed (compilation error)
// multi::array<double, 2> D = Z;  // or D{Z};
```

Another case is illustrated by `std::complex<float>` and `std::complex<double>`; 
in one direction, the conversion can be implicit, while in the other, it can only be explicit.
This behavior is reflected in the corresponding arrays:
```cpp
multi::array<std::complex<float>>  C;
multi::array<std::complex<double>> Z = C;  // implicit conversion ok
multi::array<std::complex<float>>  C2{Z};  // explicit conversion is allowed
// multi::array<std::complex<float>>  C3 = Z;  // implicit conversion is disallowed (compilation error)
```

Implicit conversions are generally considered harmful, but inconsistent conversions are worst; therefore, the library allows them when appropriate.
The main drawback of implicit conversions in this context is that they might incur unexpected (e.g. costly) data conversions when passing arguments to functions.

```cpp
void fun(multi::array<std::complex<double>> Z) { ... };
...
multi::array<double, 2> D({10, 10});
fun(D);  // real elements are converted to complex silently here
```
In many instances, specially in generic code, it might still be a desirable behavoir.

To prevent implicit conversions, use element types with no implicit conversions when possible.

Finally, arrays of unrelated element types are prevented from producing direct conversions, resulting in compilation errors.
Element-wise transformations can be used instead.
For example, to convert an array of integers to an array of text strings:

```cpp
	multi::array<int, 2> const A = {{1, 2}, {3, 4}};

	auto to_string = [](int e) {return std::to_string(e);};
	multi::array<std::string, 2> B = A.element_transformed(to_string);
	assert( B[1][1] == "4" );
```

## Const-correctness

Const-correctness refers to the property of a program to disallow mutation of certain objects when it is undesired or logically incorrect.
Honoring the const-ness declaration is fundamental not only to avoid bugs and typos but also for thread safety and generic programming.
The library goes to great lengths to ensure const-correctness for the whole or parts of any object.

Arrays are resizable, and their elements can be mutated unless declared constant (using the keyword `const`).

A reference array or subarray is never resizable, but its elements are mutable if not declared `const`.
The design ensures that the const-ness of references and values propagates to subarrays (views) and, ultimately, their elements.

```cpp
template<class Array1D>
void print(Array1D const& coll) {
//  *coll.begin() = 99;  // doesn't compile, "assignment of read-only location"

	for(auto const& e : coll) {std::cout<< e <<", ";}
	std::cout << std::endl;
}

int main() {
	multi::array<int, 1> const coll1 = {0, 8, 15, 47, 11, 42};

	print( coll1 );  // prints "0, 8, 15, 47, 11, 42"
	print( coll1({0, 3}) );  // prints "0, 8, 15"
}
```

As a general rule for passing generic arrays as arguments, pass them as `Array const&` (in the context of `template<class Array>`);
unless mutation is expected, in which case take arguments as `Array&&` (note the double ampersand, i.e., universal/forwarding reference).
Analogously, subarrays can be locally *named* into "constant language references" using `auto const&` and, if mutation is desired, `auto&&` should be used.
Regular references `Array&` or `auto&` in general do not have the expected behavior for views.

```cpp
template<class Array1D>
void fill_99(Array1D&& coll) {
	for(auto& e : coll) { e = 99; }
}

int main() {
	multi::array<int, 1> coll1 = {0, 8, 15, 47, 11, 42};

	fill_99( coll1 );
	fill_99( coll1({0, 3}) );

	auto&& coll1_take3 = coll1({0, 3});
	fill_99( coll1_take3 );

	auto const& coll2 = coll1;
//  fill_99( coll2 );  // doesn't compile because coll2 is const
//  fill_99( coll2({0, 3}) );  // similar to coll2 | take(3) doesn't compile

	auto const& coll1_take3_const = coll1({0, 3});
//  fill_99( coll1_take3_const );  // doesn't compile because coll1_take3_const is const
}
```

## Compile-time evaluation (constexpr-all-the-things)

With certain limitations imposed by the language, arrays can be declared in contexts with compile-time evaluation.

```cpp
constexpr auto trace() {
	multi::array<int, 2> arr = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	arr[2][2] = 10;
	return std::accumulate(arr.diagonal().begin(), arr.diagonal().end());
}

static_assert( trace() == 4 + 2 + 10 );
```
https://godbolt.org/z/Porre3z8s

## Broadcast (infinite views)

Broadcasting is a technique by which arrays are reinterpreted as having a higher dimension by repeating elements.
The technique allows the reuse of operations designed for high dimensionality and effectively apply them to arrays of lower dimensionality.
The result is generally an economy in the number of distinct operations that need to be provided in exchange for understanding how and where to exploit the broadcast operations.

Broadcasting is popular in array-based languages, such as Julia and NumPy, and the broadcast is generally applied automatically to match the dimension expected by the operation and other operation inputs.
The library provides a basic form of broadcasting with certain limitations.

Here is an example of an algorithm designed for two 2D arrays to obtain the row-by-row inner product.

```cpp
    auto row_by_row_dot = [](auto const& A2D, auto const& B2D, auto& results) {
        std::transform(A2D.begin(), A2D.end(), B2D.begin(), results.begin(),
            [](auto const& Arow, auto const& Brow) {return std::inner_product(Arow.begin(), Arow.end(), Brow.begin(), 0);}
        );
    };

    auto A = multi::array<int, 2>{{ 0,  1}, { 2,  3}, { 4,  5}};
    auto B = multi::array<int, 2>{{10, 11}, {12, 13}, {14, 15}};

    auto dots = multi::array<int, 1>({A.size()});

    row_by_row_dot(A, B, dots);
```

If, for some reason, we want to obtain the inner product against a _single_ right-hand vector instead of several (a single 1D array of two elements), we would need to (re)write the function (or copy the repeated vector into the 2D `B` array, which is not ideal.)
Broadcasting can help reuse the same function without changes.

```cpp
    multi::array<int, 1> b = {10, 11};

    row_by_row_dot(A, b.broadcasted(), dots);
```

The alternative, not using broadcast, is to write a very similar function,

```cpp
    auto row_fixed_dot = [](auto const& A2D, auto const& b1D, auto& results) {
        std::transform(A2D.begin(), A2D.end(), results.begin(),
            [&b1D](auto const& Arow) {return std::inner_product(Arow.begin(), Arow.end(), b1D.begin(), 0);}
        );
    };

    row_fixed_dot(A, b, dots3);
```
(https://godbolt.org/z/9ndvfKqhc)

Broadcasted arrays do not behave like normal array views in several aspects:
First, broadcasted arrays are infinite in the broadcasted dimension; iteration will never reach the end position, and calling `.size()` is undefined behavior.
Explicit loops or algorithms that depend on reaching `.end()` from `.begin()` will effectively be non-terminating.
Second, these array views are strictly read-only and alias their element addresses, e.g. `&b.broadcasted()[1][0] == &b.broadcasted()[2][0]` (since internal layouts' strides can be zero).

For illustration purposes only, `fill` here is replaced by `copy`; problematic uses are highlighted:

```cpp
    multi::array<double, 2> B({10, 2});
    std::fill(B.begin(), B.end(), b);                                       // canonical way
    std::copy_n(b.broadcasted().begin(), v.size(), v.end());                // equivalent, using broadcast

    std::copy_n(b.broadcasted().begin(), b.broadcasted().size(), v.end());  // incorrect, undefined behavior, no useful size()
    std::copy  (b.begin(), b.end(), v.begin());                             // incorrect, undefined behavior, non-terminating loop
	B = b.broadcasted();                                                    // incorrect, undefined behavior
```

Unlike popular languages, broadcasting is not automatic in the library and is applied to the leading dimension only, one dimension at a time.
Broadcasting in non-leading dimensions can be achieved by transpositions and index rotation.

Abuse of broadcast can make it harder to reason about operations; its primary use is to reuse existing efficient implementations of algorithms when implementations for a specific lower dimensions are not available.
These algorithms need to be compatible with broadcasted views (e.g., no explicit use of `.size()` or infinite loops stemming from problematic use of `.begin()/end()`.)

As a final example, consider a function that computes the elements-by-element product of two 2D arrays,

```cpp
    auto hadamard = [](auto const& A, auto const& B, auto&& C) {
        auto const [is, js] = C.extensions();
        for(auto i : is) for(auto j : js) C[i][j] = A[i][j]*B[i][j];
    };
```

As it is, this function can be reused to calculate the outer product of two 1D arrays:

```cpp
    auto outer = [&]<typename T>(auto const& a, auto const& b, T&& C) {
        return hadamard(~(a.broadcasted()), b.broadcasted(), std::forward<T>(C));
    };
```
(https://godbolt.org/z/5o95qGdKz)

Note that the function acting on 2D arrays, doesn't use the undefined (infinite) sizes (second dimension of `A` and first dimension of `B`).

## Partially formed elements

The library can take advantage of types with [partially formed](https://marcmutz.wordpress.com/tag/partially-formed-state/) state when
elements are trivial to construct (e.g., built-in types).
In such cases, `multi::array` does not initialize individual elements unless specified.
If trivial construction is unavailable, the library uses the default constructor.

For example, after construction, the values of the six elements of this array are unspecified (partially formed).
```cpp
multi::array<int, 2> A2({2, 3});
```

No behavior of the program should depend on these values. 
(Address sanitizers and memory checkers can detect this.)
This design is a slight departure from the STL's design, which [immediatelly initializes elements in containers](https://lemire.me/blog/2012/06/20/do-not-waste-time-with-stl-vectors/).

For types that afford partially formed states, elements can be later specified via assignment or assigning algorithms (e.g., copy or transform destination).
Initialization can be enforced by passing a value argument after the extensions.
```cpp
multi::array<int, 2> A2({2, 3}, 0);  // generically multi::array<T, 2>({2, 3}, T{}); or multi::array<T, 2>({2, 3}, {})
```

This design is particularly advantageous for *numeric* types for which external low-level libraries can fill values.
(or when data sits in GPUs, where the initialization step would require an expensive kernel launch and subsequent synchronization).

Unfortunately, regarding the numeric types, STL's `std::complex<double>` was standardized as not-trivially constructible.
A workaround is possible by forcing a particular flag on the client code in global scope, for example, immediately after including the library:
```cpp
#include<multi/array.hpp>
...
template<> inline constexpr
bool multi::force_element_trivial_default_construction<std::complex<double>> = true;
```

With this line, `std::complex<double>` elements inside arrays will be left uninitialized unless a value is specified.
The rule will only apply to this library's containers (`multi::array`, etc), and not to other containers (such as `std::vector`) or individual `std::complex` variables.

## Type Requirements

Thelibrary design tries to impose the minimum possible requirements over the types that parameterize the arrays.
Array operations assume that the contained type (element type) are regular (i.e. different element represent disjoint entities that behave like values).
Pointer-like random access types can be used as substitutes of built-in pointers.
(Therefore pointers to special memory and fancy-pointers are supported.)

### Linear Sequences: Pointers

An `array_ref` can reference an arbitrary random access linear sequence (e.g. memory block defined by pointer and size).
This way, any linear sequence (e.g. `raw memory`, `std::vector`, `std::queue`) can be efficiently arranged as a multidimensional array.

```cpp
std::vector<double> buffer(100);
multi::array_ref<double, 2> A({10, 10}, buffer.data());
A[1][1] = 9.0;

assert( buffer[11] == 9.0 );  // the target memory is affected
```
Since `array_ref` does not manage the memory associated with it, the reference can be simply dangle if the `buffer` memory is reallocated (e.g. by vector-`resize` in this case).

### Special Memory: Pointers and Views

`array`s managate their memory behind the scenes through allocators, which can be specified at construction.
It can handle special memory, as long as the underlying types behave coherently, these include [fancy pointers](https://en.cppreference.com/w/cpp/named_req/Allocator#Fancy_pointers) (and fancy references).
Associated fancy pointers and fancy reference (if any) are deduced from the allocator types.
Another use of fancy pointer is to create by-element "projections".

#### Allocators and Fancy Pointers

Specific uses of fancy memory are file-mapped memory or interprocess shared memory.
This example illustrates memory persistency by combining with Boost.Interprocess library. 
The arrays support their allocators and fancy pointers (`boost::interprocess::offset_ptr`).

```cpp
#include <boost/interprocess/managed_mapped_file.hpp>
using namespace boost::interprocess;
using manager = managed_mapped_file;
template<class T> using mallocator = allocator<T, manager::segment_manager>;
decltype(auto) get_allocator(manager& m) {return m.get_segment_manager();}

template<class T, auto D> using marray = multi::array<T, D, mallocator<T>>;

int main() {
{
	manager m{create_only, "mapped_file.bin", 1 << 25};
	auto&& arr2d = *m.construct<marray<double, 2>>("arr2d")(marray<double, 2>::extensions_type{1000, 1000}, 0.0, get_allocator(m));
	arr2d[4][5] = 45.001;
}
// imagine execution restarts here, the file "mapped_file.bin" persists
{
	manager m{open_only, "mapped_file.bin"};
	auto&& arr2d = *m.find<marray<double, 2>>("arr2d").first;
	assert( arr2d[7][8] == 0. );
	assert( arr2d[4][5] == 45.001 );
	m.destroy<marray<double, 2>>("arr2d");
}
}
```

(See also, examples of interactions with the CUDA Thrust library to see more uses of special pointer types to handle special memory.)

#### Transformed views

Another kind of fancy-pointer is one that transforms the underlying values.
These are useful to create "projections" or "views" of data elements.
In the following example a "transforming pointer" is used to create a conjugated view of the elements.
In combination with transposed view, it can create a hermitic (transposed-conjugate) view of the matrix (without copying elements).
We can adapt the library type `boost::transform_iterator` to save coding, but other libraries can be used also.
The hermitized view is read-only, but with additional work a read-write view can be created (see `multi::blas::hermitized` in multi-adaptors).

```cpp
constexpr auto conj = [](auto const& c) -> auto const {return std::conj(c);};

template<class T> struct conjr : boost::transform_iterator<decltype(conj), T*> {
	template<class... As> conjr(As const&... as) : boost::transform_iterator<decltype(conj), T*>{as...} {}
};

template<class Array2D, class Complex = typename Array2D::element_type>
auto hermitized(Array2D const& arr) {
	return arr
		.transposed() // lazily tranposes the array
		.template static_array_cast<Complex, conjr<Complex>>(conj)  // lazy conjugate elements
	;
}

int main() {
    using namespace std::complex_literals;
	multi::array A = {
		{ 1.0 + 2.0i,  3.0 +  4.0i},
		{ 8.0 + 9.0i, 10.0 + 11.0i}
	};

	auto const& Ah = hermitized(A);

	assert( Ah[1][0] == std::conj(A[0][1]) );
}
```

To simplify this bolier plate, the library provides the `.element_transformed(F)` method that will apply a transformation `F` to each element of the array.
In this example the original arrays is transformed into a transposed array with duplicated elements.

```cpp
    multi::array<double, 2> A = {
        {1.0, 2.0},
        {3.0, 4.0},
    };

    auto const scale = [](auto x) { return x * 2.0; };

    auto B = + A.rotated().element_transformed(scale);
	assert( B[1][0] == A[0][1] * 2 );
```
([live](https://godbolt.org/z/b7E56Mjc8))

Since `elements_transformed` is view (reference) to the original data, it is important to understand the semantics of evaluation an possible allocations associated with it.
As mentioned in other sections using `auto` and/or `+` appropriately can lead to simple and efficient expressions.

| Construction    | Allocation of `T`s | Initialization (of `T`s) | Evaluation (of `fun`) | Notes |
| -------- | ------- | ------- | ------- | ------- |
| `multi::array<T, D> [const] B = A.element_transformed(fun);` | Yes        | No  | Yes | Implicit conversion to `T` if result is different, dimensions must match   |
| `multi::array<T, D> [const] B = + A.element_transformed(fun);` | Yes (and move, or might allocate twice if types don't match)  | No  | Yes | Not recommended | 
| `multi::array<T, D> [const] B{A.element_transformed(fun)};` | Yes        | No  | Yes | Explicit conversion to `T` if result is different, dimensions must match   |
| `auto [const] B = + A.elements_transformed(fun);`           | Yes         | No  | Yes | Types and dimension are deduced, result is contiguous, preferred |
| `auto [const] B = A.element_transformed(fun);`               | Yes         | No  | No (delayed) | Result is effective a reference, may dangle with `A`, usually `const`, not recommended   |
| `auto[&&\|const&] B = A.elements_transformed(fun);`           | Yes         | No  | No (delayed) | Result is effective a reference, may dangle with `A`, usually `const&`, preferred way  |

| Assigment    | Allocation of `T`s | Initialization (of `T`s) | Evaluation (of `fun`) | Notes |
| -------- | ------- | ------- | ------- | ------- |
| `B = A.elements_transformed(fun);`           | No, if sizes match | Possibly (when `B` was initialized)  | Yes | `B` can't be declared `const`, it can be a writable subarray, preferred  |
| `B = + A.elements_transformed(fun);`           | Yes | Possibly (when `B` was initialized)  | Yes | Not recommended |

# Interoperability with other software

## STL (Standard Template Library)

The fundamental goal of the library is that the arrays and iterators can be used with STL algorithms out-of-the-box with a reasonable efficiency.
The most dramatic example of this is that `std::sort` works with array as it is shown in a previous example.

Along with STL itself, the library tries to interact with other existing quality C++ libraries listed below.

### Ranges (C++20)

[Standard ranges](https://en.cppreference.com/w/cpp/ranges), together with the array constructors provided by the library, enable a functional programming style;
this allows to work with immutable variables in many cases in place of mutable imperative code.

```cpp
    multi::array<double, 2> const A = {{...}};
    multi::array<double, 1> const V = {...};

    multi::array<double, 1> const R = std::views::zip_transform(std::plus<>{}, A[0], V);
	// multi::array<double, 1> R(V.size());  // in the alternative imperative code, R is created...
    // for(auto i : R.extension()) {R[i] = A[0][i] + V[i];}  // ...then mutated
```
[(live)](https://godbolt.org/z/M84arKMnT)

The library also works well with Ranges-v3 which is approximately a superset of STL ranges (see below).

### Polymorphic Memory Resources

In addition to supporting classic allocators (`std::allocator` by default), the library is compatible with C++17's [polymorphic memory resources (PMR)](https://en.cppreference.com/w/cpp/header/memory_resource) which allows using advanced allocation strategies, including preallocated buffers.
For example, this code uses a buffer as memory for two arrays; this buffer ends up containing the data of the arrays `"aaaabbbbbbXX"`.

```cpp
#include <memory_resource>  // for polymorphic memory resource, monotonic buffer

int main() {
	char buffer[13] = "XXXXXXXXXXXX";  // a small buffer on the stack
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::pmr::array<char, 2> A({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> B({3, 2}, 'b', &pool);

	assert( buffer == std::string{"aaaabbbbbbXX"} );
}
```
`multi::pmr::array<T, D>` is a synonym for `multi::array<T, D, std::pmr::polymorphic_allocator<T>>`.
In this particular example, the technique can be used to avoid dynamic memory allocations of small local arrays. [(live)](https://godbolt.org/z/fP9P5Ksvb)

The library also supports allocators from other libraries, including those returning special pointer types (see [CUDA Thrust](#cuda-thrust) Thurst section).

### Substitutability with standard vector and span

The one-dimensional case `multi::array<T, 1>` is special and overlaps functionality with other dynamic array implementations, such as `std::vector`.
Indeed, both types of containers are similar and usually substitutable, with no or minor modifications.
For example, both can be constructed from a list of elements (`C c = {x0, x2, ...};`) or from a size `C c(size);`, where `C` is either type.

Both values are assignable, have the same element access patterns and iterator interface, and implement all (lexical) comparisons.

They differ conceptually in their resizing operations: `multi::array<T, 1>` doesn't insert or push elements and resizing works differently.
The difference is that the library doesn't implement *amortized* allocations; therefore, these operations would be of a higher complexity cost than the `std::vector`.
For this reason, `resize(new_size)` is replaced with `reextent({new_size})` in `multi::array`, whose primary utility is for element preservation when necessary.

With the appropriate specification of the memory allocator, `multi::array<T, 1, Alloc>` can refer to special memory not supported by `std::vector`.

Finally, an array `A1D` can be copied by `std::vector<T> v(A1D.begin(), A1D.end());` or `v.assign(A1D.begin(), A1D.end());` or viseversa.
Without copying, a reference to the underlying memory can be created `auto&& R1D = multi::array_ref<double, 1>(v.data(), v.size());` or conversely `std::span<T>(A1D.data_elements(), A1D.num_elements());`. 
(See examples [here](https://godbolt.org/z/n4TY998o4).)

The `std::span` (C++20) has not a well defined reference- or pointer-semantics; it doesn't respect `const` correctness in generic code.
This behavior is contrary to the goals of this library;
and for this reason, there is no single substitute for `std::span` for all cases.
Depending on how it is used, either `multi::array_ref<T, 1> [const& | &&]` or `multi::array_ptr<T [const], 1>` may replace the features of `std::span`.
The former typically works when using it as function argument.

## Comparison to other array libraries (mdspan, Eigen, etc)

The C++23 standard is projected to provide `std::mdspan`, a non-owning _multidimensional_ array.
So here is an appropriate point to compare the two libraries.
Although the goals are similar, the two libraries differ in their generality and approach; in a few words: 

The Multi library concentrates on _well-defined value- and reference-semantics of arbitrary memory types with regularly arranged elements_ (distributions described by strides and offsets) and _extreme compatibility with STL algorithms_ (via iterators) and other fundamental libraries.

`mdspan` concentrates on _arbitrary layouts_ for non-owning memory of a single type (CPU raw pointers).
Due to the priority of arbitrary layouts, the `mdspan` research team didn't find efficient ways to introduce iterators into the library. 
Therefore, its compatibility with the rest of the STL is lacking.
The ultimate reason is that arbitrary layouts do not compose well across subdimensions, and, in turn, this imposes certain limitations in `mdspan`, such as ad-hoc slicing and subarray.

[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) is a very popular matrix linear algebra library, and as such, it only handles the special 2D (and 1D) array case.
Instead, the Multi library is dimension-generic and doesn't make any algebraic assumptions for arrays or contained elements (but still can be used to _implement_ dense linear algebra algorithms.)

Here is a table comparing with `mdspan`, R. Garcia's [Boost.MultiArray](https://www.boost.org/doc/libs/1_82_0/libs/multi_array/doc/user.html) and Eigen. 
[(online)](https://godbolt.org/z/555893MqW).


|                             | Multi                                                           | mdspan                                                                          | Boost.MultiArray (R. Garcia)                                                                         | Inria's Eigen                                                                           |
|---                          | ---                                                             | ---                                                                             | ---                                                                                                  | ---                                                                                     |
| No external Deps            | **yes** (only Standard Library C++17)                           | **yes** (only Standard Library)                                                 | **yes** (only Boost)                                                                                 | **yes**                                                                                 |
| Arbritary number of dims    | **yes**, via positive dimension (compile-time) parameter `D`    | **yes**                                                                         | **yes**                                                                                              | no  (only 1D and 2D)                                                                    |
| Non-owning view of data     | **yes**, via `multi::array_ref<T, D>(ptr, {n1, n2, ..., nD})`   | **yes**, via `mdspan m{T*, extents{n1, n2, ..., nD}};`                          | **yes**, via `boost::multi_array_ref<T, D>(T*, boost::extents[n1][n2]...[nD])` | **yes**, via `Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>(ptr, n1, n2)` |
| Compile-time dim size       | no                                                              | **yes**, via template paramaters `mdspan{T*, extent<16, dynamic_extents>{32} }` | no                                                                             | **yes**, via `Eigen::Array<T, N1, N2>` |
| Array values (owning data)  | **yes**, via `multi::array<T, D>({n1, n2, ..., nD})`            | no, (planned `mdarray`)                                   | **yes**, via `boost::multi_array<T, D>(boost::extents[n1][n2]...[nD])` | **yes**, via `Eigen::Array<T>(n1, n2)` |
| Value semantic (Regular)    | **yes**, via cctor, mctor, assign, massign, auto decay of views | no, and not planned                                       | partial, assigment on equal extensions  | **yes** (?) |
| Move semantic               | **yes**, via mctor and massign                                  | no                                                        | no (C++98 library)                      | **yes** (?) |
| const-propagation semantics | **yes**, via `const` or `const&`                                | no, const mdspan elements are assignable!                 | no, inconsistent                        | (?) |
| Element initialization      | **yes**, via nested init-list                                   | no                                                        | no                                      | no, only delayed init via `A << v1, v2, ...;` |
| References w/no-rebinding   | **yes**, assignment is deep                                     | no, assignment of mdspan rebinds!                         | **yes**                                 | **yes** (?) |
| Element access              | **yes**, via `A(i, j, ...)` or `A[i][j]...`                     | **yes**, via `A(i, j, ...)`                               | **yes**, via `A[i][j]...`               | **yes**, via `A(i, j)` (2D only) |
| Partial element access      | **yes**, via `A[i]` or `A(i, multi::all)`                       | **yes**, via `submdspan(A, i, full_extent)`               | **yes**, via `A[i]`                     | **yes**, via `A.row(i)` |
| Subarray views              | **yes**, via `A({0, 2}, {1, 3})` or `A(1, {1, 3})`              | **yes**, via `submdspan(A, std::tuple{0, 2}, std::tuple{1, 3})` | **yes**, via `A[indices[range(0, 2)][range(1, 3)]]` | **yes**, via `A.block(i, j, di, dj)` |
| Subarray with lower dim     | **yes**, via `A(1, {1, 3})`                                     | **yes**, via `submdspan(A, 1, std::tuple{1, 3})`          | **yes**, via `A[1][indices[range(1, 3)]]`                    | **yes**, via `A(1, Eigen::placeholders::all)` |
| Subarray w/well def layout  | **yes** (strided layout)                                        | no                                                        | **yes** (strided layout)                      | **yes** (strided) |
| Recursive subarray          | **yes** (layout is stack-based and owned by the view)           | **yes** (?)                                               | no (subarray may dangle layout, design bug?)  | **yes** (?) (1D only) |
| Custom Alloctors            | **yes**, via `multi::array<T, D, Alloc>`                        | no (no allocation or ownership)                           | **yes** (stateless?)                          | no | 
| PMR Alloctors               | **yes**, via `multi::pmr::array<T, D>`                          | no (no allocation or ownership)                           |   no     | no |
| Fancy pointers / references | **yes**, via `multi::array<T, D, FancyAlloc>` or views          | no                                                        |   no     | no |
| Strided Layout              | **yes**                                                         | **yes**                                                   |  **yes** | **yes** |
| Fortran-ordering            | **yes**, only for views, e.g. resulted from transposed views    | **yes** (only views are supported)                        |  **yes** | **yes** |
| Zig-zag / Hilbert ordering  | no                                                              | **yes**, via arbitrary layouts (no inverse or flattening) | no       | no |
| Arbitrary layout            | no                                                              | **yes**, possibly inneficient, no efficient slicing       | no       | no |
| Flattening of elements      | **yes**, via `A.elements()` range (efficient representation)    | **yes**, but via indices roundtrip (inefficient)          | no, only for allocated arrays | no, not for subblocks (?) |
| Iterators                   | **yes**, standard compliant, random-access-iterator             | no, or very limited                                       | **yes**, limited | no |
| Multidimensional iterators (cursors) | **yes** (experimental)                                 | no                                                        | no               | no |         
| STL algorithms or Ranges    | **yes**                                                         | no, limited via `std::cartesian_product`                  | **yes**, some do not work | no |
| Compatibility with Boost    | **yes**, serialization, interprocess  (see below)               | no                                                        | no | no |
| Compatibility with Thrust or GPUs | **yes**, via flatten views (loop fusion), thrust-pointers/-refs | no                                                   | no          | no |
| Used in production          | [QMCPACK](https://qmcpack.org/), [INQ](https://gitlab.com/npneq/inq)  | (?) , experience from Kokkos incarnation             | **yes** (?) | [**yes**](https://eigen.tuxfamily.org/index.php?title=Main_Page#Projects_using_Eigen) |


## Serialization

The capability of serializing arrays is important to save/load data to/from disk and also to communicate values via streams or networks (including MPI).
The C++ language does not give any facilities for serialization and unfortunately the standard library doesn't either.

However there are a few libraries that offer a certain common protocol for serialization,
such as [Boost.Serialization](https://www.boost.org/doc/libs/1_76_0/libs/serialization/doc/index.html) and [Cereal](https://uscilab.github.io/cereal/).
The Multi library is compatible with both of them, and yet it doesn't depend on any of them.
The user can choose one or the other, or none if serialization is not needed.
The generic protocol is such that variables are (de)serialized using the (`>>`)`<<` operator with the archive; operator `&` can be used to have single code for both.
Serialization can be binary (efficient) or text-based (human readable).

Here it is a small implementation of save and load functions for array to JSON format with Cereal.
The example can be easily adapted to other formats or libries (XML with Boost.Serialization are commented on the right).

```cpp
#include <multi/array.hpp>  // this library

#include <cereal/archives/json.hpp>                // #include <boost/archive/xml_iarchive.hpp>
                                                   // #include <boost/archive/xml_oarchive.hpp>
// for serialization of array elements (in this case strings)
#include <cereal/types/string.hpp>                 // #include <boost/serialization/string.hpp>

#include<fstream>  // saving to files in example

using input_archive  = cereal::JSONInputArchive ;  // boost::archive::xml_iarchive;
using output_archive = cereal::JSONOutputArchive;  // boost::archive::xml_oarchive;
using cereal::make_nvp;                            // boost::serialization::make_nvp;

namespace multi = boost::multi;

template<class Element, multi::dimensionality_type D, class IStream> 
auto array_load(IStream&& is) {
	multi::array<Element, D> value;
	input_archive{is} >> make_nvp("value", value);
	return value;
}

template<class Element, multi::dimensionality_type D, class OStream>
void array_save(OStream&& os, multi::array<Element, D> const& value) {
	output_archive{os} << make_nvp("value", value);
}

int main() {
	multi::array<std::string, 2> const A = {{"w", "x"}, {"y", "z"}};
	array_save(std::ofstream{"file.string2D.json"}, A);  // use std::cout to print serialization to the screen

	auto const B = array_load<std::string, 2>(std::ifstream{"file.string2D.json"});
	assert(A == B);
}
```
[(online)](https://godbolt.org/z/9j9avjh8M)

These templated functions work for any dimension and element type (as long as the element type is serializable in itself; all basic types are serializable by default).
However, note that the user must ensure that data is serialized and deserialized into the same type;
the underlying serialization libraries only do minimal consistency checks for efficiency reasons and don't try to second-guess file formats or contained types.
Serialization is a relatively low-level feature for which efficiency and economy of bytes is a priority.
Cryptic errors and crashes can occur if serialization libraries, file formats, or C++ types are mixed between writes and reads.

References to subarrays (views) can also be serialized; however, size information is not saved in such cases.
The reasoning is that references to subarrays cannot be resized in their number of elements if there is a size mismatch during deserialization.
Therefore, array views should be deserialized as other array views, with matching sizes.

The output JSON file of the previous example looks like this.
(The XML would have a similar structure.)

```json
{
    "value": {
        "cereal_class_version": 0,
        "extensions": {
            "cereal_class_version": 0,
            "extension": {
                "cereal_class_version": 0,
                "first": 0,
                "last": 2
            },
            "extension": {
                "first": 0,
                "last": 2
            }
        },
        "elements": {
            "cereal_class_version": 0,
            "item": "w",
            "item": "x",
            "item": "y",
            "item": "z"
        }
    }
}
```

Large datasets tend to be serialized slowly for archives with heavy formatting.
Here it is a comparison of speeds when (de)serializing a 134 MB 4-dimensional array of with random `double`s.

| Archive format (Library)     | file size     | speed (read - write)           | time (read - write)   |
| ---------------------------- | ------------- | ------------------------------ |-----------------------|
| JSON (Cereal)                | 684 MB        |    3.9 MB/sec -    8.4 MB/sec  |  32.1 sec - 15.1  sec |
| XML (Cereal)                 | 612 MB        |    2.  MB/sec -    4.  MB/sec  |  56   sec  - 28   sec |
| XML (Boost)                  | 662 MB        |   11.  MB/sec -   13.  MB/sec  |  11   sec  -  9   sec |
| YAML ([custom archive)](https://gitlab.com/correaa/boost-archive-yml)             | 702 MB        |   10.  MB/sec -    4.4 MB/sec  |  12   sec  - 28   sec |
| Portable Binary (Cereal)     | 134 MB        |  130.  MB/sec -  121.  MB/sec  |  9.7  sec  - 10.6 sec |
| Text (Boost)                 | 411 MB        |   15.  MB/sec -   16.  MB/sec  |  8.2  sec  - 7.6  sec |
| Binary (Cereal)              | 134 MB        |  134.4 MB/sec -  126.  MB/sec  |  0.9  sec  -  0.9 sec |
| Binary (Boost)               | 134 MB        | 5200.  MB/sec - 1600.  MB/sec  |  0.02 sec -   0.1 sec |
| gzip-XML (Cereal)            | 191 MB        |    2.  MB/sec -    4.  MB/sec  | 61    sec  - 32   sec |
| gzip-XML (Boost)             | 207 MB        |    8.  MB/sec -    8.  MB/sec  | 16.1  sec  - 15.9 sec |

## Range-v3

The library works out of the box with Eric Niebler's Range-v3 library.
The library helps removing explicit iterators (e.g. `begin`, `end`) from the code when possible.

Every Multi array object can be regarded as range.
Every subarray references (and array values) are interpreted as range views.

For example for a 2D array `d2D`, `d2D` itself is interpreted as a range of rows.
Each row, in turn, is interpreted as a range of elements.
In this way, `d2D.transposed()` is interpreted as a range of columns (of the original array), and each column a range of elements (arranged vertically in the original array).

```cpp
#include <range/v3/all.hpp>
int main(){

	multi::array<int, 2> const d2D = {
		{ 0,  1,  2,  3}, 
		{ 5,  6,  7,  8}, 
		{10, 11, 12, 13}, 
		{15, 16, 17, 18}
	};
	assert( ranges::inner_product(d2D[0], d2D[1], 0.) == 6+2*7+3*8 );
	assert( ranges::inner_product(d2D[0], rotated(d2D)[0], 0.) == 1*5+2*10+15*3 );

	static_assert(ranges::RandomAccessIterator<multi::array<double, 1>::iterator>{});
	static_assert(ranges::RandomAccessIterator<multi::array<double, 2>::iterator>{});
}
```

In this other [example](https://godbolt.org/z/MTodPEnsr), a 2D Multi array (or subarray) is modified such that each element of a column is subtracted the mean value of such column.

```cpp
#include<multi/array.hpp>
#include<range/v3/all.hpp>

template<class MultiArray2D>
void subtract_mean_columnwise(MultiArray2D&& arr) {
    auto&& tarr = arr.transposed();
    auto const column_mean = 
        tarr
        | ranges::views::transform([](auto const& row) {return ranges::accumulate(row, 0.0)/row.size();})
        | ranges::to<multi::array<double, 1>>
    ;

    ranges::transform(
        arr.elements(),
        column_mean | ranges::views::cycle,
        arr.elements().begin(),
        [](auto const elem, auto const mean) {return elem - mean;}
    );
}
```

## Boost.Interprocess

Using Interprocess allows for shared memory and for persistent mapped memory.

```cpp
#include <boost/interprocess/managed_mapped_file.hpp>
#include "multi/array.hpp"
#include<cassert>

namespace bip = boost::interprocess;
using manager = bip::managed_mapped_file;
template<class T> using mallocator = bip::allocator<T, manager::segment_manager>;
auto get_allocator(manager& m){return m.get_segment_manager();}

namespace multi = boost::multi;
template<class T, int D> using marray = multi::array<T, D, mallocator<T>>;

int main(){
{
	manager m{bip::create_only, "bip_mapped_file.bin", 1 << 25};
	auto&& arr2d = *m.construct<marray<double, 2>>("arr2d")(std::tuple{1000, 1000}, 0., get_allocator(m));
	arr2d[4][5] = 45.001;
	m.flush();
}
{
	manager m{bip::open_only, "bip_mapped_file.bin"};
	auto&& arr2d = *m.find<marray<double, 2>>("arr2d").first;
	assert( arr2d[4][5] == 45.001 );
	m.destroy<marray<double, 2>>("arr2d");//    eliminate<marray<double, 2>>(m, "arr2d");}
}
}
```

(Similarly works with [LLNL's Meta Allocator](https://github.com/llnl/metall))

## CUDA Thrust (and HIP Thrust)

The library works out-of-the-box in combination with the Thrust library.

```cpp
#include <multi/array.hpp>  // this library

#include <thrust/device_allocator.h>  // from CUDA or ROCm distributions

namespace multi = boost::multi;

int main() {
	multi::array<double, 2, thrust::device_allocator<double>> A({10,10});
	multi::array<double, 2, thrust::device_allocator<double>> B({10,10});
	A[5][0] = 50.0;

	thrust::copy(A.rotated()[0].begin(), A.rotated()[0].end(), B.rotated()[0].begin());  // copy row 0
	assert( B[5][0] == 50.0 );
}
```
[(live)](https://godbolt.org/z/e7bjKqh69)

which uses the default Thrust device backend (i.e. CUDA when compiling with `nvcc`, HIP/ROCm when compiling with a HIP/ROCm compiler, or OpenMP or TBB in other cases).
Universal memory (accessible from normal CPU code) can be used with `thrust::universal_allocator` (from `<thrust/universal_allocator.h>`) instead.

More specific allocators can be used ensure CUDA backends, for example CUDA managed memory:

```cpp
#include <thrust/system/cuda/memory.h>
...
	multi::array<double, 2, thrust::cuda::universal_allocator<double>> A({10,10});
```

In the same way, to *ensure* HIP backends please replace the `cuda` namespace by the `hip` namespace, and in the directory name `<thrust/system/hip/memory.h>`.
`<thrust/system/hip/memory.h>` is provided by the ROCm distribution (in `/opt/rocm/include/thrust/system/hip/`, and not by the NVIDIA distribution.)

Multi doesn't have a dependency on Thrust (or viseversa);
they just work well together, both in terms of semantics and efficiency.
Certain "patches" (to improve Thrust behavior) can be applied to Thrust to gain extra efficiency and achieve near native speed by adding the `#include<multi/adaptors/thrust.hpp>`.

Multi can be used on existing memory in a non-invasive way via (non-owning) reference arrays:

```cpp
	// assumes raw_pointer was allocated with cudaMalloc or hipMalloc
	using gpu_ptr = thrust::cuda::pointer<double>;  // or thrust::hip::pointer<double> 
	multi::array_ref<double, 2, gpu_ptr> Aref({n, n}, gpu_ptr{raw_pointer});
```

Finally, the element type of the device array has to be device-friendly to work correctly; 
this includes all build in types, and classes with basic device operations, such as construction, destruction, and assigment.
They notably do not include `std::complex<T>`, in which can be replaced by the device-friendly `thrust::complex<T>` can be used as replacement.

### OpenMP via Thrust

In an analog way, Thrust can also handle OpenMP (omp) allocations and multi-threaded algorithms of arrays.
The OMP backend can be enabled by the compiler flags `-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_BACKEND_OMP` or by using the explicit `omp` system types: 

```cpp
#include <multi/array.hpp>
#include <thrust/system/omp/memory.h>

namespace multi = boost::multi;

int main() {
	multi::array<double, 2, thrust::omp::allocator<double>> A({10,10});
	multi::array<double, 2, thrust::omp::allocator<double>> B({10,10});

	A[5][0] = 50.0;

    // copy row 0
	thrust::copy(A.rotated()[0].begin(), A.rotated()[0].end(), B.rotated()[0].begin());

	assert( B[5][0] == 50.0 );

	auto C = B;  // uses omp automatically for copying behind the scenes
}
```
https://godbolt.org/z/e3cGbY87r

Compilation might need to link to an omp library, `-fopenmp -lgomp`.

### Thrust memory resources

GPU memory is relative expensive to allocate, therefore any application that allocates and deallocates arrays often will suffer performance issue.
This is where special memory management is important, for example for avoiding real allocations when possible by caching and reusing memory blocks.

Thrust implements both polymorphic and non-polymorphic memory resources via `thrust::mr::allocator<T, MemoryResource>`;
Multi supports both.

```cpp
auto pool = thrust::mr::disjoint_unsynchronized_pool_resource(
	thrust::mr::get_global_resource<thrust::universal_memory_resource>(),
	thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
);

// memory is handled by pool, not by the system allocator
multi::array<int, 2, thrust::mr::allocator<int, decltype(pool)>> arr({1000 - i%10, 1000 + i%10}, &pool);  // or multi::mr::array<int, 2, decltype(pool)> for short
```

The associated pointer type for the array data is deduced from the _upstream_ resource; in this case, `thrust::universal_ptr<int>`.

As as quick recipe to improve performance in many cases, here it is a recipe for a `caching_allocator` which uses a global (one per thread) memory pool.
The requested memory resides in GPU (managed) memory (`thrust::cuda::universal_memory_resource`) while the cache _bookkeeping_ is held in CPU memory (`new_delete_resource`).

```cpp
template<class T, class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>>>
struct caching_allocator : Base_ {
	caching_allocator() : 
		Base_{&thrust::mr::tls_disjoint_pool(
			thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
			thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
		)} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}  // all caching allocator are equal
	template<class U> struct rebind {using other = caching_allocator<U>;};
};
...
int main() {
	...
	using array2D = multi::array<double, 2, caching_allocator<double>>;

	for(int i = 0; i != 10; ++i) { array2D A({100, 100}); ... use A ...}
}
```

In the example, most of the memory requests are handled by reutilizing the memory pool avoiding expensive system allocations.
More targeted usage patterns may require locally (non-globally) defined memory resources.

## CUDA C++

CUDA is a dialect of C++ that allows writing pieces of code directly for GPU execution, known as "CUDA kernels".
CUDA code is generally "low level" (less abstracted) but it can be used in combination with CUDA Thrust or the CUDA runtime library, specially to implement algorithm that are hard to implement otherwise.
Although code inside kernels has certain restrictions, most Multi expressions can be used. 
(Most functions in Multi, except those involving memory allocations, are marked `__device__` to allow this.)

Calling kernels involves a special syntax (`<<< ... >>>`), and they cannot take arguments by reference (or by values that are not trivial, e.g. not entirely contained in the stack).
Since arrays are usually passed by reference (e.g. `multi::array<double, 2>&` or `Array&&`), a different idiom needs to be used.
(Large arrays are not passed by value to avoid copies, but even if a copy would be fine, kernel arguments cannot allocate memory themselves.)
Iterators (e.g. `.begin()/.end()`) and "cursors" (e.g. `.home()`) are "trivial to copy" and can be passed by value and represent a "proxy" to an array, including allowing the normal index syntax and other transformations.

Cursors are a generalization of iterators for multiple dimensions.
They are cheaply copied (like iterators) and they allow indexing.
Also, they have no associated `.size()` or `.extensions()`, but this is generally fine for kernels.

Here it is an example implementation for matrix multiplication, in combination with Thrust and Multi,

```cpp
#include <multi/array.hpp>  // from https://gitlab.com/correaa/boost-multi
#include <thrust/system/cuda/memory.h>  // for thrust::cuda::allocator

template<class ACursor, class BCursor, class CCursor>
__global__ void Kernel(ACursor A, BCursor B, CCursor C, int N) {
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	typename CCursor::element_type value{0.0};
	for (int k = 0; k != N; ++k) { value += A[y][k] * B[k][x]; }
	C[y][x] = value;
}

namespace multi = boost::multi;

int main() {
	int N = 1024;

	// declare 3 square arrays
	multi::array<double, 2, thrust::cuda::allocator<double>> A({N, N}); A[0][0] = ...;
	multi::array<double, 2, thrust::cuda::allocator<double>> B({N, N}); B[0][0] = ...;
	multi::array<double, 2, thrust::cuda::allocator<double>> C({N, N});

	// kernel invocation code
	assert(N % 32 == 0);
	dim3 dimBlock(32, 32);
	dim3 dimGrid(N/32, N/32);
	Kernel<<<dimGrid, dimBlock>>>(A.home(), B.home(), C.home(), N);
	cudaDeviceSynchronize();

    // now C = A x B
}
```
[(live)](https://godbolt.org/z/eKbeosrWa)

## TotalView

TotalView visual debugger (commercial), popular in HPC environments, can display arrays in human-readable form (for simple types, like `double` or `std::complex`).
To use it, simply `#include "multi/adaptors/totalview.hpp"` and link to the TotalView libraries, compile and run the code with the TotalView debugger.

# Technical points

### What's up with the multiple bracket notation? 

The chained bracket notation (`A[i][j][k]`) allows to refer to elements and subarrays lower dimensional subarrays in a consistent and _generic_ manner and it is the recommended way to access the array objects.
It is a frequently raised question whether the chained bracket notation is good for performance, since it appears that each utilization of the bracket leads to the creation of a temporary object which in turn generates a partial copy of the layout.
Moreover, this goes against [historical recommendations](https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op).

It turns out that modern compilers with a fair level of optimization (`-O2`) can elide these temporary objects, so that `A[i][j][k]` generates identical assembly code as `A.base() + i*stride1 + j*stride2 + k*stride3` (+offsets not shown).
In a subsequent optimization, constant indices can have their "partial stride" computation removed from loops. 
As a result, these two loops lead to the [same machine code](https://godbolt.org/z/ncqrjnMvo):

```cpp
	// given the values of i and k and accumulating variable acc ...
    for(long j = 0; j != M; ++j) {acc += A[i][j][k];}
```
```cpp
    auto* base = A.base() + i*std::get<0>(A.strides()) + k*std::get<2>(A.strides());
    for(long j = 0; j != M; ++j) {acc += *(base + j*std::get<1>(A.strides()));}
```

Incidentally, the library also supports parenthesis notation with multiple indices `A(i, j, k)` for element or partial access, but it does so as part of a more general syntax to generate sub-blocks.
In any case `A(i, j, k)` is expanded to `A[i][j][k]` internally in the library when `i, j, k` are normal integer indices.
Additionally, array coordinates can be directly stored in tuple-like data structures, allowing this functional syntax:

```cpp
std::array p = {2, 3, 4};
std::apply(A, p) = 234;  // same as A(2, 3, 4) = 234; and same as A[2][3][4] = 234;
```

### Customizing recursive operations: SCARY iterators

A level of customization can be achieved by intercepting internal recursive algorithms.
Multi iterators are [SCARY](http://www.open-std.org/jtc1/sc22/WG21/docs/papers/2009/n2980.pdf). 
SCARY means that they are independent of any container and can be accessed generically through their dimension and underlying pointer types:

For example, `boost::multi::array_iterator<double, 2, double*> it` is a row (or column) iterator of an array of dimension 2 or higher, whose underlying pointer type is `double*`.
This row (or column) and subsequent ones can be accessed by the normal iterator(pointer) notation `*it` and `it[n]` respectively.
Indirection `it->...` is supported (even for iterators if high dimension). 
The base pointer, the strides and the size of the arrow can be accessed by `base(it)`, `stride(it)`, `it->size()`.

The template arguments of the iterator can be used to customize operations that are recursive (and possibly inefficient in certain context) in the library:

```cpp
namespace boost{namespace multi{
template<class It, class T>  // custom copy 1D (aka strided copy)
void copy(It first, It last, multi::array_iterator<T, 1, fancy::ptr<T> > dest){
	assert( stride(first) == stride(last) );
	std::cerr<<"1D copy(it1D, it1D, it1D) with strides "<< stride(first) <<" "<< stride(dest) <<std::endl;
}

template<class It, class T> // custom copy 2D (aka double strided copy)
void copy(It first, It last, multi::array_iterator<T, 2, fancy::ptr<T> > dest){
	assert( stride(first) == stride(last) );
	std::cerr<<"2D copy(It, It, it2D) with strides "<< stride(first) <<" "<< stride(dest) <<std::endl;
}
}}
```

For example, if your custom pointers refers a memory type in which 2D memory copying (strided copy) is faster than sequencial copying, that kind of instruction can be ejecuted when the library internally calls `copy`.
This customization must be performed (unfortunately) in the `boost::multi` namespace (this is where the Multi iterators are defined) and the customization happens through matching the dimension and the pointer type.
