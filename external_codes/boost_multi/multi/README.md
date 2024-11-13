<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.] Multi

> **Disclosure: This is not an official or accepted Boost library and is unrelated to the std::mdspan proposal. It is in the process of being proposed for inclusion in [Boost](https://www.boost.org/).**

_© Alfredo A. Correa, 2018-2024_

_Multi_ is a modern C++ library that provides manipulation and access of data in multidimensional arrays for both CPU and GPU memory.

Multidimensional array data structures are fundamental to several branches of computing, such as data analysis, image processing, and scientific simulations, and in combination with GPUs to Artificial Intelligence and Machine Learning.
This library offers array containers and subarrays in arbitrary dimensions with well-behaved value semantics,
featuring logical access recursively across dimensions and to elements through indices and iterators.

The internal data structure layout is stride-based, which makes it compatible with low-level C libraries.

The library interface is designed to be compatible with standard algorithms and ranges (STL) and special memory (including GPUs) and follows modern C++ design principles.

Features of this library that aim to facilitate the manipulation of multidimensional arrays include:

* Value semantics of multidimensional array containers and well-defined referential semantics to avoid unnecessary copies if possible.
* Availability of different access patterns to the elements in the multidimensional structure, as nested sequences or as a single sequence of elements.
A D-dimensional array can be interpreted either as an (STL-compatible) sequence of (D-1)-dimensional subarrays or as a flattened one-dimensional (also STL-compatible) sequence of elements.
* Interoperability with both legacy C and modern C++ libraries (e.g., STL, ranges, Thrust --CUDA and AMD GPUs--, Boost).
* Memory management and allocation to exploit modern memory spaces, including GPU memory, mapped memory, and fancy pointers.

Do not confuse this library with [Boost.MultiArray](https://www.boost.org/doc/libs/1_69_0/libs/multi_array/doc/index.html) 
or with the standard MDSpan proposal `std::mdspan`.
This library shares some of their goals and is compatible with them, but it is designed at a different level of generality and with other priorities (such as the features listed above).
The code is entirely independent and has fundamental implementation and semantics differences.

The library's primary concern is with the storage and logic structure of data;
it doesn't make algebraic or geometric assumptions about the arrays and their elements.
(It is still a good building block for implementing mathematical algorithms, such as representing algebraic dense matrices in the 2D case.)

The library does not throw exceptions and provides basic guarantees (such as no memory leaks) in their presence (e.g., thrown from allocations).
Indexing and other logical errors result in undefined behavior, which this library attempts to reflect via assertions.

The library requires C++17 or higher.

**Contents:**

[[_TOC_]]

## Using the library, installation and tests

Before using the library, you can try it [online](https://godbolt.org/z/dvacqK8jE).

_Multi_ doesn't require installation; a single header `#include <multi/array.hpp>` is enough to use the entire core library.
_Multi_ has no dependencies (except for the standard C++ library) and can be used immediately after downloading.

```bash
git clone https://gitlab.com/correaa/boost-multi.git
```

Although installation is unnecessary, the library can still be installed with CMake.
The header (and CMake) files will be installed in the chosen prefix location (by default, `/usr/local/include/multi` and `/usr/local/share/multi`).

```bash
cd boost-multi
mkdir -p build && cd build
cmake ..  # --install-prefix=$HOME/.local
cmake --install .  # or sudo ...
```

_Testing_ the library requires the Boost.Test library, installed for example, via `sudo apt install cmake git g++ libboost-test-dev make` or `sudo dnf install boost-devel cmake gcc-c++ git`.
A CMake build system is provided to compile and run basic tests.

```bash
cmake --build .
ctest
```

Once installed, other CMake projects (targets) can depend on Multi by adding a simple `add_subdirectory(my_multi_path)` or by `find_package`:

```cmake
find_package(multi)  # see https://gitlab.com/correaa/boost-multi#using-the-library-installation-and-tests
```

As an alternatively to using `find_package`, the library can be fetched on demand:
```cmake
include(FetchContent)
FetchContent_Declare(multi GIT_REPOSITORY https://gitlab.com/correaa/boost-multi.git)
FetchContent_MakeAvailable(multi)
...
target_link_libraries(my_target PUBLIC multi)
```

The code requires compilers with standard C++17 support; for reference any of:
LLVM's       `clang` [(5.0+)](https://godbolt.org/z/51E1hjfnn) (`libc++` and `libstdc++`),
GNU's        `g++` [(7.1+)](https://godbolt.org/z/1nGEbKc5a),
Nvidia's    [`nvcc`](https://godbolt.org/z/abdT73PqM) (11.4+) 
and 
            [`nvc++`](https://godbolt.org/z/6z39PjT47) (22.7+),
Intel's      `icpx` (2022.0.0+),
Baxter's    [`circle`](https://www.circle-lang.org/) (build 202+),
[Zig](https://zig.news/kristoff/compile-a-c-c-project-with-zig-368j) in [c++ mode (v0.9.0+)](https://godbolt.org/z/cKGebsWMG),
Edison Desing's [EDG]() [(6.5+)](https://godbolt.org/z/693fxPedx)
and
Microsoft's [MSVC](https://visualstudio.microsoft.com/vs/features/cplusplus/) ([+14.1](https://godbolt.org/z/vrfh1fxWK)).

(Multi code inside CUDA kernel can be compiled with `nvcc` and with [`clang` (in CUDA mode)](https://godbolt.org/z/7dTKdPTxc).
Inside HIP code, it can be compiled with AMD's clang rocm (5.0+).)

Optional "adaptor" sublibraries (included in `multi/adaptors/`) have specific dependencies: fftw, , lapack, thurst, or CUDA
(all of them can be installed with `sudo apt install libfftw3-dev lib64-dev liblapack64-dev libthrust-dev nvidia-cuda-dev` or `sudo dnf install -devel fftw-devel ...`.)

In the following sections we present basic and advanced uses of the libraries. 
Feel free to jump to the "Reference of fundamental types" section to explore a more exhaustive description of the classes provided by the library.

## Basic Usage

The following code declares an array by specifying the element type and the dimensions;
individual elements can be initialized from a nested rectangular list.
```cpp
multi::array<double, 2> A = {
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0},
};

auto const [n, m] = A.sizes();

assert( n == 2 );  // or std::get<0>(A.sizes()) == 2
assert( m == 3 );  // or std::get<1>(A.sizes()) == 3

assert( A.size() == 2 );  // size in first dimension, same as std::get<0>(A.sizes())
assert( A.num_elements() == 6 );  // total number of elements
```

The value of an array can be copied, (moved,) and compared;
copies are equal but independent (disjoint).

```cpp
std::array<double, 2> B = A;
assert(  B       ==  A                 );  // copies are equal
assert( extensions(B) == extensions(A) );  // extensions (sizes) are equal
assert(  B[0][1] ==  A[0][1]           );  // all elements are equal
assert( &B[0][1] != &A[0][1]           );  // elements are independent (dfferent addresses)
```

Individual elements can be accessed by the multidimensional indices, either with square brackets (one index at a time, as above) or with parenthesis (comma separated).

```cpp
assert( &A(1, 2) ==  &A[1][2] );
```

An array can be initialized from its sizes alone, in which case the element values are defaulted (possibly uninitialized):

```cpp
multi::array<double, 3> C({3, 4, 5});
assert( num_elements(C) == 3*4*5 );   // 60 elements with unspecified values
```

Arrays can be passed by value or by reference.
Most of the time, arguments should be passed through generic parameters to also allow functions to work with parts (subblocks, slices, etc.) of an array.
The most useful functions work on the _concept_ of an array rather than on a concrete type, for example:

```cpp
template<class ArrayDouble2D>  // instead of the overspecific argument std::array<double, 2>
auto element_1_1(ArrayDouble2D const& m) -> double const& { return m[1][1]; }
...
assert( &element_1_1(A) == &A[1][1] );
```

The function expects any array or subarray of dimension 2 and returns an element with type `double`.

The generic function template arguments that are not intended to be modified are passed by `const&`; otherwise, they are passed by forward-reference `&&`.
In this way, the functions can be applied to subblocks of larger matrices.

```cpp
assert( &element_1_1(C3D[0]) == &C3D[0][1][1] );
```

(Although most of the examples use numeric elements for conciseness, the library is designed to hold general types (e.g. non-numeric, non-trivial types, like `std::string`, other containers or, in general, user-defined value-types.)

## Advanced Usage

In this example, we are going to use memory that is not managed by the library and manipulate the elements.
We can create a static C-array of `double`s, and refer to it via a bidimensional array `multi::array_ref<double, 2>`.

```cpp
#include <multi/array.hpp>

#include <algorithm>  // for sort
#include <iostream>  // for print

namespace multi = boost::multi;

int main() {
	double d_data[20] = {
		150.0, 16.0, 17.0, 18.0, 19.0,
		 30.0,  1.0,  2.0,  3.0,  4.0,
		100.0, 11.0, 12.0, 13.0, 14.0,
		 50.0,  6.0,  7.0,  8.0,  9.0
	};  // block of 20 elements ...
	multi::array_ref<double, 2> d2D_ref{&d_data[0], {4, 5}};  // .. interpreted as a 4 by 5 array
	...
```

Next, we print the elements in a way that corresponds to the logical arrangement:

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
Presumably, if one can sort over a range, one can perform any other standard algorithm.

```cpp
		...
		std::stable_sort( d2D_ref.begin(), d2D_ref.end() );
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

(Note that `std::sort` cannot be applied directly to a multidimensional C-array or to other libraries, such as Boost.MultiArray.
The arrays implemented by this library are, to the best of my knowledge, the only ones that support all STL algorithms directly.)

If we want to order the matrix on a per-column basis, we need to "view" the matrix as a range of columns.
This is done in the bidimensional case, by accessing the matrix as a range of columns:

```cpp
		...
		std::stable_sort( rotated(d2D_ref).begin(), rotated(d2D_ref).end() );
	}
```

The `rotate` operation rotates indices, providing a new logical view of the original array without modifying it.

In this case, the original array will be transformed by sorting the matrix into:

> ```
> 1 2 3 4 30
> 6 7 8 9 50
> 11 12 13 14 100
> 16 17 18 19 150
> ```

By combining index rotations and transpositions, an array of dimension `D` can be viewed simultaneously as `D!` (D-factorial) different ranges of different "transpositions" (rotation/permutation of indices.)

## Initialization

`array_ref` is initialized from a preexisting contiguous range, the index extensions should be compatible with the total number of elements.

```cpp
double* dp = new double[12];
multi::array_ref<double, 2> A({3, 4}, dp);
multi::array_ref<double, 2> B({2, 6}, dp);
...
delete[] dp;
```

Array references do not own memory and, just as language references, can not be rebinded (i.e. resized or "reseated") to refer to a different memory location.
Since `array_ref` is an array reference, it can "dangle" if the original memory is deallocated.

Array objects (`multi::array`), in contrast, own the elements they contain and can be resized later.
An `array` is initialized by specifying the index extensions and, optionally, a default value).

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

## Copy and assigment (and aliasing)

The library offers value semantics for the `multi::array<T, D>` family of classes.
Constructing or assigning from an existing array generates a copy of the original object, independent of the original one but equal in value.

```cpp
auto B2 = A2;  // same as multi::array<double, 2> B2 = A2; (A2 is defined above)

assert(  B2       ==  A2       );  // copies have the same value (and also the same shape)
assert(  B2[0][0] ==  A2[0][0] )
assert( &B2[0][0] != &A2[0][0] );  // but they are independent
```

A (mutable) array can be assigned at any moment, independently of the previous state or shape (extensions).
The dimensionalities must match.
```cpp
B2 = A2;  // both have dimensionality 2
```

Sometimes it is necessary to generate copies from views or subblocks.
```cpp
multi::array<double, 3> C2 = A2( {0, 2}, {0, 2} );
```
or equivalently,
```cpp
auto C2 = + A2( {0, 2}, {0, 2} );
```
Note the use of the prefix `+` as an indicator that a copy must be created (it has no arithmetic implications).
Due to a language limitation, omitting the `+`` will create another non-independent reference view of the left-hand side, which is generally undesired.

Subarray-references can also assigned, but only if the shapes of the left-hand side (LHS) and right-hand side (RHS) match.
Otherwise, the behavior is undefined (in debug mode, the program will fail an assertion).

```cpp
C2( {0, 2}, {0, 2} ) = A2( {0, 2}, {0, 2} );  // both are 2x2 views of arrays, *elements* are copied
```

Using the same or overlapping arrays in the RHS and LHS of assignment produces undefined behavior in general (and the library doesn't check).
Notably, this instruction does not transpose the array but produces an undefined result:

```cpp
A2 = A2.transposed();
```

While this below instead does produce a transposition, at the cost of making one copy (implied by `+`) of the transposed array first and assigning (or moving) it back to the original array.

```cpp
A2 = + A2.transposed();
```

This is an instance of the problem of _data aliasing_, which describes a common situation in which a data location in memory can be accessed through different parts of an expression or function call.
Within the confines of the library interface, this pitfall can only occur on assignment as illustrated above.

However the problem of aliasing can persist when taking mutable array-references in function arguments.
The most general solution to this problem is to make copies or directly work with completely disjoint objects; 
but other case-by-case solutions might be possible.

For example, in-place transposition is an active subject of research;
_optimal_ speed and memory transpositions might require specially designed libraries.

Finally, arrays can be efficiently moved by transferring ownership of the internal data.

```cpp
auto B2 = std::move(A2);  // A2 is empty after this
```

Subarrays do not own the data; therefore they cannot directly take advantage of this feature.
However, individual elements of a view can be moved; this is particularly useful if the elements are expensive to copy.
A "moved" subview is simply another kind of view of the elements.

```cpp
multi::array<std::vector<double>, 2> A({10, 10}, std::vector<double>(1000));
multi::array<std::vector<double>, 2> B({10, 10});
...
B[1] = A[2].element_moved();
```

Each of the 10 *elements* of the third row of A is moved into the second row of B. A[2] still has 10 (moved-from) empty vectors


## Change sizes (extents)

Arrays can change their size while _preserving elements_ with the `reextent` method.

```cpp
multi::array<int, 2> A = {
	 {1, 2, 3},
	 {4, 5, 6}
};

A.rextent({4, 4});

assert( A[0][0] == 1 );
```

Arrays can be emptied (to zero-size) with `.clear()` (equivalent to `.rextent({0, 0, ...})`).

The main purpose of `reextent` is element preservation.
Allocations are not amortized; 
except for trivial cases, all calls to reextend allocate and deallocate memory.
If element preservation is not desired, a simple assignment (move) from a new array expresses the intention better and it is more efficient since it doesn't need to copy preexisting elements.

```cpp
A = multi::array<int, 2>({4, 4});  // like A.rextent({4, 4}) but elements are not preserved.
```

An alternative syntax, `.rextent({...}, value)` sets _new_ (not preexisting) elements to a specific value.

Subarrays or views cannot change their size or be emptied (e.g. `A[1].rextent({4})` or `A[1].clear()` will not compile).
For the same reason, subarrays cannot be assigned from an array or another subarray of a different size.

Changing the size of arrays by `reextent`, `clear`, or assignment generally invalidates existing iterators and ranges/views.

## Iteration (range-based loops vs iterators)

Historically, iteration over arrays has been done with `for`-loops, where each nesting level is associated with each dimension.
The valid range of indices in all the dimensions of an array is extracted with `.extensions()`.
In the 2D case, `.extensions()` can be conveniently decomposed into two ranges, one for each dimension.

```
	multi::array<int, 2> A = {
		{1, 2, 3},
		{4, 5, 6}
	};

	auto [is, js] = A.extensions();
	for(auto i : is) {  // is == {0, 1} (range from 0 to 2, not included)
		for(auto j : js) {  // ij = {0, 1, 2} (range from 0 to 3, not included)
			A[i][j] *= 2;
		}
	}
```

The elements of the 2D array can be accessed directly without intermediate indices:

```cpp
	for(auto&& row : A) {
		for(int& e: row) {
			e *= 2;
		}
	}
```

However, in some cases it is better to use the iterator-based interface.
The iterator-based interface is more convenient to express and interact with generic algorithms, which in turn can be parallelized and less prone to index errors (such as off-by-one, and out-of-range access.)

Array (and subarray-references) provide a members `.begin()` and `.end()` that produce random-access iterators that access the multidimensional structure through the first dimension (leftmost index).
Accessing arrays by iterators (`begin`/`end`) enables the use of many iterator-based algorithms (see the sort example above).
`begin(A)/end(A)` (or equivalently `A.begin()/A.end()`) gives iterators that are linear and random access in the leading dimension.

As an alternative the elements can be iterated in a flat manner, using the `.elements()` member.
This flattening is done in a canonical order (rightmost index changes fastest) and it is provided whether the elements are contiguous or not in memory.
This "elements" range also provides the begin and end iterators (`.elements().begin()`).

Other non-leading dimensions can be obtained by "rotating" indices first.
`A.rotated().begin()/.end()` gives access to a range of subarrays in the second dimension number (the first dimension is put at the end).
(`.cbegin()/.cend()` give constant (read-only) access.)

As an example, this function allows printing arrays of arbitrary dimensionality into a linear comma-separated form.

```cpp
void recursive_print(double const& d) { cout<<d; };  // terminating overload

template<class Array>
void recursive_print(Array const& ma) {
	cout << "{";
	if(! ma.empty()) {
		flat_print(*ma.begin());  // first element
		std::for_each(ma.begin() + 1, ma.end(), [](auto const& e) { cout<<", "; flat_print(e);});  // rest
	}
	cout << "}";
}
...
recursive_print(A);
```
> ```
> {{{1.2, 1.1}, {2.4, 1}}, {{11.2, 3}, {34.4, 4}}, {{15.2, 99}, {32.4, 2}}}
> ```

Except for those corresponding to the one-dimensional case, dereferencing iterators generally produce "proxy"-references (i.e. objects that behave in a large degree like language references).
These references can be given a name; using `auto` can be misleading since the resulting variable does not have value semantics.

```cpp
auto row = *A.begin();  // accepted by the language but misleading, row is not an independent value
```

In my experience, however, the following usage pattern produces a more consistent idiom for generating references (still without copying elements):

```cpp
auto&&       row0 = *A.begin() ;  // same as decltype(A)::      reference  row0 = * begin(A);
auto const& crow0 = *A.cbegin();  // same as decltype(A)::const_reference crow0 = *cbegin(A);

auto&&       row1 =               A [1];  // same as decltype(A)::      reference  row1 =               A [1];
auto const& crow1 = std::as_const(A)[1];  // same as decltype(A)::const_reference crow0 = std::as_const(A)[1];
```

If a new value is desired, these (equivalent) options express the intention more explicitly:

```cpp
decltype(A)::value_type row =   *begin(A);  // there is a real copy of the row
                   auto row = + *begin(A);  // there is another copy, note the use of '+' (unary plus)
```

In the examples above all elements are accessed in a nested way, recursively down the dimensions.
To iterate over all the elements regardless of the multidimensional structure the following function can print all the elements.

```cpp
template<class Array>
void flat_print(Array const& ma) {
	cout << "[";
	std::for_each(ma.elements().begin(), ma.elements().end(), [](auto&& e) { cout<< e << ", ";});
	cout << "]";
}
...
recursive_print(A);
```
> ```
> [1.2, 1.1, 2.4, 1, 11.2, 3, 34.4, 4, 15.2, 99, 32.4, 2]
> ```

This feature allows to view the array as a flat sequence using the `.elements()` range, which also has `.begin()`/`.end()` and indexing.
For example array element at indices 1,1 is the same as the element 


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
However, in the library's abstraction, `A[2]` references an existing part of the original array, i.e. it is a "library reference", whose "library address" can be obtained with the `&` operator. 
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

Other notations are available, for example this is equivalent to `A(multi::_ , {1, 3, /*every*/2})` or `~(~A)({1, 3, 2})`.
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

Broadcasting is popular in array-based languages, such as Julia and NumPy, and the broadcast operation is generally applied automatically to match the dimension expected by the operation and other operation inputs.
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

<!-- For illustration purposes only, `fill` here is replaced by `copy`; problematic uses are highlighted:

```cpp
multi::array<double, 2> B({10, 2});
std::fill  (B.begin(), B.end(), b);                                       // canonical way
std::fill_n(B.begin(), B.size(), b);                                      // canonical way

std::copy_n(b.broadcasted().begin(), B.size(), B.begin());                // equivalent, using broadcast

std::copy_n(b.broadcasted().begin(), b.broadcasted().size(), B.begin());  // incorrect, undefined behavior, no useful size()
std::copy  (b.broadcasted().begin(), b.broadcasted().end(), B.begin());   // incorrect, undefined behavior, non-terminating loop (end is not reacheable)
B = b.broadcasted();                                                      // incorrect, undefined behavior, B would be of infinite allocated size
``` -->

Unlike in popular languages, broadcasting is not automatic in the library and is applied to the leading dimension only, one dimension at a time.
Broadcasting in non-leading dimensions can be achieved by transpositions and index rotation.

Abuse of broadcast can make it harder to reason about operations;
its primary use is to reuse existing efficient implementations of algorithms when implementations for a specific lower dimensions are not available.
These algorithms need to be compatible with broadcasted views (e.g., no explicit use of `.size()` or infinite loops stemming from problematic use of `.begin()/end()`.)

(In STL, algorithms ending with `_n` should be friendly to broadcast arrays, unfortunately `std::copy_n` is sometimes internally implemented in terms of `std::copy` causing a problematic iterator arithmetic on infinite arrays.
NB: `thrust::copy_n` can be used instead.)

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

Note that the function `hadamard`, acting on 2D arrays, doesn't use the undefined (infinite) sizes (second dimension of `A` and first dimension of `B`).

## Uninitialized vs. initialized elements

If available, the library can take advantage of trivial initialization for the specific element type.
These types can be primitive or user-defined and come with "trivial default constructors". In simple terms, these constructors are not specified and do nothing, not even set values.

When used in the stack, these types can be declared with no initialization (e.g., `double x;`, the initial value is not well defined or partially-formed) or with initialization (e.g., `double x{};`, same as `double x = 0.0;`).
Analogously, `multi::array` does not initialize individual elements of this kind of type unless specified.

For example, after this construction of the array, the values of the six elements of this array are unspecified (partially-formed).
```cpp
multi::array<int, 2> A2({2, 3});  // A2 elements have unspecified value
```

No behavior of the program should depend on these values. 
(Address sanitizers and memory checkers can detect use of uninitialized values.)
This design is a slight departure from the STL's design, which [eagerly initializes elements in containers](https://lemire.me/blog/2012/06/20/do-not-waste-time-with-stl-vectors/).

If trivial construction is unavailable, the library uses the default initialization.
```cpp
multi::array<std::string, 2> A2({2, 3});  // A2 elements have specified value, the empty value std::string{}
```

For types that afford this partially formed states, elements can be later specified via assignment or assigning algorithms (e.g., copy or transform destination).

Initialization can be enforced by passing a single value argument after the extensions.
```cpp
multi::array<int, 2> A2({2, 3}, 0);  // generically multi::array<T, 2>({2, 3}, T{}); or multi::array<T, 2>({2, 3}, {})
```

This design is particularly advantageous for *numeric* types for which external low-level libraries can fill values.
(or when data sits in GPUs, where the initialization step would require an expensive kernel launch and subsequent synchronization).

Unfortunately, regarding the numeric types, STL's `std::complex<double>` was standardized as not-trivially constructible.
A workaround built-in this library is available by forcing a particular flag on the client code in global scope, for example, immediately after including the library:
```cpp
#include<multi/array.hpp>

template<> inline constexpr
bool multi::force_element_trivial_default_construction<std::complex<double>> = true;  // should be defined as early as possible
```

With this line, `std::complex<double>` elements inside arrays will be left uninitialized unless a value is specified.
The rule will only apply to this library's containers (`multi::array`, etc), and not to other containers (such as `std::vector`) or individual `std::complex` variables.

## Type Requirements

The library design tries to impose the minimum possible requirements over the types that parameterize the arrays.
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

`array`s manage their memory behind the scenes through allocators, which can be specified at construction.
It can handle special memory, as long as the underlying types behave coherently, these include [fancy pointers](https://en.cppreference.com/w/cpp/named_req/Allocator#Fancy_pointers) (and fancy references).
Associated fancy pointers and fancy reference (if any) are deduced from the allocator types.

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
([live](https://godbolt.org/z/oeTss3s35))

(See also, examples of interactions with the CUDA Thrust library to see more uses of special pointer types to handle special memory.)

#### Transformed views

Another kind of use of the internal pointer-like type is to transform underlying values.
These are useful to create "projections" or "views" of data elements.
In the following example a "transforming pointer" is used to create a conjugated view of the elements.
In combination with a transposed view, it can create a hermitic (transposed-conjugate) view of the matrix (without copying elements).
We can adapt the library type `boost::transform_iterator` to save coding, but other libraries can be used also.
The hermitized view is read-only, but with additional work, a read-write view can be created (see `multi::::hermitized` in multi-adaptors).

```cpp
constexpr auto conj = [](auto const& c) {return std::conj(c);};

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

To simplify this boilerplate, the library provides the `.element_transformed(F)` method that will apply a transformation `F` to each element of the array.
In this example, the original array is transformed into a transposed array with duplicated elements.

```cpp
	multi::array<double, 2> A = {
		{1.0, 2.0},
		{3.0, 4.0},
	};

	auto const scale = [](auto x) { return x * 2.0; };

	auto B = + A.transposed().element_transformed(scale);
	assert( B[1][0] == A[0][1] * 2 );
```

([live](https://godbolt.org/z/TYavYEG1T))

Since `element_transformed` is a reference-like object (transformed view) to the original data, it is important to understand the semantics of evaluation and possible allocations incurred.
As mentioned in other sections using `auto` and/or `+` appropriately can lead to simple and efficient expressions.

| Construction    | Allocation of `T`s | Initialization (of `T`s) | Evaluation (of `fun`) | Notes |
| -------- | ------- | ------- | ------- | ------- |
| `multi::array<T, D> const B = A.element_transformed(fun);` | Yes        | No  | Yes | Implicit conversion to `T` if result is different, dimensions must match. B can be mutable.   |
| `multi::array<T, D> const B = + A.element_transformed(fun);` | Yes (and move, or might allocate twice if types don't match)  | No  | Yes | Not recommended | 
| `multi::array<T, D> const B{A.element_transformed(fun)};` | Yes        | No  | Yes | Explicit conversion to `T` if result is different, dimensions must match   |
| `auto const B = + A.elements_transformed(fun);`           | Yes         | No  | Yes | Types and dimension are deduced, result is contiguous, preferred |
| `auto const B = A.element_transformed(fun);`               | No         | No  | No (delayed) | Result is effective a reference, may dangle with `A`, usually `const`, not recommended   |
| `auto const& B = A.elements_transformed(fun);`           | No         | No  | No (delayed) | Result is effective a reference, may dangle with `A`. Preferred way.  |
| `multi::array<T, D> B(A.extensions()); B = A.element_transformed(fun);`           | Yes         | Yes (during construction)  | Yes | "Two-step" construction. `B` is mutable. Not recommended  |

| Assigment    | Allocation of `T`s | Initialization (of `T`s) | Evaluation (of `fun`) | Notes |
| -------- | ------- | ------- | ------- | ------- |
| `B = A.elements_transformed(fun);`           | No, if sizes match | Possibly (when `B` was initialized)  | Yes | `B` can't be declared `const`, it can be a writable subarray, preferred  |
| `B = + A.elements_transformed(fun);`           | Yes | Possibly (when `B` was initialized)  | Yes | Not recommended. |

## Reference documentation of fundamental types

The library interface presents several closely related C++ types (classes) representing arrays.
The fundamental types represent multidimensional containers (called `array`), references that can refer to subsets of these containers (called `subarray`), and iterators.
In addition, there are other classes for advanced uses, such as multidimensional views of existing buffers (called `array_ref`) and non-resizable owning containers (called `static_array`).

When using the library, it is simpler to start from `array`, and other types are rarely explicitly used, especially if using `auto`;
however, it is convenient for documentation to present the classes in a different order since the classes `subarray`, `array_ref`, `static_array`, and `array` have an *is-a* relationship (from left to right). 
For example, `array_ref` has all the methods available to `subarray`, and `array` has all the operations of `array_ref`.
Furthermore, the *is-a* relationship is implemented through C++ public inheritance, so, for example, a reference of type `subarray<T, D>&` can refer to a variable of type `array<T, D>`.

### class `multi::subarray<T, D, P = T* >`

A subarray-reference is part (or a whole) of another larger array.
It is important to understand that `subarray`s have referential semantics, their elements are not independent of the values of the larger arrays they are part of.
To recover value semantics, the elements and structure of a `subarray` can be copied into a new array (type `array`, see later).
An instance of this class represents a subarray with elements of type `T` and dimensionality `D`, stored in memory described by the pointer type `P`.
(`T`, `D`, and `P` initials are used in this sense across the documentation unless indicated otherwise.)

Instances of this class have reference semantics and behave like "language references" as much as possible.
As references, they cannot be rebinded or resized; assignments are always "deep".
They are characterized by a size that does not change.
They are usually the result of indexing over other `subarray`s and `array`s (generally of higher dimensions);
therefore, the library doesn't expose constructors for this class.
The whole object can be invalidated if the original array is destroyed.

| Member types      |                           |
|---                |---                        |
| `value_type`      | `multi::array<T, D - 1, P >` or, for `D == 1`, `iterator_traits<P>::value_type` (usually `T`)   
| `reference`       | `multi::subarray<T, D-1, P >` or, for `D == 1`, `pointer_traits<P>::reference` (usually `T&`) 
| `const_reference` | `multi::const_subarray<T, D-1, P >` or, for `D == 1`, `pointer_traits<P>::rebind<T const>::reference` (usually `T const&`)
| `index`           | indexing type in the leading dimension (usually `std::diffptr_t`)
| `size_type`       | describe size (number of subarrays) in the leading dimension (signed version of pointer size type, usually `std::diffptr_t`)
| `index_range`     | describe ranges of indices, constructible from braced indices types or from an `extension_type`. Can be continuous (e.g. `{2, 14}`) or strided (e.g. `{2, 14, /*every*/ 3}`)
| `extesion_type`   | describe a contiguous range of indices, constructible from braced index (e.g. `{0, 10}`) or from a single integer size (e.g. 10, equivalent to `{0, 10}`). 
| `difference_type` | describe index differences in leading dimension (signed version of pointer size type, usually `std::diffptr_t`)
| `pointer`         | `multi::subarray_ptr<T, D-1, P > or, for `D == 1, `P` (usually `T*`)
| `const_pointer`   | `multi::const_subarray_ptr<T, D-1, P >` or, for `D == 1, `pointer_traits<P>::rebind<T const>` (usually `T const*`)
| `iterator`        | `multi::array_iterator_t<T, D-1, P >`
| `const_iterator`  | `multi::const_array_iterator_t<T, D-1, P >`

| Member fuctions   |    |
|---                |--- |
| (constructors)    | not exposed; copy constructor is not available since the instances are not copyable; destructors are trivial since it doesn't own the elements |
| `size`            | returns the size of the leading dimension |
| `extension`       | returns a range that generates valid indices for the leading dimension, for example `{0, ... size() - 1}` |
| `sizes`           | returns a tuple with the sizes in all dimensions, `std::get<0>(A.sizes()) == A.size()` |
| `extensions`      | returns a tuple of ranges in all dimensions, `std::get<0>(A.extensions()) == A.extension()` |
| `operator=`       | assigns the elements from the source; the sizes must match |

It is important to note that assignments in this library are always "deep," and reference-like types cannot be rebound after construction.
(Reference-like types have corresponding pointer-like types that provide an extra level of indirection and can be rebound (just like language pointers);
these types are `multi::array_ptr` and `multi::subarray_ptr` corresponding to `multi::array_ref` and `multi::subarray` respectively.)

| Relational fuctions       |    |
|---                        |--- |
| `operator==`/`operator!=` | Tells if elements of two `subarray` are equal (and if extensions of the subarrays are the same)
| `operator<`/`operator<=`  | Less-than/less-or-equal      lexicographical comparison (requires elements to be comparable)
| `operator>`/`operator>=`  | Greater-than/grater-or-equal lexicographical comparison (requires elements to be comparable)

It is important to note that, in this library, comparisons are always "deep".
Lexicographical order is defined recursively, starting from the first dimension index and from left to right.
For example, `A < B` if `A[0] < B[0]`, or `A[0] == B[0]` and `A[1] < B[1]`, or ..., etc.
Lexicographical order applies naturally if the extensions of `A` and `B` are different; however, their dimensionalities must match.
(See sort examples).

| Element access    |    |
|---                |--- |
|`operator[]`       | access specified element by index (single argument), returns a `reference` (see above), for `D > 1` it can be used recursively |
|`front`            | access first element (undefined result if array is empty). Takes no argument.
|`back`             | access last element  (undefined result if array is empty). Takes no argument.
|`operator()`       | When used with zero arguments, it returns a `subarray` representing the whole array. When used with one argument, access a specified element by index (return a `reference`) or by range (return a `subarray` of equal dimension). For more than one, arguments are positional and reproduce expected array access syntax from Fortran or Matlab: |

- `subarray::operator()(i, j, k, ...)`, as in `S(i, j, k)` for indices `i`, `j`, `k` is a synonym for `A[i][j][k]`, the number of indices can be lower than the total dimension (e.g., `S` can be 4D).
Each index argument lowers the dimension by one.
- `subarray::operator()(ii, jj, kk)`, the arguments can be indices or ranges of indices (`index_range` member type).
This function allows positional-aware ranges.
Each index argument lowers the rank by one.
A special range is given by `multi::_`, which means "the whole range" (also spelled `multi::all`).
For example, if `S` is a 3D of sizes 10-by-10-by-10, `S(3, {2, 8}, {3, 5})` gives a reference to a 2D array where the first index is fixed at 3, with sizes 6-by-2 referring the subblock in the second and third dimension.
Note that `S(3, {2, 8}, {3, 5})` (6-by-2) is not equivalent to `S[3]({2, 8})({3, 5})` (2-by-10).
- `operator()()` (no arguments) gives the same array but always as a subarray type (for consistency), `S()` is equivalent to `S(S.extension())` and, in turn to `S(multi::_)` or `S(multi::all)`.

| Structure access  | (Generally used for interfacing with C-libraries)   |
|---                |--- |
| `base`            | direct access to underlying memory pointer (`S[i][j]... == S.base() + std::get<0>(S.strides())*i + std::get<1>(S.strides())*j + ...`)
| `stride`          | return the stride value of the leading dimension, e.g `(&A[1][0][0]... - &A[0][0]...)`
| `strides`         | returns a tuple with the strides defining the internal layout
| `layout`          | returns a single layout object with stride and size information |

| Iterators         |    |
|---                |--- |
| `begin/cbegin`    | returns (const) iterator to the beginning
| `end/cend`        | returns (const) iterator to the end

| Capacity          |    |
|---                |--- |
| `sizes`           | returns a tuple with the sizes in each dimension
| `extensions`      | returns a tuple with the extensions in each dimension
| `size`            | returns the number of subarrays contained in the first dimension |
| `extension`       | returns a contiguous index range describing the set of valid indices
| `num_elements`    | returns the total number of elements

| Creating views        | (these operations do not copy elements or allocate)    |
|---                    |---  |
| `broadcasted`         | returns a view of dimensionality `D + 1` obtained by infinite repetition of the original array. (This returns a special kind of subarray with a degenerate layout and no size operation. Takes no argument.)
| `dropped`             | (takes one integer argument `n`) returns a subarray with the first n-elements (in the first dimension) dropped from the original subarray. This doesn't remove or destroy elements or resize the original array 
| `element_transformed` | creates a view of the array, where each element is transformed according to a function (first and only argument) |
| `elements`            | a flatted view of all the elements rearranged canonically. `A.elements()[0] -> A[0][0]`, `A.elements()[1] -> A[0][1]`, etc. The type of the result is not a subarray but a special kind of range. Takes no argument.
| `rotated/unrotated`   | a view (`subarray`) of the original array with indices (un)rotated from right to left (left to right), for `D = 1` returns the same `subarray`. For given `i`, `j`, `k`, `A[i][j][k]` gives the same element as `A.rotated()[j][k][i]` and, in turn the same as `A.unrotated()[k][i][j])`. Preserves dimension. The function is cyclic; `D` applications will give the original view. Takes no argument. |
| `transposed` (same as `operator~`) | a view (`subarray`) of the original array with the first two indices exchanged, only available for `D > 1`; for `D = 2`, `rotated`, `unrotated` and `transposed` give same view. Takes no argument.  |
| `sliced`              | (takes two index arguments `a` and `b`) returns a subarray with elements from index `a` to index `b` (non-inclusive) `{S[a], ... S[b-1]}`. Preserves the dimension.
| `strided`             | (takes one integer argument `s`) returns a subarray skipping `s` elements. Preserves the dimension.

| Creating views by pointer manipulation     |     |
|---                                         |---  |
| `static_cast_array<T2, P2 = T2*>(args...)` | produces a view where the underlying pointer constructed by `P2{A.base(), args...}`. Usually, `args...` is empty. Non-empty arguments are useful for stateful fancy pointers, such as transformer iterators.
| `reinterpret_cast_array<T2>`               | underlying elements are reinterpreted as type T2, element sizes (`sizeof`) have to be equal; `reinterpret_cast_array<T2>(n)` produces a view where the underlying elements are interpreted as an array of `n` elements of type `T2`.

| Creating arrays                     |     |
|---                                  |---  |
| `decay` (same as prefix unary `operator+`) | creates a concrete independent `array` with the same dimension and elements as the view. Usually used to force a value type (and forcing a copy of the elements) and avoid the propagation of a reference type in combination with `auto` (e.g., `auto A2_copy = + A[2];`).

A reference `subarray` can be invalidated when its origin array is invalidated or destroyed.
For example, if the `array` from which it originates is destroyed or resized.

### class `multi::array_ref<T, D, P = T* >`

A D-dimensional view of the contiguous pre-existing memory buffer.
This class doesn't manage the elements it contains, and it has reference semantics (it can't be rebound, assignments are deep, and have the same size restrictions as `subarray`)

Since `array_ref` is-a `subarray`, it inherits all the class methods and types described before and, in addition, it defines these members below.

| Member types      | same as for `subarray` |
|---                |---                        |

| Member functions  | same as for `subarray` plus ... |
|---                |--- |
| (constructors)    | `array_ref::array_ref({e1, e2, ...}, p)` constructs a D-dimensional view of the contiguous range starting at p and ending at least after the size size of the multidimensional array (product of sizes). The default constructor and copy constructor are not exposed. Destructor is trivial since elements are not owned or managed. |

| Element access    | same as for `subarray` |
|---                |--- |

| Structure access  | same as for `subarray` |
|---                |--- |

| Iterators         | same as for `subarray`   |
|---                |--- |

| Capacity          | same as for `subarray`   |
|---                |--- |

| Creating views    | same as for `subarray`  |
|---                |---  |

| Creating arrays   | same as for `subarray`  |
|---                |---  |

| Relational functions   |  same as for `subarray`  |
|---                |--- |

An `array_ref` can be invalidated if the original buffer is deallocated.

### class `multi::static_array<T, D, Alloc = std::allocator<T> >`

A D-dimensional array that manages an internal memory buffer.
This class owns the elements it contains; it has restricted value semantics because assignments are restricted to sources with equal sizes.
Memory is requested by an allocator of type Alloc (standard allocator by default).
It supports stateful and polymorphic allocators, which are the default for the special type `multi::pmr::static_array`.

The main feature of this class is that its iterators, subarrays, and pointers do not get invalidated unless the whole object is destroyed.
In this sense, it is semantically similar to a C-array, except that elements are allocated from the heap.
It can be useful for scoped uses of arrays and multi-threaded programming and to ensure that assignments do not incur allocations.
The C++ coreguiles proposed a similar (albeith one-dimensional) class, called [`gsl::dyn_array`](http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#gslowner-ownership-pointers).

For most uses, a `multi::array` should be preferred instead.

| Member types      | same as for `array_ref` |
|---                |---                        |

| Member fuctions   | same as for `array_ref` plus ... |
|---                |--- |
| (constructors)    | `static_array::static_array({e1, e2, ...}, T val = {}, Alloc = {})` constructs a D-dimensional array by allocating elements. `static_array::static_array(std::initializer_list<...>` constructs the array with elements initialized from a nested list.
| (destructor)      | Destructor deallocates memory and destroy the elements |
| `operator=`       | assigns the elements from the source, sizes must match.

| Element access    | same as for `array_ref` |
|---                |--- |

| Structure access  | same as for `array_ref` |
|---                |--- |

| Iterators         | same as for `array_ref`   |
|---                |--- |

| Capacity          | same as for `array_ref`   |
|---                |--- |

| Creating views    | same as for `array_ref`  |
|---                |---  |

| Creating arrays   | same as for `array_ref`  |
|---                |---  |

| Relational fuctions   |  same as for `array_ref`  |
|---                |--- |

### class `multi::array<T, D, Alloc = std::allocator<T> >`

An array of integer positive dimension D has value semantics if element type T has value semantics.
It supports stateful and polymorphic allocators, which is implied for the special type `multi::pmr::array<T, D>`.

| Member types      | same as for `static_array` |
|---                |---                         |

| Member fuctions   |    |
|---                |--- |
| (constructors)    | `array::array({e1, e2, ...}, T val = {}, Alloc = {})` constructs a D-dimensional array by allocating elements;`array::array(It first, It last)` and `array::array(Range const& rng)`, same for a range of subarrays. `static_array::static_array(std::initializer_list<...>, Alloc = {})` constructs the array with elements initialized from a nested list.
| (destructor)      | Destructor deallocates memory and destroy the elements |
| `operator=`       | assigns for a source `subarray`, or from another `array`. `array`s can be moved |

| Element access    | same as for `static_array` |
|---                |--- |

| Structure access  | same as for `static_array` |
|---                |--- |

| Iterators         | same as for `static_array`   |
|---                |--- |

| Capacity          | same as for `static_array`  |
|---                |--- |

| Creating views    | same as for `static_array`  |
|---                |---  |

| Creating arrays   | same as for `static_array`  |
|---                |---  |

| Relational fuctions   |  same as for `static_array`  |
|---                |--- |

| Manipulation      |     |
|---                |---  |
| `clear`           | Erases all elements from the container. The array is resized to zero size. |
| `reextent`        | Changes the size of the array to new extensions. `reextent({e1, e2, ...})` elements are preserved when possible. New elements are initialized with a default value `v` with a second argument `reextent({e1, e2, ...}, v)`. The first argument is of `extensions_type`, and the second is optional for element types with a default constructor. 

### class `multi::subarray<T, D, P >::(const_)iterator`

A random-access iterator to subarrays of dimension `D - 1`, that is generally used to interact with or implement algorithms.
They can be default constructed but do not expose other constructors since they are generally created from `begin` or `end`, manipulated arithmetically, `operator--`, `operator++` (pre and postfix), or random jumps `operator+`/`operator-` and `operator+=`/`operator-=`.
They can be dereferenced by `operator*` and index access `operator[]`, returning objects of lower dimension `subarray<T, D, ... >::reference` (see above).
Note that this is the same type for all related arrays, for example, `multi::array<T, D, P >::(const_)iterator`.

`iterator` can be invalidated when its original array is invalidated, destroyed or resized.
An `iterator` that stems from `static_array` becomes invalid only if the original array was destroyed or out-of-scope.

# Interoperability with other software

## STL (Standard Template Library)

The fundamental goal of the library is that the arrays and iterators can be used with STL algorithms out-of-the-box with a reasonable efficiency.
The most dramatic example of this is that `std::sort` works with array as it is shown in a previous example.

Along with STL itself, the library tries to interact with other existing quality C++ libraries listed below.

### Ranges (C++20)

[Standard ranges](https://en.cppreference.com/w/cpp/ranges) extend standard algorithms, reducing the need for iterators, in favor of more composability and a less error-prone syntax.

In this example, we replace the values of the first row for which the sum of the elements is odd:

```cpp
	static constexpr auto accumulate = [](auto const& R) {return std::ranges::fold_left(R, 0, std::plus<>{});};

	auto arr = multi::array<int, 2>{
		{2, 0, 2, 2},
		{2, 7, 0, 2},  // this row adds to an odd number
		{2, 2, 0, 4},
	};

	auto const row = std::ranges::find_if(arr, [](auto const& r) { return accumulate(r) % 2 == 1; });
	if(row != arr.end()) std::ranges::fill(*row, 9);

	assert(arr[1][0] == 9 );
```
[(live)](https://godbolt.org/z/cT9WGffM3)

Together with the array constructors, the ranges library enables a more functional programming style;
this allows us to work with immutable variables in many cases.

```cpp
	multi::array<double, 2> const A = {{...}};
	multi::array<double, 1> const V = {...};

	multi::array<double, 1> const R = std::views::zip_transform(std::plus<>{}, A[0], V);

	// Alternative imperative mutating code:
	// multi::array<double, 1> R(V.size());  // R is created here...
	// for(auto i : R.extension()) {R[i] = A[0][i] + V[i];}  // ...and then mutated here
```
[(live)](https://godbolt.org/z/M84arKMnT)


The "pipe" (`|`) notation of standard ranges allows one-line expressions.
In this example, the expression will yield the maximum value of the rows sums:
[`std::ranges::max(arr | std::views::transform(accumulate))`](https://godbolt.org/z/hvqnsf4xb)

Like in classic STL, standard range algorithms acting on sequences operate in the first dimension by default,
for example, lexicographical sorting on rows can be performed with the `std::ranges::sort` algorithm.

```cpp
	auto A = multi::array<char, 2>{
		{'S', 'e', 'a', 'n', ' ', ' '},
		{'A', 'l', 'e', 'x', ' ', ' '},
		{'B', 'j', 'a', 'r', 'n', 'e'},
	};
	assert(!std::ranges::is_sorted(A));

	std::ranges::sort(A);  // will sort on rows

	assert( std::ranges::is_sorted(A));

	assert(
		A == multi::array<char, 2>{
			{'A', 'l', 'e', 'x', ' ', ' '},
			{'B', 'j', 'a', 'r', 'n', 'e'},
			{'S', 'e', 'a', 'n', ' ', ' '},
		}
	);
```

To operate on the second dimension (sort by columns), use `std::ranges::sort(~A)` (or `std::ranges::sort(A.transposed())`).

### Execution policies (parallel algorithms)

Multi's iterators can exploit parallel algorithms by specifying execution policies.
This code takes every row of a two-dimensional array and sums its elements, putting the results in a one-dimensional array of compatible size.
The execution policy (`par`) selected is passed as the first argument.

```cpp
    multi::array<double, 2> const A = ...;
    multi::array<double, 1> v(size(A));

    std::transform(std::execution::par, arr.begin(), arr.end(), vec.begin(), [](auto const& row) {return std::reduce(row.begin(), row.end());} );
```
[(live)](https://godbolt.org/z/63jEdY7zP)

For an array of 10000x10000 elements, the execution time decreases to 0.0288 sec, compared to 0.0526 sec for the non-parallel version (i.e. without the `par` argument).

Note that parallelization is, in this context, inherently one-dimensional.
For example, parallelization happens for the transformation operation, but not to the summation.

The optimal way to parallelize specific operations strongly depends on the array's size and shape.
Generally, straightforward parallelization without exploiting the n-dimensional structure of the data has a limited pay-off;
and nesting parallelization policies usually don't help either.

Flattening the n-dimensional structure for certain algorithms might help, but such techniques are beyond the scope of this documentation.

Some member functions internally perform algorithms and that can benefit from execution policies;
in turn, some of these functions have the option to pass a policy.
For example, this copy construction can initialize elements in parallel from the source:

```cpp
    multi::array<double, 2> const A = ...;
    multi::array<double, 1> const B(std::execution::par, A);  // copies A into B, in parallel, same effect as multi::array<double, 1> const B(A); or ... B = A;
```

Execution policies are not limited to STL;
Thrust and oneAPI also offer execution policies that can be used with the corresponding algorithms.

Execution policies and ranges can be mixed (`x` and `y` can be 1D dimensional arrays, of any arithmetic element type)
```cpp
template <class X1D, class Y1D>
auto dot_product(X1D const& x, Y1D const& y) {
	assert(x.size() == y.size());
	auto const& z = std::ranges::views::zip(x, y)
		| std::ranges::views::transform([](auto const& ab) { auto const [a, b] = ab;
			return a * b;
		})
	;
	return std::reduce(std::execution::par_unseq, z.begin(), z.end());
}
```
[(live)](https://godbolt.org/z/cMq87xPvb)

### Polymorphic Memory Resources

In addition to supporting classic allocators (`std::allocator` by default), the library is compatible with C++17's [polymorphic memory resources (PMR)](https://en.cppreference.com/w/cpp/header/memory_resource), which allows using advanced allocation strategies, including preallocated buffers.
This example code uses a buffer as memory for two arrays; 
in it, a predefined buffer will contain the arrays' data (something like `"aaaabbbbbbXX"`).

```cpp
#include <memory_resource>  // for polymorphic memory resource, monotonic buffer

int main() {
	char buffer[13] = "XXXXXXXXXXXX";  // a small buffer on the stack
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::pmr::array<char, 2> A({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> B({3, 2}, 'b', &pool);

	assert( buffer != std::string{"XXXXXXXXXXXX"} );  // overwritten w/elements, implementation-dependent (libstd consumes from left, and libc++, from the right)
}
```

`multi::pmr::array<T, D>` is a synonym for `multi::array<T, D, std::pmr::polymorphic_allocator<T>>`.
In this particular example, the technique can be used to avoid dynamic memory allocations of small local arrays. [(live)](https://godbolt.org/z/fP9P5Ksvb)

The library also supports memory resources from other libraries, including those returning special pointer types (see the [CUDA Thrust](#cuda-thrust) section and the Boost.Interprocess section).

### Substitutability with standard vector and span

The one-dimensional case `multi::array<T, 1>` is special and overlaps functionality with other dynamic array implementations, such as `std::vector`.
Indeed, both types of containers are similar and usually substitutable, with no or minor modifications.
For example, both can be constructed from a list of elements (`C c = {x0, x2, ...};`) or from a size `C c(size);`, where `C` is either type.

Both values are assignable, have the same element access patterns and iterator interface, and implement all (lexical) comparisons.

They differ conceptually in their resizing operations: `multi::array<T, 1>` doesn't insert or push elements and resizing works differently.
The difference is that the library doesn't implement *amortized* allocations; therefore, these operations would be of a higher complexity cost than the `std::vector`.
For this reason, `resize(new_size)` is replaced with `reextent({new_size})` in `multi::array`, whose primary utility is for element preservation when necessary.

In a departure from standard containers, elements are left initialized if they have trivial constructor.
So, while `multi::array<T, 1> A({N}, T{})` is equivalent to `std::vector<T> V(N, T{})`, `multi::array<T, 1> A(N)` will leave elements `T` uninitialized if the type allows this (e.g. built-ins), unlike `std::vector<T> V(N)` which will initialize the values.
RAII types (e.g. `std::string`) do not have trivial default constructor, therefore they are not affected by this rule.

With the appropriate specification of the memory allocator, `multi::array<T, 1, Alloc>` can refer to special memory not supported by `std::vector`.

Finally, an array `A1D` can be copied by `std::vector<T> v(A1D.begin(), A1D.end());` or `v.assign(A1D.begin(), A1D.end());` or vice versa.
Without copying, a reference to the underlying memory can be created `auto&& R1D = multi::array_ref<double, 1>(v.data(), v.size());` or conversely `std::span<T>(A1D.data_elements(), A1D.num_elements());`. 
(See examples [here](https://godbolt.org/z/n4TY998o4).)

The `std::span` (C++20) has not a well defined reference- or pointer-semantics; it doesn't respect `const` correctness in generic code.
This behavior is contrary to the goals of this library;
and for this reason, there is no single substitute for `std::span` for all cases.
Depending on how it is used, either `multi::array_ref<T, 1> [const& | &&]` or `multi::array_ptr<T [const], 1>` may replace the features of `std::span`.
The former typically works when using it as function argument.

Multi-dimensinal arrays can interoperate with C++23's non-owning `mdspan`.
[Preliminarily](https://godbolt.org/z/aWW3vzfPj), Multi's subarrays (arrays) can be converted (viewed as) `mdspan`.

A detailed comparison with other array libraries (mspan, Boost.MultiArray, Eigen) is explained in an Appendix.

## Serialization

The ability of serializing arrays is important to save/load data to/from disk and also to communicate values via streams or networks (including MPI).
The C++ language does not give any facilities for serialization and unfortunately the standard library doesn't either.

However there are a few libraries that offer a certain common protocol for serialization,
such as [Boost.Serialization](https://www.boost.org/doc/libs/1_76_0/libs/serialization/doc/index.html) and [Cereal](https://uscilab.github.io/cereal/).
The Multi library is compatible with both of them (and yet it doesn't depend on any of them).
The user can choose one or the other, or none, if serialization is not needed.
The generic protocol is such that variables are (de)serialized using the (`>>`)`<<` operator with the archive; operator `&` can be used to have single code for both.
Serialization can be binary (efficient) or text-based (human readable).

Here it is a small implementation of save and load functions for array to JSON format with Cereal.
The example can be easily adapted to other formats or libraries (XML with Boost.Serialization are commented on the right).

```cpp
#include<multi/array.hpp>  // this library

#include<cereal/archives/json.hpp>  // or #include<cereal/archives/xml.hpp>   // #include <boost/archive/xml_iarchive.hpp>
                                                                              // #include <boost/archive/xml_oarchive.hpp>
// for serialization of array elements (in this case strings)
#include<cereal/types/string.hpp>                                             // #include <boost/serialization/string.hpp>

#include<fstream>  // saving to files in example

using input_archive  = cereal::JSONInputArchive ;  // or ::XMLInputArchive ;  // or boost::archive::xml_iarchive;
using output_archive = cereal::JSONOutputArchive;  // or ::XMLOutputArchive;  // or boost::archive::xml_oarchive;

using cereal::make_nvp;                                                       // or boost::serialization::make_nvp;

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
Some formats are human-readable, but not particularly pretty for showing as output (see section on Formatting on how to print to the screen).

References to subarrays (views) can also be serialized; however, size information is not saved in such cases.
The reasoning is that references to subarrays cannot be resized in their number of elements if there is a size mismatch during deserialization.
Therefore, array views should be deserialized as other array views, with matching sizes.

The output JSON file created by Cereal in the previous example looks like this.

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
(The [Cereal XML](https://godbolt.org/z/de814Ycar) and Boost XML output would have a similar structure.)

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

The library works out of the box with Eric Niebler's Range-v3 library, a precursor to the standard Ranges library (see above).
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

Multi doesn't have a dependency on Thrust (or vice versa);
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

In an analogous way, Thrust can also handle OpenMP (omp) allocations and multi-threaded algorithms of arrays.
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

GPU memory is relative expensive to allocate, therefore any application that allocates and deallocates arrays often will suffer performance issues.
This is where special memory management is important, for example for avoiding real allocations when possible by caching and reusing memory blocks.

Thrust implements both polymorphic and non-polymorphic memory resources via `thrust::mr::allocator<T, MemoryResource>`;
Multi supports both.

```cpp
auto pool = thrust::mr::disjoint_unsynchronized_pool_resource(
	thrust::mr::get_global_resource<thrust::universal_memory_resource>(),
	thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
);

// memory is handled by pool, not by the system allocator
multi::array<int, 2, thrust::mr::allocator<int, decltype(pool)>> arr({1000, 1000}, &pool);
```

The associated pointer type for the array data is deduced from the _upstream_ resource; in this case, `thrust::universal_ptr<int>`.

As as quick way to improve performance in many cases, here it is a recipe for a `caching_allocator` which uses a global (one per thread) memory pool that can replace the default Thrust allocator.
The requested memory resides in GPU (managed) memory (`thrust::cuda::universal_memory_resource`) while the cache _bookkeeping_ is held in CPU memory (`new_delete_resource`).

```cpp
template<class T, class Base_ = thrust::mr::allocator<T, thrust::mr::memory_resource<thrust::cuda::universal_pointer<void>>>>
struct caching_allocator : Base_ {
	caching_allocator() : 
		Base_{&thrust::mr::tls_disjoint_pool(
			thrust::mr::get_global_resource<thrust::cuda::universal_memory_resource>(),
			thrust::mr::get_global_resource<thrust::mr::new_delete_resource>()
		)} {}
	caching_allocator(caching_allocator const&) : caching_allocator{} {}  // all caching allocators are equal
	template<class U> struct rebind { using other = caching_allocator<U>; };
};
...
int main() {
	...
	using array2D = multi::array<double, 2, caching_allocator<double>>;

	for(int i = 0; i != 10; ++i) { array2D A({100, 100}); /*... use A ...*/ }
}
```
https://godbolt.org/z/rKG8PhsEh

In the example, most of the frequent memory requests are handled by reutilizing the memory pool avoiding expensive system allocations.
More targeted usage patterns may require locally (non-globally) defined memory resources.

## CUDA C++

CUDA is a dialect of C++ that allows writing pieces of code for GPU execution, known as "CUDA kernels".
CUDA code is generally "low level" (less abstracted) but it can be used in combination with CUDA Thrust or the CUDA runtime library, specially to implement custom algorithms.
Although code inside kernels has certain restrictions, most Multi features can be used. 
(Most functions in Multi, except those involving memory allocations, are marked `__device__` to allow this.)

Calling kernels involves a special syntax (`<<< ... >>>`), and they cannot take arguments by reference (or by values that are not trivial).
Since arrays are usually passed by reference (e.g. `multi::array<double, 2>&` or `Array&&`), a different idiom needs to be used.
(Large arrays are not passed by value to avoid copies, but even if a copy would be fine, kernel arguments cannot allocate memory themselves.)
Iterators (e.g. `.begin()/.end()`) and "cursors" (e.g. `.home()`) are "trivial to copy" and can be passed by value and represent a "proxy" to an array, including allowing the normal index syntax and other transformations.

Cursors are a generalization of iterators for multiple dimensions.
They are cheaply copied (like iterators) and they allow indexing.
Also, they have no associated `.size()` or `.extensions()`, but this is generally fine for kernels.
(Since `cursors` have minimal information for indexing, they can save stack/register space in individual kernels.)

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

Expressions such as `A.begin()` (iterators) can also be passed to kernels, but they could unnecessarely occupy more kernel "stack space" when size information is not needed (e.g. `A.begin()->size()`).

## SYCL

The SYCL library promises the unify CPU, GPU and FPGA code.
At the moment, the array containers can use the Unified Shared Memory (USM) allocator, but no other tests have been investigated.

```cpp
    sycl::queue q;

    sycl::usm_allocator<int, sycl::usm::alloc::shared> q_alloc(q);
    multi::array<int, 1, decltype(q_alloc)> data(N, 1.0, q_alloc);

    //# Offload parallel computation to device
    q.parallel_for(sycl::range<1>(N), [=,ptr = data.base()] (sycl::id<1> i) {
        ptr[i] *= 2;
    }).wait();
```
https://godbolt.org/z/8WG8qaf4s

Algorithms are expected to work with oneAPI execution policies as well (not tested)

```cpp
    auto policy = oneapi::dpl::execution::dpcpp_default;
    sycl::usm_allocator<int, sycl::usm::alloc::shared> alloc(policy.queue());
    multi::array<int, 1, decltype(alloc)> vec(n, alloc);

    std::fill(policy, vec.begin(), vec.end(), 42);
```

## Formatting ({fmt} pretty printing)

The library doesn't have a "pretty" printing facility to display arrays.
Although not ideal, arrays can be printed and formatting by looping over elements and dimensions, as shown as examples in this documentation,
Fortunatelly the library automatically works with the external library [{fmt}](https://fmt.dev/latest/index.html), both for arrays and subarrays.
The fmt library is not a dependency of the Multi library; they simply work well together using the "ranges" part of the formatting library.

This example prints a 2-dimensional subblock of a larger array.

```cpp
#include "fmt/ranges.h"
...
    multi::array<double, 2> A2 = {
        {1.0, 2.0,      3.0}, 
        /*-subblock-**/
        {3.0, 4.0, /**/ 5.0},
        {6.0, 7.0, /**/ 8.0},
    };

    fmt::print("A2 subblock = {}", A2({1, 3}, {0, 2}));  // second and third row, first and second column
```
with the "flat" output `A2 subblock = [[3, 4], [6, 7]]`

For 2 or more dimensions the output can be conveniently structured in different lines using the `fmt::join` facility:

```cpp
    fmt::print("{}\n", fmt::join(A2({1, 3}, {0, 2}), "\n"));  // first dimension rows are printer are in different lines
```
with the output:

> ```
> [3, 4]
> [6, 7]
> ```
https://godbolt.org/z/vjc6n1ove


When saving arrays to files, consider using serialization (see section) instead.

## Legacy libraries (C-APIs)

Multi-dimensional array data structures exist in all languages, whether implicitly defined by its strides structure or explicitly at the language level.
Functions written in C tend to receive arrays by pointer arguments (e.g., to the "first" element) and memory layout (sizes and strides).

A C-function taking a 2D array with a concrete type might look like this in the general case:
```c
void fun(double* data, int size1, int size2, int stride1, int stride2);
```
such a function can be called from C++ on Multi array (`arr`), by extracting the size and layout information,
```cpp
fun(arr.base(), std::get<0>(arr.sizes()), std::get<1>(arr.sizes()), std::get<0>(arr.strides()), std::get<1>(arr.strides());
```
or
```cpp
auto const [size1, size2] = arr.sizes();
auto const [stride1, stride2] = arr.strides();

fun(arr.base(), size1, size2, stride1, stride2);
```

Although the recipe can be applied straightforwardly, different libraries make various assumptions about memory layouts (e.g.,  2D arrays assume that the second stride is 1), and some might take stride information in a different way (e.g., FFTW doesn't use strides but stride products).
Furthermore, some arguments may need to be permuted if the function expects arrays in column-major (Fortran) ordering.

For these reasons, the library is accompanied by a series of adaptor libraries to popular C-based libraries, which can be found in the `include/multi/adaptors/` subdirectory:

- ##### [BLAS/cuBLAS Adator 🔗](include/boost/multi/adaptors/blas/README.md)

Interface for BLAS-like linear algebra libraries, such as openblas, Apple's Accelerate, MKL and hipBLAS/cuBLAS (GPUs).
Simply `#include "multi/adaptors/blas.hpp"` (and link your program with `-lblas` for example).

- ##### Lapack

Interface for Lapack linear solver libraries.
Simply `#include "multi/adaptors/lapack.hpp"` (and link your program with `-llapack` for example).

- ##### FFTW/cuFFT

Interface for FFTW libraries, including FFTW 3, MKL, cuFFT/hipFFT (for GPU).
Simply `#include "multi/adaptors/fftw.hpp"` (and link your program with `-lfftw3` for example).

- ##### [MPI Adaptor 🔗](include/boost/multi/adaptors/mpi/README.md)

Use arrays (and subarrays) as messages for distributed interprocess communication (GPU and CPU) that can be passed to MPI functions through datatypes.
Simply `#include "multi/adaptors/mpi.hpp"`.

- ##### TotalView: visual debugger (commercial)

Popular in HPC environments, can display arrays in human-readable form (for simple types, like `double` or `std::complex`).
Simply `#include "multi/adaptors/totalview.hpp"` and link to the TotalView libraries, compile and run the code with the TotalView debugger.

# Technical points

### What's up with the multiple bracket (vs. parenthesis) notation? 

The chained bracket notation (`A[i][j][k]`) allows to refer to elements and subarrays lower dimensional subarrays in a consistent and _generic_ manner and it is the recommended way to access the array objects.
It is a frequently raised question whether the chained bracket notation is good for performance, since it appears that each utilization of the bracket leads to the creation of a temporary object which in turn generates a partial copy of the layout.
Moreover, this goes against [historical recommendations](https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op).

It turns out that modern compilers with a fair level of optimization (`-O2`) can elide these temporary objects, so that `A[i][j][k]` generates identical machine code as `A.base() + i*stride1 + j*stride2 + k*stride3` (+offsets not shown).
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

Incidentally, the library also supports parenthesis notation with multiple indices `A(i, j, k)` for element or partial access;
it does so as part of a more general syntax to generate sub-blocks.
In any case `A(i, j, k)` is expanded to `A[i][j][k]` internally in the library when `i`, `j`, `k` are normal integer indices.
For this reason, `A(i, j, k)`, `A(i, j)(k)`, `A(i)(j)(k)`, `A[i](j)[k]` are examples of equivalent expressions.

Sub-block notation, when at least one argument is an index ranges, e.g. `A({i0, i1}, j, k)` has no equivalent square-bracket notation.
Note also that `A({i0, i1}, j, k)` is not equivalent to `A({i0, i1})(j, k)`; their resulting sublocks have different dimensionality.

Additionally, array coordinates can be directly stored in tuple-like data structures, allowing this functional syntax:

```cpp
std::array<int, 3> p = {2, 3, 4};
std::apply(A, p) = 234;  // same as assignment A(2, 3, 4) = 234; and same as A[2][3][4] = 234;
```

### Iteration past-end in the abstract machine

It's crucial to grasp that pointers are limited to referencing valid memory in the strict C abstract machine, such as allocated memory.
This understanding is key to avoiding undefined behavior in your code.
Since the library iteration is pointer-based, the iterators replicate these restrictions.

There are three cases to consider; the first two can be illustrated with one-dimensional arrays, and one is intrinsic to multiple dimensions.

The first case is that of strided views (e.g. `A.strided(n)`) whose stride value are not divisors of original array size.
The second case is that or negative strides in general.
The third case is that of iterators of transposed array.

In all these cases, the `.end()` iterator may point to invalid memory. 
It's important to note that the act of constructing certain iterators, even if the element is never dereferenced, is undefined in the abstract machine.
This underscores the need for caution when using such operations in your code.

A thorough description of the cases and workaround is beyond the scope of this section.

# Appendix: Comparison to other array libraries (mdspan, Boost.MultiArray, etc)

The C++23 standard provides `std::mdspan`, a non-owning _multidimensional_ array.
So here is an appropriate point to compare the two libraries.
Although the goals are similar, the two libraries differ in their generality and approach.

The Multi library concentrates on _well-defined value- and reference-semantics of arbitrary memory types with regularly arranged elements_ (distributions described by strides and offsets) and _extreme compatibility with STL algorithms_ (via iterators) and other fundamental libraries.
While `mdspan` concentrates on _arbitrary layouts_ for non-owning memory of a single type (CPU raw pointers).
Due to the priority of arbitrary layouts, the `mdspan` research team didn't find efficient ways to introduce iterators into the library. 
Therefore, its compatibility with the rest of the STL is lacking.
[Preliminarily](https://godbolt.org/z/aWW3vzfPj), Multi array can be converted (viewed as) `mdspan`.

[Boost.MultiArray](https://www.boost.org/doc/libs/1_82_0/libs/multi_array/doc/user.html) is the original multidimensional array library shipped with Boost.
This library can replace Boost.MultiArray in most contexts, it even fulfillis the concepts of `boost::multi_array_concepts::ConstMultiArrayConcept` and `...::MutableMultiArrayConcept`.
Boost.MultiArray has technical and semantic limitations that are overcome in this library, regarding layouts and references;
it doesn't support value-semantics, iterator support is limited and it has other technical problems.

[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) is a very popular matrix linear algebra framework library, and as such, it only handles the special 2D (and 1D) array case.
Instead, the Multi library is dimension-generic and doesn't make any algebraic assumptions for arrays or contained elements (but still can be used to _implement_, or in combination, with dense linear algebra algorithms.)

Other frameworks includes the OpenCV (Open Computing Vision) framework, which is too specialized to make a comparison here.

Here is a table comparing with `mdspan`, R. Garcia's [Boost.MultiArray](https://www.boost.org/doc/libs/1_82_0/libs/multi_array/doc/user.html) and Eigen. 
[(online)](https://godbolt.org/z/555893MqW).


|                             | Multi                                                           | mdspan/mdarray                                                                          | Boost.MultiArray (R. Garcia)                                                                         | Inria's Eigen                                                                           |
|---                          | ---                                                             | ---                                                                             | ---                                                                                                  | ---                                                                                     |
| No external Deps            | **yes** (only Standard Library C++17)                           | **yes** (only Standard Library C++17/C++26)                                                 | **yes** (only Boost)                                                                                 | **yes**                                                                                 |
| Arbritary number of dims    | **yes**, via positive dimension (compile-time) parameter `D`    | **yes**                                                                         | **yes**                                                                                              | no  (only 1D and 2D)                                                                    |
| Non-owning view of data     | **yes**, via `multi::array_ref<T, D>(ptr, {n1, n2, ..., nD})`   | **yes**, via `mdspan m{T*, extents{n1, n2, ..., nD}};`                          | **yes**, via `boost::multi_array_ref<T, D>(T*, boost::extents[n1][n2]...[nD])` | **yes**, via `Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>(ptr, n1, n2)` |
| Compile-time dim size       | no                                                              | **yes**, via template paramaters `mdspan{T*, extent<16, dynamic_extents>{32} }` | no                                                                             | **yes**, via `Eigen::Array<T, N1, N2>` |
| Array values (owning data)  | **yes**, via `multi::array<T, D>({n1, n2, ..., nD})`            | yes for `mdarray`                                               | **yes**, via `boost::multi_array<T, D>(boost::extents[n1][n2]...[nD])` | **yes**, via `Eigen::Array<T>(n1, n2)` |
| Value semantic (Regular)    | **yes**, via cctor, mctor, assign, massign, auto decay of views | yes (?) for `mdarray` planned                                   | partial, assigment on equal extensions  | **yes** (?) |
| Move semantic               | **yes**, via mctor and massign                                  | yes (?) for `mdarray` (depends on adapted container)            | no (C++98 library)                      | **yes** (?) |
| const-propagation semantics | **yes**, via `const` or `const&`                                | no, const mdspan elements are assignable!                       | no, inconsistent                        | (?) |
| Element initialization      | **yes**, via nested init-list                                   | no (?)                                                          | no                                      | no, only delayed init via `A << v1, v2, ...;` |
| References w/no-rebinding   | **yes**, assignment is deep                                     | no, assignment of mdspan rebinds!                               | **yes**                                 | **yes** (?) |
| Element access              | **yes**, via `A(i, j, ...)` or `A[i][j]...`                     | **yes**, via `A[i, j, ...]`                                     | **yes**, via `A[i][j]...`               | **yes**, via `A(i, j)` (2D only) |
| Partial element access      | **yes**, via `A[i]` or `A(i, multi::all)`                       | no, only via `submdspan(A, i, full_extent)`                     | **yes**, via `A[i]`                     | **yes**, via `A.row(i)` |
| Subarray views              | **yes**, via `A({0, 2}, {1, 3})` or `A(1, {1, 3})`              | **yes**, via `submdspan(A, std::tuple{0, 2}, std::tuple{1, 3})` | **yes**, via `A[indices[range(0, 2)][range(1, 3)]]` | **yes**, via `A.block(i, j, di, dj)` |
| Subarray with lower dim     | **yes**, via `A(1, {1, 3})`                                     | **yes**, via `submdspan(A, 1, std::tuple{1, 3})`                | **yes**, via `A[1][indices[range(1, 3)]]`                    | **yes**, via `A(1, Eigen::placeholders::all)` |
| Subarray w/well def layout  | **yes** (strided layout)                                        | no                                                              | **yes** (strided layout)                      | **yes** (strided) |
| Recursive subarray          | **yes** (layout is stack-based and owned by the view)           | **yes** (?)                                                     | no (subarray may dangle layout, design bug?)  | **yes** (?) (1D only) |
| Custom Alloctors            | **yes**, via `multi::array<T, D, Alloc>`                        | yes(?) through `mdarray`'s adapted container                    | **yes** (stateless?)                          | no | 
| PMR Alloctors               | **yes**, via `multi::pmr::array<T, D>`                          | yes(?) through `mdarray`'s adapted container                    |   no                          | no |
| Fancy pointers / references | **yes**, via `multi::array<T, D, FancyAlloc>` or views          | no                                                              |   no                          | no |
| Stride-based Layout         | **yes**                                                         | **yes**                                                   |  **yes**                      | **yes** |
| Fortran-ordering            | **yes**, only for views, e.g. resulted from transposed views    | **yes**                                                   |  **yes**.                     | **yes** |
| Zig-zag / Hilbert ordering  | no                                                              | **yes**, via arbitrary layouts (no inverse or flattening) | no                            | no |
| Arbitrary layout            | no                                                              | **yes**, possibly inneficient, no efficient slicing       | no                            | no |
| Flattening of elements      | **yes**, via `A.elements()` range (efficient representation)    | **yes**, but via indices roundtrip (inefficient)          | no, only for allocated arrays | no, not for subblocks (?) |
| Iterators                   | **yes**, standard compliant, random-access-iterator             | no                                                        | **yes**, limited | no |
| Multidimensional iterators (cursors) | **yes** (experimental)                                 | no                                                        | no               | no |         
| STL algorithms or Ranges    | **yes**                                                         | no, limited via `std::cartesian_product`                  | **yes**, some do not work | no |
| Compatibility with Boost    | **yes**, serialization, interprocess  (see below)               | no                                                        | no | no |
| Compatibility with Thrust or GPUs | **yes**, via flatten views (loop fusion), thrust-pointers/-refs | no                                                  | no          | no |
| Used in production          | [QMCPACK](https://qmcpack.org/), [INQ](https://gitlab.com/npneq/inq)  | (?) , experience from Kokkos incarnation            | **yes** (?) | [**yes**](https://eigen.tuxfamily.org/index.php?title=Main_Page#Projects_using_Eigen) |

# Appendix: Multi for FORTRAN programmers

This section summarizes simple cases translated from FORTRAN syntax to C++ using the library.
The library strives to give a familiar feeling to those who use arrays in FORTRAN, including the manipulation of arrays with arbitrary dimensionality (e.g. 1D, 2D, 3D, etc.)
Arrays can be indexes arrays using square-brakers or using parenthesis, which would be more familiar to FORTRAN syntax.
The most significant difference is that in FORTRAN, array indices start at `1` by default, while in Multi, they start at `0` by default, following C++ conventions.
Like in FORTRAN, for simple types (e.g., numeric), Multi arrays are not initialized automatically; such initialization needs to be explicit.

|                             | FORTRAN                                          | C++ Multi                                            |
|---                          | ---                                              | ---                                                  |
| Declaration/Construction 1D | `real, dimension(2) :: numbers` (at top)         | `multi::array<double, 1> numbers(2);` (at scope)     |
| Initialization (2 elements) | `real, dimension(2) :: numbers = [ 1.0, 2.0 ]`   | `multi::array<double, 1> numbers = { 1.0, 2.0 };`    |
| Element assignment          | `numbers(2) = 99.0`                              | `numbers(1) = 99.0;` (or `numbers[1]`)               |
| Element access (print 2nd)  | `Print *, numbers(2)`                            | `std::cout << numbers(1) << '\n';`                   |
| Initialization              | `DATA numbers / 10.0 20.0 /`                     | `numbers = {10.0, 20.0};`                            |

In the more general case for the dimensionality, we have the following correspondance:

|                              | FORTRAN                                          | C++ Multi                                            |
|---                           | ---                                              | ---                                                  |
| Construction 2D (3 by 3)     | `real*8 :: A2D(3,3)` (at top)                    | `multi::array<double, 2> A2D({3, 3});` (at scope)    |
| Construction 2D (2 by 2)     | `real*8 :: B2D(2,2)` (at top)                    | `multi::array<double, 2> B2D({2, 2});` (at scope)    |
| Construction 1D (3 elements) | `real*8 :: v1D(3)`   (at top)                    | `multi::array<double, 2> v1D({3});` (at scope)       |
| Assign the 1st column of A2D | `v1D(:) = A2D(:,1)`                              | `v1( _ ) = A2D( _ , 0 );`                            | 
| Assign the 1st row of A2D    | `v1D(:) = A2D(1,:)`                              | `v1( _ ) = A2D( 0 , _ );`                            |
| Assign upper part of A2D     | `B2D(:,:) = A2D(1:2,1:2)`                        | `B2D( _ , _ ) = A2D({0, 2}, {0, 2});`                |

Note that these correspondences are notationally logical.
Internal representation (memory ordering) can still be different, and this could affect operations that interpret 2D arrays as contiguous elements in memory.

Range notation such as `1:2` is replaced by `{0, 2}`, which takes into account both the difference in the start index and the half-open interval notation in the C++ conventions.
Stride notation such as `1:10:2` (i.e. from first to tenth included, every 2 elements), is replaced by `{0, 10, 2}`.
Complete range interval (single `:` notation) is replaced by `multi::_`, which can be used simply as `_` after the declaration `using multi::_;`.
These rules extend to higher dimesionality, such as 3 dimensions.

Unlike FORTRAN, Multi doesn't provide algebraic operators, using algorithms is encouraged instead.
For example a FORTRAN statement like `A = A + B` is translated as this in the one-dimensional case:

```cpp
std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::plus{});  // valid for 1D arrays only
```

In the general dimensionality case we can write:

```cpp
auto&&      Aelems = A.elements();
auto const& Belems = B.elements();
std::transform(Aelems.begin(), A.elems.end(), Belems.begin(), Aelems.begin(), std::plus{});  // valid for arbitrary dimension
```

or
```
std::ranges::transform(Aelems, Belems, Aelems.begin(), std::plus{});  // alternative using C++20 ranges
```

A FORTRAN statement like `C = 2.0*C` is rewritten as `std::ranges::transform(C.elements(), C.elements().begin(), [](auto const& e) {return 2.0*e;});`.

It is possible to use C++ operator overloading for functions such as `operartor+=` (`A += B;`) or `operator*=` (`C *= 2.0;`); however, this possibility can become unwindenly complicated beyond simple cases.
Also it can become inefficient if implemented naively.

Simple loops can be mapped as well, taking into account indexing differences:
```fortran
do i = 1, 5         ! for(int i = 0; i != 5; ++i) {
  do j = 1, 5       !   for(int j = 0; j != 5; ++j) {
    D2D(i, j) = 0   !     D2D(i, j) = 0;
  end do            !   }
end do              ! }
```

However, algorithms like `transform`, `reduce`, `transform_reduce`and `for_each`, and offer a higher degree of control over operations, including memory allocations if needed, and even enable parallelization, providing a higher level of flexibility.


[(live)](https://godbolt.org/z/77onne46W)


> Thanks to Joaquín López Muñoz and Andrzej Krzemienski for the critical reading of the documentation and to Matt Borland for his help integrating Boost practices in the testing code.
