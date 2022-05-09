<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->
# [Boost.]Multi

(not an official Boost library)

_© Alfredo A. Correa, 2018-2022_

`Multi` provides multidimensional array access to contiguous or regularly contiguous memory (or ranges).
It shares the goals of [Boost.MultiArray](https://www.boost.org/doc/libs/1_69_0/libs/multi_array/doc/index.html),
although the code is completely independent and the syntax has slight differences or has been extended.
`Multi` and `Boost.MultiArray` types can be used interchangeably for the most part, they differ in the semantics of reference and value types.

Multi aims to simplify the semantics of Boost.MultiArray and make it more compatible with the Standard (STL) Algorithms and special memory.
It requires, at least, C++17. (It is C++20 ready.)

Some features:

* Arbitrary pointer types (minimal requirements)
* Simplified implementation (~4000 lines)
* Fast access of subarrays (view) types
* Value semantics of multi-dimensional array container
* Better semantics of subarray (view) types
* Interoperability with other libraries, STL, ranges, 

(Do not confuse this library with Boost.MultiArray or Boost.MultiIndex.)

## Contents
[[_TOC_]]

## Installation and Tests

`Multi` doesn't require instalation, single file `#include<multi/array.hpp>` is enough to use the full core library.
`Multi`'s _only_ dependecy is the standard C++ library.

It is important to compile programs that use the library with a decent level of optimization, specially if element-access is intensively used.
For example, when testing speed, please make sure that you are compiling in release mode (`-DNDEBUG`) and with optimizations (`-O3`),
if your test involves mathematical operations add arithmetic optimizations (`-Ofast`) to compare with Fortran code.

A CMake build system is provided to automatically run basic tests. (Test do depend on the Boost.Test library.)

```bash
git clone https://gitlab.com/correaa/boost-multi.git multi
cd multi
mkdir -p build
cd build
cmake ..
make -j
make test -j
```

### Dependecies and compiler requirements

The core of the library doesn't have dependencies (other than the standard library).

Compiling and running the tests depends on Boost.Test
(which can be installed with `sudo apt install libboost-test-dev` in Debian-like systems.)

"Adaptor" sublibraries (included in `multi/adaptors/`) have specific dependencies, Boost.Serialization, fftw, blas, lapack, thurst, CUDA
(which can be installed with `sudo apt install libboost-serialization-dev libfftw3-dev libblas64-dev liblapack64-dev libthrust-dev libcudart11.0` or indiviudually.)

The code is developed for several compilers with standard C++17 support, for reference:
LLVM's `clang` (9+), GNU's `g++` (7+), Nvidia's `nvcc` (11.3+) and `nvc++` (20.7-21.3+), Intel's `icpc` (2021.2.0+) and `icpx` (2022.0.2+) and Baxter's [`circle`](https://www.circle-lang.org/) (build 168). 
For detailed compilation instructions of test please inspect the Continuous Integration (CI) [definition file](https://gitlab.com/correaa/boost-multi/-/blob/master/.gitlab-ci.yml).

## Types

* `multi::array<T, D, A>`: 
Array of dimension `D`, it has value semantics if `T` has value semantics. 
Memory is requested by allocator of type `A`, supports stateful and polymorphic allocators.
* `multi::array_ref<T, D, P = T*>`: 
Array interpretation of a random access range, usually a contiguous memory block. 
It has reference semantics. 
Thanks to (non-virtual) inheritance an `array<T, D, A>` is-a `array_ref<T, D, A::pointer>`.
* Other derived "unspecified types" fulfil (a still loosely defined) `MultiArrayView` concept, for example by taking partial indices or rotations (transpositions).
These reference types cannot be stored except through life-time extensions `auto&&` or `auto const&`,
and the can decay to value types.
* `MultiArrayView<T,D,P>::(const_)iterator`:
Iterator to subarrays of dimension `D - 1`. For `D == 1` this is an iterator to an element. This types are generated by `begin` and `end` functions.
* `MultiArrayView<T, D, P>::(const_)reference`: 
Reference to subarrays of dimension `D - 1`. For `D > 1` this are not true C++-references but types emulate them (with reference semantics), therefore `auto` is not well behaved.
For `D == 1` this is a true C++ reference to an elements. These types are generated by dereferencing iterators, e.g. `*begin(MA)`.

## Basic Usage

Declare an array by specifying the element type and the dimension;
indiviudual elements can be input with nested braced notation.
```cpp
multi::array<double, 2> A = {
	{1., 2., 3.}
	{4., 5., 6.}
};

assert( A.size() == 2 );
assert( A.num_elements() == 6 );

assert( std::get<0>(A.sizes()) == 3 );
assert( std::get<1>(A.sizes()) == 3 );
```

(The size is automatically deduced; the first (leading) dimension are the (two) "rows" above.)

The value of an array can be copied, (moved,) and compared;
copies are equal but independent.

```cpp
std::array<double, 2> B = A;
assert( extensions(B) == extensions(A) );
assert(  B       ==  A                 );
assert(  B[0][1] ==  A[0][1]           );
assert( &B[0][1] != &A[0][1]           );
```

Individual elements can be accessed by the multidimensional indices, either with square bracket (one index at a time) or with parenthesis.

```
assert(  A(1, 2)  ==  A[1][2] );
```

Arrays can be initialized from its sizes alone, in which case the element values are default constructed:

```cpp
multi::array<double, 3> C({3, 4, 5});
assert( num_elements(C) == 3*4*5 );  // 60 elements
```

Arrays can be passed by value or by reference.
Most of the time, arguments should be passed through generic parameters to allow functions to work with parts (subblocks, slices, etc.) of an array.
Most useful function work on the concept of array rather than on a concrete type.

```cpp
template<class ArrayDouble2D>  // instead of the over specific argument std::array<double, 2>
auto element_1_1(ArrayDouble2D const& m) -> double const& {return m[1][1];}
...
assert( element_1_1(A) == A[1][1] );
```

The generic function arguments that are not intended to be modified are passed by `const&`; otherwise pass by forward-reference `&&`.
In this way the functions can be called on subblocks of larger matrices.

```cpp
assert( &element_1_1(C3D[0]) == &C3D[0][1][1] );
```

## Advanced Usage

We create a static C-array of `double`s, and refer to it via a bidimensional array `multi::array_ref<double, 2>`.

```cpp
#include "../array.hpp"

#include<algorithm>  // for sort
#include<iostream>  // for print

namespace multi = boost::multi;
using std::cout;

int main() {
	double d_data[20] = {
		150., 16., 17., 18., 19.,
		 30.,  1.,  2.,  3.,  4.,
		100., 11., 12., 13., 14.,
		 50.,  6.,  7.,  8.,  9.
	};  // block of 20 elements ...
	multi::array_ref<double, 2> d2D_ref{&d_data[0], {4, 5}};  // interpreted as a 4 by 5 array
	...
```

Note that the syntax of creating a reference array involves passing the pointer to a memory block (20 elements here) and the logical dimensions of that memory block (4 by 5 here).

Next we print the elements in a way that corresponds to the logical arrangement:

```cpp
	...
	auto [ix, jx] = d2D_ref.extensions();
	for(auto i : is) {
		for(auto j : js) {
			cout<< d2D_ref[i][j] <<' ';
		}
		cout <<'\n';
	}
	...
```

This will output:

> ```cpp
> 150 16 17 18 19  
> 30 1 2 3 4  
> 100 11 12 13 14  
> 50 6 7 8 9
> ```

It is sometimes said (by Sean Parent) that the whole of STL algorithms can be seen as intermediate pieces to implement`std::stable_sort`. 
Pressumably if one can sort over a range, one can perform any other standard algorithm.

```cpp
		...
		std::stable_sort( begin(d2D_ref), end(d2D_ref) );
		...
```

If we print this we will get

> ```cpp
> 30 1 2 3 4  
> 50 6 7 8 9  
> 100 11 12 13 14  
> 150 16 17 18 19
> ```


The array has been changed to be in row-based lexicographical order.
Since the sorted array is a reference to the original data, the original array has changed.

(Note that `std::*sort` cannot be applied directly to a multidimensional C-array or to Boost.MultiArray types.
The arrays implemented by this library are, to the best of my knowledge, the only ones that support STL algorithms directly.)

If we want to order the matrix in a per-column basis we need to "view" the matrix as range of columns.
This is done in the bidimensional case, by accessing the matrix as a range of columns:

```cpp
		...
		std::stable_sort( rotated(d2D_ref).begin(), rotated(d2D_ref).end() );
	}
```

Which will transform the (original) matrix into:

> ```cpp
> 1 2 3 4 30  
> 6 7 8 9 50  
> 11 12 13 14 100  
> 16 17 18 19 150 
> ```

In other words, a matrix of dimension `D` can be viewed simultaneously as `D` different ranges of different "transpositions".

## Initialization

`array_ref` is initialized from a preexisting contiguous range, the index extensions should compatible with the total number of elements.

```cpp
double* dp = new double[12];
multi::array_ref<double, 2> A({3, 4}, dp);
multi::array_ref<double, 2> B({2, 6}, dp);
...
delete[] dp;
```
Array references do not own memory and can not be resized or "reseated" to refer to a different location.
Since `array_ref` is an array reference, it can "dangle" if the the original memory is deallocated.

Array objects own the elements it contains and can be resized;
`array` is initialized by specifying the index extensions (and optionally a default value).

```cpp
multi::array<double, 1> A1({3}      , 11.);  // {11., 11., 11.}

multi::array<double, 2> A2({2, 3}   , 22.);  // { {22., 22., 22.}, {22., 22., 22.} }

multi::array<double, 3> A3({3, 2, 2}, 33.);  // { { { 33., ...}, { ... }, ... } }
```
... or alternatively from a rectangular list.

```cpp
multi::array<double, 1> A1 = {1., 2., 3.};
assert( num_elements(A1)==3 );

multi::array<double, 2> A2 {
	 { 1., 2., 3.},
	 { 4., 5., 6.}
};

assert( num_elements(A2) == 2*3);

multi::array<double, 3> const A3 = {
    {{ 1.2,  0.}, { 2.4, 1.}},
    {{11.2,  3.}, {34.4, 4.}},
    {{15.2, 99.}, {32.4, 2.}}
};

assert( num_elements(A3) == 3*2*2 );
```

In all cases constness (`const` declaration) is honored in the expected way.

## Changing extents (sizes)

Arrays can change their size preserving elements with `reextents`.

```cpp
multi::array<double, 2> A {
	 {1., 2., 3.},
	 {4., 5., 6.}
};

A.reextents({4, 4});

assert( A[0][0] = 1. );
```

Subarrays, or views cannot change their size. `A[1].reextents({4})`.
The main utility of `reextents` is element preservation.
If element preservation is not desired, simple (move) assignment from a new array expresses the intention better as it is more efficient.

```cpp
A = multi::array<double, 2>({4, 4});
```

## Iteration

Accessing arrays by iterators (`begin`/`end`) enables the use of many iterator based algorithms (see the sort example above).
`begin(A)/end(A)` (or equivalently `A.begin()/A.end()`) gives iterators that are linear and random access in the leading dimension.

`cbegin/cend` give constant (read-only access).

Other non-leading dimensions can be obtained by "rotating" indices first.

`A.rotated().begin()/.end()` gives access to a range of subarrays in second dimension number (first dimension is put at the end).

For example in a three dimensional array,

```cpp
	(cbegin(A)+1)->operator[](1).begin()[0] = 342.4;  // error, read-only
	( begin(A)+1)->operator[](1).begin()[0] = 342.4;  // assigns to A[1][1][0]
	assert( ( begin(A)+1)->operator[](1).begin()[0] == 342.4 );
```

As an example, this function allows printing arrays of arbitrary dimension into a linear comma-separated form.

```cpp
void print(double const& d){cout<<d;};
template<class MultiArray>
void print(MultiArray const& ma) {  // note the recursion in the template function `print`
	cout<<"{";
	if(not ma.empty()) {
		print(*cbegin(ma));
		std::for_each(cbegin(ma)+1, cend(ma), [](auto&& e) {cout<<","; print(e);});
	}
	cout<<"}";
}
...
print(A);
```
> ```
> {{{1.2,1.1},{2.4,1}},{{11.2,3},{34.4,4}},{{15.2,99},{32.4,2}}}
> ```


Except for those corresponding to the one-dimensional case, derreferencing iterators generally produce "proxy"-references (i.e. objects that behave in a large degree like language references).
These references can be given a name; using `auto` can be misleading since the resulting variable does not have value semantics.

```cpp
auto row = *begin(A);  // accepted by the language but misleading, row is not an independent value
```

In my experience, however, the following usage patter a more consistent idiom for generating references (still without copying elements):

```cpp
auto&&       row0 = * begin(A);  // same as decltype(A)::      reference  row0 = * begin(A);
auto const& crow0 = *cbegin(A);  // same as decltype(A)::const_reference crow0 = *cbegin(A);

auto&&       row0 =               A [1];  // same as decltype(A)::      reference  row0 =               A [1];
auto const& crow1 = std::as_const(A)[1];  // same as decltype(A)::const_reference crow0 = std::as_const(A)[1];
```


If a new value is desired, these (equivalent) options express the intention more explicitly:

```cpp
decltype(A)::value_type row =   *begin(A);  // there is a real copy of the row
                   auto row = + *begin(A);  // there is another copy, note the use of '+' (unary plus)
```

## Indexing

Arrays provide random access to elements or subviews.
Many algorithms on arrays are oriented to linear algebra,
which are ubiquitously implemented in terms of multidimensional index access.

Iterator access and index access are two alternatives for accessing elements.
For example `*(begin(A) + n)` and `A[n]` are semantically equivalent
and the range defined by the pair `begin(A), end(A)` is `A(extension(A))` (even for multidimensional `A`).
The syntax can be combined in arbitrary ways, for example `*begin(A[n])` is equivalent to `A[n][0]` (if the dimensionality of `A` is equal or greater than two).

### Element access and partial access

Index access mimics that of C-fixed sizes arrays. 
For example, a 2-dimensional array will access to an element by specifying two indices `A[1][2]`,
which can be used for direct write and read operations; 
while _partial_ index arguments `A[1][2]` generate a view 1-dimensional object (reference).

```cpp
A        // is a 2D value array
A[0]     // is a 1D "reference"/"view" array
A[0][0]  // is a an element reference, zero-D
```

Transpositions are also multi-dimensional arrays views in which the index are *logically* rearranged, for example `rotated(m)[2][3][1] == m[1][2][3]`.
(`rotate` refers to the fact that the logical indices are _rotated_ to the left.)

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
multi::array<double, 2> A = {{-3., 2., -4.},{0., 1., 2.},{2., 4., 5.}};
multi::array<double, 1> y = {12.,5.,2.}; //(M); assert(y.size() == M); iota(y.begin(), y.end(), 3.1);
gj_solve(A, y);
```

and also to a combination of `MultiArrayView`-type objects (including standard vectors):

```cpp
multi::array<double, 2> A({6000, 7000}); std::iota(A.data_elements(), A.data_elements() + A.num_elements(), 0.1);
std::vector<double> y(3000); std::iota(y.begin(), y.end(), 0.2);  // could be also a multi::array<double, 1> y({3000});
gj_solve(A({1000, 4000}, {0, 3000}), y);
```

### Slices and strides

Given an array, a slice in the first dimension can be taken with the `sliced` function. `sliced` takes two arguments, the first index of the slice and the last index (not included) of the slice. For example,

```cpp
multi::array<double, 2> A({4, 5});  // A is a value
assert( std::get<0>(A) == 2 );
assert( std::get<1>(A) == 5 );

auto&& A_sliced = A.sliced(1, 3); // {{d2D[1], d2D[2]}}
assert( std::get<0>(A_sliced) == 2 );
assert( std::get<1>(A_sliced) == 5 );
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
assert( d2D_slicedstrided.size(0) == 1 and d2D_slicedstrided.size(1) == 5 );
```

For convenience, `A.sliced(a, b, c)` is the same as `A.sliced(a, b).strided(c)`.

By combining `rotated`, `sliced` and `strided` one can take sub arrays at any dimension.
For example in a two dimensional array one can take a subset of columns by defining.

```cpp
auto&& subA = A.rotated().strided(1, 3).sliced(2).unrotated();
```

Other notations are available, but when in doubt, the `rotated/strided/sliced/rotated` and combinations of them provides the most control over the subview operations.
(At the moment the `strided` argument has to divide the total size of the slice (or matrix), otherwise the behavior is undefined.)

Blocks (slices) in multidimensions can be obtained but pure index notation using parentheses `()` (`.operator()`):

```cpp
auto        A = multi::array<double, 2>({6, 7});  // 2D value array

auto&&      A_block1 = A({1, 4}, {2, 4});  // 2D subarray reference (modifiable)
auto const& A_block2 = A({1, 4}, {2, 4});  // 2D subarray reference (non-modifiable)

auto        A_block3 = A({1, 4}, {2, 4});  // disabled
```

Note that the last case gives a compilation error, the library prevents the use of this references as if they are were values.
Some times copies are necessary, specifically from a subarray block, this can be done by constructing a new array. 
The value array can be deduces by using `auto` and the `decay` member, which in turn is equivalent to the unary `+` operator.

```cpp
multi::array<double, 2> block_value_1 =   A({1, 4}, {2, 4})        ;
auto                    block_value_2 =   A({1, 4}, {2, 4}).decay();
auto                    block_value_3 = + A({1, 4}, {2, 4})        ;
```


## Concept Requirements

The design tries to impose the minimum possible requirements over the  types that parameterize the arrays.
Array operations assume that the contained type (element type) are regular (i.e. different element represent disjoint entities that behave like values).
Pointer-like random access types can be used as substitutes of built-in pointers.
Therefore pointers to special memory (fancy-pointers) are supported.

```cpp
int main() {
	double* buffer = new double[100];
	multi::array_ref<double, 2, minimal::ptr<double> > CC(minimal::ptr<double>{buffer}, {10, 10});
	CC[2]; // requires operator+ 
	CC[1][1]; // requires operator*
	CC[1][1] = 9;
	assert(CC[1][1] == 9);
	delete[] buffer;
}
```

### Linear Sequences: Pointers

An `array_ref` can reference to an arbitrary random access sequence (e.g. memory block defined by pointer and size).
This way, any linear (random access) sequence (e.g. `raw memory`, `std::vector`, `std::queue`) can be efficiently arranged as a multidimensional array.

```cpp
std::vector<double> buffer(100);
multi::array_ref<double, 2> A({10, 10}, buffer.data());
A[1][1] = 9;

assert(buffer[11]==9);  // the target memory is affected
```
Since `array_ref` does not manage the memory associated with it, the reference can be simply dangle if the `buffer` memory is reallocated (e.g. by `resize` in this case).

### Special Memory: Allocators and Fancy Pointers

`array`'s manages its memory behind the scenes through allocators, which can be specified at construction.
It can handle special memory, as long as the underlying types behave coherently, these include fancy pointers and fancy references.
Associated fancy pointers and fancy reference (if any) are deduced from the allocator types.

The behavior regarding memory managament of the [fancy pointers](https://en.cppreference.com/w/cpp/named_req/Allocator#Fancy_pointers) can be customized (if necessary) by specializations of some or all of these functions:

```cpp
destroy(alloc, first, last);
destroy_n(alloc, first, n) -> last
uninitialized_copy_n(alloc, first, n, dest) -> last;
uninitialized_fill_n(alloc, first, n, value) -> last;
uninitialized_default_construct_n(alloc, first, n) -> last;
uninitialized_value_construct_n(alloc, first, n) -> last;
```

where `alloc` is the special allocator, `n` is a size (usually the number of elements), `first`, `last` and `dest` are fancy pointers.

Copying underlying memory can be customized by specializing:

```cpp
copy_n(first, n, dest)
fill_n(first, n, value)
```

Specific cases of fancy memory are file-mapped memory or interprocess shared memory.
This example illustrates memory persistency by combining with Boost.Interprocess library. 
The arrays support their allocators and fancy pointers (`boost::interprocess::offset_ptr`).

```cpp
#include <boost/interprocess/managed_mapped_file.hpp>
using namespace boost::interprocess;
using manager = managed_mapped_file;
template<class T> using mallocator = allocator<T, manager::segment_manager>;
decltype(auto) get_allocator(manager& m){return m.get_segment_manager();}

template<class T, auto D> using marray = multi::array<T, D, mallocator<T>>;

int main(){
{
	manager m{create_only, "mapped_file.bin", 1 << 25};
	auto&& arr2d = *m.construct<marray<double, 2>>("arr2d")(std::tuple{1000, 1000}, 0.0, get_allocator(m));
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

# Interoperability with other software

## STL (Standard Template Library)

The fundamental goal of the library is that the arrays and iterators can be used with STL algorithms out-of-the-box with a reasonable efficiency.
The most dramatic example of this is that `std::sort` works with array as it is shown in a previous example.

Along with STL itself, the library tries to interact with other existing quality C++ libraries described below.

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
#include <multi/array.hpp>  // our library
#include<fstream>  // saving to files in example
#include <cereal/archives/json.hpp>                // #include <boost/archive/xml_iarchive.hpp>
                                                   // #include <boost/archive/xml_oarchive.hpp>
// for serialization of array elements (in this case strings)
#include <cereal/types/string.hpp>                 // #include <boost/serialization/string.hpp>
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

These templated functions work for any dimension and element type (as long as the element type is serializable in itself; all basic types are serializable by default).
However note that it is responsibility of the user to make sure that data is serialized and deserialized into the same type and also assuming the same format.
This is because the underlying serialization library only do minimal consistency checks for efficiency reasons and doesn't try to second guess file formats or contained types.
Serialization is a relatively low level feature for which efficiency and economy of bytes is priority.
Cryptic errors and crashes can occur if serialization libraries, file formats or C++ types are mixed between writes and reads.
On top of serialization checks can be added by the user before and after loading a file.

References to subarrays can also be serialized, however, in such case size information is not saved.
The reason is that references to subarrays cannot be resized in their number of elements if there is size mismatch during deserialization.

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

## Polymorphic Memory Resources

The library is compatible with C++17's polymorphic memory resources (PMR) which allows using preallocated buffers as described in this example. 
This enables the use of stack memory, with many performance advantaneges.
For example, this code uses a buffer to allocate memory for two arrays, we will see how this buffer ends up containing the data of the arrays `"aaaabbbbbbXX"`.

```cpp
#include <memory_resource>  // polymorphic memory resource, monotonic buffer, needs C++17

int main() {
	char buffer[13] = "XXXXXXXXXXXX";  // a small buffer on the stack
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> A({2, 2}, 'a', &pool);
	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> B({3, 2}, 'b', &pool);

	assert( buffer == std::string{"aaaabbbbbbXX"} );
}
```

The library supports classic allocators (`std::allocator` by default) and also allocators from other libraries (see Thurst section).

## Range v3

```cpp
#include <range/v3/all.hpp>
int main(){

	multi::array const d2D = {
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
	m.destroy<marray<double, 2>>("arr2d");//	eliminate<marray<double, 2>>(m, "arr2d");}
}
}
```

(Similarly works with [LLNL's Meta Allocator](https://github.com/llnl/metall))

## Cuda thrust

```cpp
#include "multi/adaptors/thrust/allocator_traits.hpp"
#include "multi/adaptors/thrust/algorithms.hpp"
#include "multi/array.hpp"

namespace multi = boost::multi;
int main(){
	multi::array<double, 2, thrust::device_allocator<double>> A2({10,10});
	multi::array<double, 2, thrust::device_allocator<double>> B2({10,10});
	A2[5][0] = 50.;
	thrust::copy(begin(rotated(A2)[0]), end(rotated(A2)[0]), begin(rotated(B2)[0]));
	assert( B2[5][0] == 50. );
}
```

## TotalView

TotalView visual debugger (commercial) can display arrays in human-readable form (for simple types, like `double` or `std::complex`).
To use it, simply `#include "multi/adaptors/totalview.hpp"` and link to the TotalView libraries, compile and run the code with the TotalView debugger.

# Technical points

### What's up with the multiple bracket notation? 

The chained bracket notation (`A[i][j][k]`) allows to refer to elements and subarrays lower dimensional subarrays in a consistent and _generic_ manner and it is the recommended way to access the array objects.
It is a frequently raised question whether the chained bracket notation is good for performance, since it appears that each utilization of the bracket leads to the creation of a temporary which in turn generates a partial copy of the layout.
Moreover, this goes against [historical recommendations](https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op).

It turns out that [modern compilers with a fair level of optimization (`-O2`)](https://godbolt.org/z/3fYd5c) can elide these temporary objects, so that `A[i][j][k]` generates identical assembly code as `A.base() + i*stride1 + j*stride2 + k*stride3` (+offsets not shown).

In a subsequent optimization, constant indices can have their "partial stride" computation removed from loops. 
As a result, these two loops lead to the [same machine code](https://godbolt.org/z/z1se74):

```cpp
    for(int j = 0; j != nj; ++j)
        ++A[i][j][k];
```
```cpp
    double* Ai_k = A.base() + i*A_stride1 + k*A_stride3;
    for(int j = 0; j != nj; ++jj)
        ++(*(Ai_k + j*A_stride2));
```

Incidentally, the library also supports parenthesis notation with multiple indices `A(i, j, k)` for element or partial access, but it does so for accidental reasons as part of a more general syntax to generate sub-blocks.
In any case `A(i, j, k)` is expanded to `A[i][j][k]` internally in the library when `i, j, k` are normal integer indices.
Additionally, array coordinates can be directly stored in tuple-like data structures, allowing this functional syntax:

```cpp
std::array p = {2,3,4};
std::apply(A, p) = 234; // A[2][3][4] = 234;
```

### Customizing recursive operations: SCARY iterators

A custom level of customization can be achieved by intercepting internal recursive algorithms.
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

