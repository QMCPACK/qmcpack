<!--
(pandoc `#--from gfm` --to html --standalone --metadata title=" " $0 > $0.html) && firefox --new-window $0.html; sleep 5; rm $0.html; exit
-->

<!--
> **⚠️ ALERT**
> This library is under active Boost review until **March 15, 2026**.
> If you are interested in reviewing the library, please send an email to the review manager, **Matt Borland** (matt AT mattborland DOT com).
-->

# Boost.Multi

Boost.Multi provides owning multidimensional arrays containers and non-owning multidimensional views for C++17. It supports slicing, layout transformations, and CPU/GPU memory.

> **Project status:** Boost.Multi is not an official or accepted Boost library. It is being proposed for inclusion in [Boost](https://www.boost.org/) and has no Boost-library dependencies.

_© Alfredo A. Correa, 2018–2026_

## Why Boost.Multi?

* **Generic element type and dimensionality**
* **Natural multidimensional access:** write `A[i][j]...` or iterate recursively over rows, subarrays or elements.
* **Composable subarray views:** slice, transpo
se, rotate indices, and broadcast arrays without copying their elements.
* **Interoperability:** use standard algorithms (iterators and ranges), allocators, and legacy libraries for strided arrays (BLAS, FFTW, CUDA, HIP, MPI, etc.)

## Quick start

```bash
git clone https://github.com/correaa/boost-multi.git  # or gitlab.com/correaa/boost-multi.git
cd boost-multi
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Read the [documentation](https://correaa.github.io/boost-multi/multi/intro.html) or explore `test/`.

## First array and view

```cpp
#include <cassert>
#include <multi/array.hpp>

namespace multi = boost::multi;

int main() {
    multi::array<int, 2> A = {
        {1, 2, 3},
        {4, 5, 6}
    };

    assert(A[1][1] == 5);  // ordinary multidimensional indexing

    using multi::_;
    auto&& second_column = A(_, 1);  // a non-owning view
    second_column[0] = 20;           // changes A[0][1]
    assert(A[0][1] == 20);
}
```

## Learn about Multi

* [Online documentation](https://correaa.github.io/boost-multi/multi/intro.html)

## Use Boost.Multi

Try the library [online](https://godbolt.org/z/occ7Yz78d) with Godbolt's Compiler Explorer.

### Header-only

Add the repository's `include/` directory to your compiler's include path and include the headers you need:

```cpp
#include <multi/array.hpp>
```

Alternatively, download the [single amalgamated header](https://correaa.gitlab.io/boost-multi/boost-multi.hpp) and include it as `#include <boost-multi.hpp>`.

### CMake: in-tree source

Add Multi from an in-tree source directory and link its interface target:

```cmake
add_subdirectory(path/to/boost-multi)
target_link_libraries(my_target PRIVATE multi)
```

### CMake: FetchContent

Fetch Multi on demand:

```cmake
include(FetchContent)
FetchContent_Declare(multi GIT_REPOSITORY https://gitlab.com/correaa/boost-multi.git)
FetchContent_MakeAvailable(multi)
target_link_libraries(my_target PRIVATE multi)
```

### Testing

Testing requires Boost.Test. For example, install `libboost-test-dev` on Debian/Ubuntu or `boost-devel` on Fedora.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Support

* File a Gitlab [issue](https://gitlab.com/correaa/boost-multi/-/issues/new?type=ISSUE) or Github [issue](https://github.com/correaa/boost-multi/issues/new/choose).
* Join the [**#boost-multi**](https://cpplang.slack.com/archives/C071VGKUA5P) discussion group at [cpplang.slack.com](https://cpplang.slack.com/)
