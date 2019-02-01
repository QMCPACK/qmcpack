Just a bogus change
# Getting and building QMCPACK

 Obtain the latest release or development copy from http://www.qmcpack.org

# Prerequisites

 * C++ 14 and C99 capable compilers. e.g. gcc 5.0, Intel 2018, clang 3.5 or newer. 
 * CMake, build utility, http://www.cmake.org
 * BLAS/LAPACK, numerical library, use platform-optimized libraries
 * Libxml2, XML parser, http://xmlsoft.org/
 * HDF5, portable I/O library, http://www.hdfgroup.org/HDF5/
 * BOOST, peer-reviewed portable C++ source libraries, http://www.boost.org , version 1.61.0 or higher.
 * FFTW, FFT library, http://www.fftw.org/

See the manual for more detailed version requirements.

# Building with CMake

 The build system for QMCPACK is based on CMake.  It will autoconfigure
 based on the detected compilers and libraries. Previously QMCPACK made
 extensive use of toolchains, but the system has since been updated to
 eliminate the use of toolchain files for most cases.  The build
 system works with GNU, Intel, and IBM XLC compilers.  Specific compile options
 can be specified either through specific environment or CMake
 variables.  When the libraries are installed in standard locations,
 e.g., /usr, /usr/local, there is no need to set environment or cmake
 variables for the packages.

 See the manuals linked at https://www.qmcpack.org/documentation or buildable via
 manual/build_manual.sh for build examples on Linux, Mac OS X etc.

## Quick build

 If you are feeling lucky and are on a standard UNIX-like system such
 as a Linux workstation:

 * Safest quick build option is to specify the C and C++ compilers
   through their MPI wrappers. Here we use Intel MPI and Intel
   compilers. Move to the build directory, run cmake and make
```
cd build
cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
make -j 8
```

 * Substitute mpicc and mpicxx or other wrapped compiler names to suit
   your system. e.g. With OpenMPI use
```
cd build
cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
make -j 8
```

* If you are feeling particularly lucky, you can skip the compiler
   specification:
```
cd build
cmake ..
make -j 8
```

 The complexities of modern computer hardware and software systems are
 such that you should check that the autoconfiguration system has made
 good choices and picked optimized libraries and compiler settings
 before doing significant production. i.e. Check the details below.

## Set the environment

 A number of environment variables affect the build.  In particular
 they can control the default paths for libraries, the default
 compilers, etc.  The list of environment variables is given below:

| Environment variable | Description |
|----------------------|-------------|
|   CXX          |    C++ compiler |
|   CC           |    C Compiler |
|   MKL_ROOT     |    Path for MKL |
|   LIBXML2_HOME |    Path for libxml2 |
|   HDF5_ROOT    |    Path for HDF5 |
|   BOOST_ROOT   |    Path for Boost |
|   FFTW_HOME    |    Path for FFTW |

## CMake options

 In addition to reading the environment variables, CMake provides a
 number of optional variables that can be set to control the build and
 configure steps.  When passed to CMake, these variables will take
 precident over the environment and default variables.  To set them
 add -D FLAG=VALUE to the configure line between the cmake command and
 the path to the source directory.

 * General build options

```
    CMAKE_C_COMPILER    Set the C compiler
    CMAKE_CXX_COMPILER  Set the C++ compiler
    CMAKE_BUILD_TYPE    A variable which controls the type of build (defaults to Release).
                        Possible values are:
                        None (Do not set debug/optmize flags, use CMAKE_C_FLAGS or CMAKE_CXX_FLAGS)
                        Debug (create a debug build)
                        Release (create a release/optimized build)
                        RelWithDebInfo (create a release/optimized build with debug info)
                        MinSizeRel (create an executable optimized for size)
    CMAKE_C_FLAGS       Set the C flags.  Note: to prevent default debug/release flags
                        from being used, set the CMAKE_BUILD_TYPE=None
                        Also supported: CMAKE_C_FLAGS_DEBUG, CMAKE_C_FLAGS_RELEASE,
                                        CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS     Set the C++ flags.  Note: to prevent default debug/release flags
                        from being used, set the CMAKE_BUILD_TYPE=None
                        Also supported: CMAKE_CXX_FLAGS_DEBUG, CMAKE_CXX_FLAGS_RELEASE,
                                        CMAKE_CXX_FLAGS_RELWITHDEBINFO
```

 * Key QMC build options

```
     QMC_CUDA            Enable CUDA and GPU acceleration (1:yes, 0:no)
     QMC_COMPLEX         Build the complex (general twist/k-point) version (1:yes, 0:no)
     QMC_MIXED_PRECISION Build the mixed precision (mixing double/float) version
                         (1:yes (GPU default), 0:no (CPU default)).
                         The CPU support is experimental.
                         Use float and double for base and full precision.
                         The GPU support is quite mature.
                         Use always double for host side base and full precision
                         and use float and double for CUDA base and full precision.
     ENABLE_TIMERS       Enable fine-grained timers (1:yes, 0:no (default)).
                         Timers are off by default to avoid potential slowdown in small
                         systems. For large systems (100+ electrons) there is no risk.
     ENABLE_SOA          (Experimental) Enable CPU optimization based on Structure-
                         of-Array (SoA) datatypes (1:yes, 0:no (default)). ```
```

 * Additional QMC options

```
     QMC_INCLUDE         Add extra include paths
     QMC_EXTRA_LIBS      Add extra link libraries
     QMC_BUILD_STATIC    Add -static flags to build
     QMC_DATA            Specify data directory for QMCPACK performance and integration tests
     QE_BIN              Location of Quantum Espresso binaries including pw2qmcpack.x
```

  * libxml2 related

```
     Libxml2_INCLUDE_DIRS  Specify include directories for libxml2
     Libxml2_LIBRARY_DIRS  Specify library directories for libxml2
```

* HDF5 related
```
     ENABLE_PHDF5        1(default)/0, enables/disable parallel collective IO.

```

  * FFTW related
```
     FFTW_INCLUDE_DIRS   Specify include directories for FFTW
     FFTW_LIBRARY_DIRS   Specify library directories for FFTW
```

## Example configure and build

In the build directory, run cmake with appropriate options, then
make.

* Using Intel compilers and their MPI wrappers. Assumes HDF5 and
libxml2 will be automatically detected.

```
cd build
cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
make -j 8
```

##  Special notes

It is recommended to create a helper script that contains the
configure line for CMake.  This is particularly useful when using
environment variables, packages are installed in custom locations,
or the configure line may be long or complex.  In this case it is
recommended to add "rm -rf CMake*" before the configure line to remove
existing CMake configure files to ensure a fresh configure each time
that the script is called.  and example script build.sh is given
below:
```
export CXX=mpic++
export CC=mpicc
export ACML_HOME=/opt/acml-5.3.1/gfortran64
export HDF5_ROOT=/opt/hdf5
export BOOST_ROOT=/opt/boost

rm -rf CMake*

cmake                                               \
  -D CMAKE_BUILD_TYPE=Debug                         \
  -D Libxml2_INCLUDE_DIRS=/usr/include/libxml2      \
  -D Libxml2_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu \
  -D FFTW_INCLUDE_DIRS=/usr/include                 \
  -D FFTW_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu    \
  -D QMC_EXTRA_LIBS="-ldl ${ACML_HOME}/lib/libacml.a -lgfortran" \
  -D QMC_DATA=/projects/QMCPACK/qmc-data            \
  ..
```

##  Additional examples:

Set compile flags manually:
```
   cmake                                                \
      -D CMAKE_BUILD_TYPE=None                          \
      -D CMAKE_C_COMPILER=mpicc                         \
      -D CMAKE_CXX_COMPILER=mpic++                      \
      -D CMAKE_C_FLAGS="  -O3 -fopenmp -malign-double -fomit-frame-pointer -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated -march=native -mtune=native" \
      -D CMAKE_CXX_FLAGS="-O3 -fopenmp -malign-double -fomit-frame-pointer -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated -march=native -mtune=native" \
      ..
```

Add extra include directories:
```
   cmake                                                \
      -D CMAKE_BUILD_TYPE=Release                       \
      -D CMAKE_C_COMPILER=mpicc                         \
      -D CMAKE_CXX_COMPILER=mpic++                      \
      -D QMC_INCLUDE="~/path1;~/path2"                  \
      ..
```

# Testing and validation of QMCPACK

Before using QMCPACK we highly encourage tests to be run.
QMCPACK includes extensive validation tests to ensure the correctness of the
code, compilers, tools, and runtime. The tests should ideally be run
each compilation, and certainly before any research use. The tests include
checks of the output against known mean-field, quantum chemistry, and
other QMC results.

While some tests are fully deterministic, due to QMCPACK's stochastic
nature some tests are statistical and can occasionally fail. We employ
a range of test names and labeling to differentiate between these, as
well as developmental tests that are known to fail. In particular,
"deterministic" tests include this in their ctest test name, while
tests known to be unstable (stochastically or otherwise) are labeled
unstable using ctest labels.

For more informaton, consult http://www.qmcpack.org and the manual.
The tests currently use up to 16 cores in various combinations of MPI
tasks and OpenMP threads. Current status for many systems can be
checked at https://cdash.qmcpack.org

Note that due to the small electron and walker counts used in the
tests, they should not be used for any performance measurements. These
should be made on problem sizes that are representative of actual
research calculations. As described in the manual, performance tests
are provided to aid in monitoring performance.

## Run the deterministic tests

From the build directory, invoke ctest specifying only tests
that are deterministic and known to be reliable.
```
ctest -R deterministic -LE unstable
```

These tests currently take a few seconds to run, and include all the
 unit tests. All tests should pass. Failing tests likely indicate a
 significant problem that should be solved before using QMCPACK
 further. This ctest invocation can be used as part of an automated
 installation verification process.
 
## Run the short (quick) tests

 From the build directory, invoke ctest specifying only tests
 including "short" to run that are known to be stable.
```
ctest -R short -LE unstable
```

 These tests currently take up to around one hour. On average, all
 tests should pass at a three sigma level of reliability. Any
 initially failing test should pass when rerun.

## Run individual tests

Individual tests can be run by specifying their name
```
ctest -R name-of-test-to-run
```

# Documentation and support

For more informaton, consult QMCPACK pages at http://www.qmcpack.org,
the manual PDF at https://docs.qmcpack.org/qmcpack_manual.pdf, 
or its sources in the manual directory.

If you have trouble using or building QMCPACK, or have questions about
its use, please post to the [Google QMCPACK group](https://groups.google.com/forum/#!forum/qmcpack) or contact a developer.

# Contributing

Contributions of any size are very welcome. Guidance for contributing
to QMCPACK is included in Chapter 1 of the manual
https://docs.qmcpack.org/qmcpack_manual.pdf . We use a git flow model
including pull request reviews. A continuous integration system runs
on pull requests. See https://github.com/QMCPACK/qmcpack/wiki for
details. For an extensive contribution, it can be helpful to discuss
on the [Google QMCPACK group](https://groups.google.com/forum/#!forum/qmcpack), to create a GitHub issue, or to talk
directly with a developer.

Contributions are made under the same UIUC/NCSA open source license
that covers QMCPACK. Please contact us if this is problematic.
