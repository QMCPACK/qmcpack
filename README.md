# Getting and building QMCPACK

 Obtain the latest release or development copy from http://www.qmcpack.org

# Prerequisites

 * C/C++ compilers
 * CMake, build utility, http://www.cmake.org
 * BLAS/LAPACK, numerical library, use platform-optimized libraries
 * Libxml2, XML parser, http://xmlsoft.org/
 * HDF5, portable I/O library, http://www.hdfgroup.org/HDF5/
 * BOOST, peer-reviewed portable C++ source libraries, http://www.boost.org
 * FFTW, FFT library, http://www.fftw.org/

 Note that the einspline library is no longer required.

# Building with CMake

 The build system for QMCPACK is based on CMake.  It will autoconfigure
 based on the detected compilers and libraries. Previously QMCPACK made
 extensive use of toolchains, but the system has since been updated to
 eliminate the use of toolchain files for most cases.  The build
 system works with GNU, Intel, and IBM XLC compilers.  Specific compile options
 can be specified either through specific environmental or CMake
 variables.  When the libraries are installed in standard locations,
 e.g., /usr, /usr/local, there is no need to set environmental or cmake
 variables for the packages.

 See the manual in manual/qmcpack_manual.pdf for build examples on Linux, Mac OS X etc.

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

 A number of enviornmental variables affect the build.  In particular
 they can control the default paths for libraries, the default
 compilers, etc.  The list of enviornmental variables is given below:

| Environment variable | Description |
|----------------------|-------------|
|   CXX          |    C++ compiler |
|   CC           |    C Compiler |
|   MKL_HOME     |    Path for MKL |
|   LIBXML2_HOME |    Path for libxml2 |
|   HDF5_ROOT    |    Path for HDF5 |
|   BOOST_ROOT   |    Path for Boost |
|   FFTW_HOME    |    Path for FFTW |

## CMake options

 In addition to reading the enviornmental variables, CMake provides a
 number of optional variables that can be set to control the build and
 configure steps.  When passed to CMake, these variables will take
 precident over the enviornmental and default variables.  To set them
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
```

 * Additional QMC options

```
     QMC_DATA            Specify data directory for QMCPACK (currently unused, but
                         likely to be used for performance tests)
     QMC_INCLUDE         Add extra include paths
     QMC_EXTRA_LIBS      Add extra link libraries
     QMC_BUILD_STATIC    Add -static flags to build
```

  * libxml

```
     Libxml2_INCLUDE_DIRS  Specify include directories for libxml2
     Libxml2_LIBRARY_DIRS  Specify library directories for libxml2
```

  * FFTW
```
     FFTW_INCLUDE_DIRS   Specify include directories for FFTW
     FFTW_LIBRARY_DIRS   Specify library directories for FFTW
```

## Configure and build

Move to build directory, run cmake and make
```
cd build
cmake ..
make -j 8
```

## Example configure and build

* Set the environments (the examples below assume bash, Intel compilers and MKL library)
```
export CXX=icpc
export CC=icc
export MKL_HOME=/usr/local/intel/mkl/10.0.3.020
export LIBXML2_HOME=/usr/local
export HDF5_ROOT=/usr/local
export BOOST_ROOT=/usr/local/boost
export FFTW_HOME=/usr/local/fftw
```
* Move to build directory, run cmake and make
```
cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make -j 8
```

##  Special notes

It is recommended to create a helper script that contains the
configure line for CMake.  This is particularly useful when avoiding
enviornmental variables, packages are installed in custom locations,
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

QMCPACK includes validation tests to ensure the correctness of the
code, compilers, tools, and runtime. The tests should ideally be run
each compilation, and certainly before any research use. The tests
check the output against known mean-field, quantum chemistry, and
other QMC results.

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

For more informaton, consult QMCPACK pages at http://www.qmcpack.org and the manual.
The tests currently use up to 16 cores in various combinations of MPI
tasks and OpenMP threads.

Note that due to the small electron and walker counts used in the
tests, they should not be used for any performance measurements. These
should be made on problem sizes that are representative of actual
research calculations.

## Run the short (quick) tests

 From the build directory, invoke ctest specifying only tests
 including "short" should be run
```
ctest -R short
```

 These tests currently take several minutes to run. All tests should pass.

## Run the long verification tests

For greater surety, the long verification tests use a far greater
number of statistical samples than the "short" tests. These take
several hours each to run.

From the build directory, invoke ctest with an increased test timeout
```
ctest --timeout 36000
```

This will run all the defined tests, "short" and "long" as well as the
unit and other tests. If you are running on a system such as a large
shared supercomputer you will likely have to run these tests from
inside a submitted job to avoid run length limits.

## Run individual tests

Individual tests can be run by specifying their name
```
ctest -R name-of-test-to-run
```

# Documentation and support

For more informaton, consult QMCPACK pages at http://www.qmcpack.org,
the linked documentation, and the local copy of the manual in
manual/qmcpack_manual.pdf

If you have trouble using or building QMCPACK, or have questions about
its use, please post to the Google QMCPACK group or contact a developer.
