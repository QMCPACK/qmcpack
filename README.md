![QMCPACK Logo](docs/figs/QMCPACK_logo.png)

[![License](https://img.shields.io/badge/License-UIUC/NCSA-blue.svg)](https://opensource.org/licenses/NCSA)
[![Documentation Status](https://readthedocs.org/projects/qmcpack/badge/?version=develop)](https://qmcpack.readthedocs.io/en/develop/?badge=develop)

[![GitHub release](https://img.shields.io/github/release/QMCPACK/qmcpack/all.svg)](https://github.com/QMCPACK/qmcpack/releases)
[![Spack Version](https://img.shields.io/spack/v/qmcpack)](https://spack.readthedocs.io/en/latest/package_list.html#qmcpack)

[![GitHub Actions CI](https://github.com/QMCPACK/qmcpack/actions/workflows/ci-github-actions.yaml/badge.svg)](https://github.com/QMCPACK/qmcpack/actions/workflows/ci-github-actions.yaml)
[![codecov-deterministic](https://codecov.io/gh/QMCPACK/qmcpack/branch/develop/graph/badge.svg?token=35D0u6GlBm)](https://codecov.io/gh/QMCPACK/qmcpack)

QMCPACK is an open-source production-level many-body ab initio Quantum Monte Carlo code for computing the electronic structure of
atoms, molecules, 2D nanomaterials and solids. The solid-state capabilities include metallic systems as well as insulators.
QMCPACK is expected to run well on workstations through to the latest generation supercomputers. Besides high performance,
particular emphasis is placed on code quality and reproducibility.

# Obtaining and installing QMCPACK

Obtain the latest release from https://github.com/QMCPACK/qmcpack/releases or clone the development source from
https://github.com/QMCPACK/qmcpack. A full installation guide and steps to perform an initial QMC calculation are given in the
[extensive online documentation for QMCPACK](https://qmcpack.readthedocs.io/en/develop/index.html).

The [CHANGELOG.md](CHANGELOG.md) describes key changes made in each release as well as any major changes to the development version.

# Documentation and support

For more information, consult QMCPACK pages at http://www.qmcpack.org, the manual at
https://qmcpack.readthedocs.io/en/develop/index.html, or its sources in the docs directory.

If you have trouble using or building QMCPACK, or have questions about its use, please post to the [Google QMCPACK
group](https://groups.google.com/forum/#!forum/qmcpack), create a GitHub issue at https://github.com/QMCPACK/qmcpack/issues or
contact a developer.

# Learning about Quantum Monte Carlo

To learn about the fundamentals of Quantum Monte Carlo through to their practical application to molecular and solid-state systems with
QMCPACK, see the [materials and tutorials from our most recent QMC workshop](https://github.com/QMCPACK/qmc_workshop_2021). These include a virtual machine
to run examples without having to install QMCPACK yourself, and slides and recorded videos of introductory talks through to spin-orbit QMC.

# Citing QMCPACK

Please cite J. Kim _et al._ J. Phys. Cond. Mat. **30** 195901 (2018), https://doi.org/10.1088/1361-648X/aab9c3, and if space allows,
P. Kent _et al._ J. Chem. Phys. **152** 174105 (2020), https://doi.org/10.1063/5.0004860 . These papers are both open access.

# Installing QMCPACK

[Complete instructions for downloading, compiling, and installing QMCPACK are given in the manual](https://qmcpack.readthedocs.io/en/develop/installation.html).
Guides are provided for [installing QMCPACK on common workstations and supercomputers](https://qmcpack.readthedocs.io/en/develop/installation.html#installation-instructions-for-common-workstations-and-supercomputers). We also provide build
scripts are for many systems in the config directory. The following sections of this README give a brief overview. 

# Installation Prerequisites

 * C++ 17 and C99 capable compilers. 
 * CMake v3.21.0 or later, build utility, http://www.cmake.org
 * BLAS/LAPACK, numerical library. Use vendor and platform-optimized libraries.
 * LibXml2, XML parser, http://xmlsoft.org/
 * HDF5 v1.10.0 or later, portable I/O library, http://www.hdfgroup.org/HDF5/
 * BOOST v1.70.0 or newer, peer-reviewed portable C++ source libraries, http://www.boost.org
 * FFTW, FFT library, http://www.fftw.org/
 * MPI, parallel library. Optional, but a near requirement for production calculations.
 * Python3. Older versions are not supported as of January 2020.
 * CUDA v11.0 or later. Optional, but required for builds with NVIDIA GPU support. Use 12.3 or newer if possible. 11.3-12.2 have
   a bug affecting multideterminant calculations. Single determinant calculations are OK.

We aim to support open source compilers and libraries released within two years of each QMCPACK release. Use of software versions
over two years old may work but is discouraged and untested. Proprietary compilers (Intel, NVHPC) are generally supported over the
same period but may require use of an exact version. We also aim to support the standard software environments on machines such as
Frontier and Summit at OLCF, Aurora and Polaris at ALCF, and Perlmutter at NERSC. Use of the most recently released compilers and
library versions is particularly encouraged for highest performance and easiest configuration.

Details of regularly tested versions are given in the [manual section on installation
prerequisites](https://qmcpack.readthedocs.io/en/develop/installation.html#prerequisites).

# Building with CMake

The build system for QMCPACK is based on CMake.  It will auto-configure based on the detected compilers and libraries. When these 
are installed in standard locations, e.g., /usr, /usr/local, there is no need to set either environment or CMake variables.

See the manual linked at https://qmcpack.readthedocs.io/en/develop/ and https://www.qmcpack.org/documentation or buildable using
sphinx from the sources in docs/. A PDF version is still available at https://qmcpack.readthedocs.io/_/downloads/en/develop/pdf/

## Quick build

On a standard UNIX-like system such as a Linux workstation:

* Safest quick build option is to specify the C and C++ compilers
  through their MPI wrappers. Here we use Intel MPI and Intel
  compilers. Move to the build directory, run CMake and make
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

* Non-MPI build:
```
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DQMC_MPI=0 ..
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
such that you should check that the auto-configuration system has made
good choices and picked optimized libraries and compiler settings
before doing significant production. i.e. Check the details below.

## Set the environment

 A number of environment variables affect the build.  In particular,
 they can control the default paths for libraries, the default
 compilers, etc. The list of environment variables is given below:

| Environment variable | Description |
|----------------------|-------------|
|   CXX          |    C++ compiler |
|   CC           |    C Compiler |
|   MKL_ROOT     |    Path for MKL |
|   HDF5_ROOT    |    Path for HDF5 |
|   BOOST_ROOT   |    Path for Boost |
|   FFTW_HOME    |    Path for FFTW |

## CMake options

 In addition to reading the environment variables, CMake provides a
 number of optional variables that can be set to control the build and
 configure steps.  When passed to CMake, these variables will take
 precedent over the environment and default variables.  To set them
 add -D FLAG=VALUE to the configure line between the CMake command and
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
    CMAKE_SYSTEM_NAME   Set value to CrayLinuxEnvironment when cross-compiling
                        in Cray Programming Environment.
    CMAKE_C_FLAGS       Set the C flags.  Note: to prevent default debug/release flags
                        from being used, set the CMAKE_BUILD_TYPE=None
                        Also supported: CMAKE_C_FLAGS_DEBUG, CMAKE_C_FLAGS_RELEASE,
                                        CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS     Set the C++ flags.  Note: to prevent default debug/release flags
                        from being used, set the CMAKE_BUILD_TYPE=None
                        Also supported: CMAKE_CXX_FLAGS_DEBUG, CMAKE_CXX_FLAGS_RELEASE,
                                        CMAKE_CXX_FLAGS_RELWITHDEBINFO
```

 * Key QMCPACK build options

```
    QMC_COMPLEX           ON/OFF(default). Build the complex (general twist/k-point) version.
    QMC_MIXED_PRECISION   ON/OFF(default). Build the mixed precision (mixing double/float) version
                          Mixed precision calculations can be signifiantly faster but should be
                          carefully checked validated against full double precision runs,
                          particularly for large electron counts.
    QMC_GPU               Semicolon-separated list of GPU features to build (openmp,cuda,hip,sycl).
                          "openmp", "cuda", "hip" and "sycl" for GPU acceleration via OpenMP offload, CUDA, HIP and SYCL.
                          Recommended values: "openmp;cuda" for NVIDIA, "openmp;hip" for AMD, "openmp;sycl" for Intel.
                          Its default value is set to the recommended value if QMC_GPU_ARCHS indicates a specific vendor
                          or left empty otherwise.
    QMC_GPU_ARCHS         Specify GPU architectures. For example, "gfx90a" targets AMD MI200 series GPUs.
                          "intel_gpu_pvc" targets Intel Data Center GPU Max 1xxx.
                          "sm_80;sm_70" creates a single executable running on both NVIDIA A100 and V100 GPUs.
                          Mixing vendor "gfx90a;sm_70" is not supported. If not set, atempt to derive it
                          from CMAKE_CUDA_ARCHITECTURES or CMAKE_HIP_ARCHITECTURES if available and then
                          atempt to auto-detect existing GPUs.
```

 * Additional QMCPACK options

```
     BUILD_AFQMC            ON/OFF(default). Build the Auxiliary-Field Quantum Monte Carlo (AFQMC) feature
     BUILD_AFQMC_WITH_NCCL  ON/OFF(default). Enable the optimized code path using NVIDIA Collective Communications Library (NCCL) in AFQMC.
                            AFQMC and CUDA features required to enable this feature.
     BUILD_AFQMC_HIP        ON/OFF(default). Enable HIP accelerated code paths in AFQMC. AFQMC feature required to enable this feature.
     QE_BIN                 Location of Quantum Espresso binaries including pw2qmcpack.x
     RMG_BIN                Location of RMG binary
     QMC_DATA               Specify data directory for QMCPACK performance and integration tests
     QMC_INCLUDE            Add extra include paths
     QMC_EXTRA_LIBS         Add extra link libraries
     QMC_BUILD_STATIC       ON/OFF(default). Add -static flags to build
     QMC_SYMLINK_TEST_FILES Set to zero to require test files to be copied. Avoids space
                            saving default use of symbolic links for test files. Useful
                            if the build is on a separate filesystem from the source, as
                            required on some HPC systems.
     ENABLE_TIMERS          ON(default)/OFF. Enable fine-grained timers. Timers are on by default but at level coarse
                            to avoid potential slowdown in tiny systems.
                            For systems beyond tiny sizes (100+ electrons) there is no risk.
     USE_OBJECT_TARGET      ON/OFF(default). Use CMake object library targets to workaround linker not being able to handle hybrid
                            binary archives which contain both host and device codes.
```

  * libxml2 related

```
     LIBXML2_INCLUDE_DIR Include directory for libxml2
     LIBXML2_LIBRARY     Libxml2 library
```

* HDF5 related
```
     HDF5_PREFER_PARALLEL 1(default for MPI build)/0, enables/disable parallel HDF5 library searching.
     ENABLE_PHDF5         1(default for parallel HDF5 library)/0, enables/disable parallel collective I/O.

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
export HDF5_ROOT=/opt/hdf5
export BOOST_ROOT=/opt/boost

rm -rf CMake*

cmake                                               \
  -D CMAKE_BUILD_TYPE=Debug                         \
  -D LIBXML2_INCLUDE_DIR=/usr/include/libxml2      \
  -D LIBXML2_LIBRARY=/usr/lib/x86_64-linux-gnu/libxml2.so \
  -D FFTW_INCLUDE_DIRS=/usr/include                 \
  -D FFTW_LIBRARY_DIRS=/usr/lib/x86_64-linux-gnu    \
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

We highly encourage tests to be run before using QMCPACK. Details are given in the [QMCPACK
manual](https://qmcpack.readthedocs.io/en/develop/index.html). QMCPACK includes extensive validation tests to ensure the
correctness of the code, compilers, tools, and runtime. The tests should ideally be run each compilation, and certainly before any
research use. The tests include checks of the output against known mean-field, quantum chemistry, and other QMC results.

While some tests are fully deterministic, due to QMCPACK's stochastic nature some tests are statistical and can occasionally fail.
We employ a range of test names and labeling to differentiate between these, as well as developmental tests that are known to
fail. In particular, "deterministic" tests include this in their ctest test name, while tests known to be unstable (stochastically
or otherwise) are labeled unstable using ctest labels.

The tests currently use up to 16 cores in various combinations of MPI tasks and OpenMP threads. Current status for many
combinations of systems, compilers, and libraries can be checked at https://cdash.qmcpack.org

Note that due to the small electron and walker counts used in the tests, they should not be used for any performance measurements.
These should be made on problem sizes that are representative of actual research calculations. As described in the manual,
performance tests are provided to aid in monitoring performance.

## Run the unit tests

From the build directory, invoke ctest specifying only the unit tests
```
ctest -j 16 -R unit --output-on-failure
```
All of these tests should pass within a few minutes. Modify the parallization setting (-j 16) to suit the core count of your system.

## Run the deterministic tests

From the build directory, invoke ctest specifying only tests
that are deterministic and known to be reliable.
```
ctest -j 16 -R deterministic -LE unstable --output-on-failure
```

These tests currently take a few minutes to run, and include all the unit tests. All tests should pass. Failing tests likely
indicate a significant problem that should be solved before using QMCPACK further. This ctest invocation can be used as part of an
automated installation verification process. Many of the tests use a multiple of 16 processes, so on large core count machines
a significant speedup can be obtained with -j 64 etc.
 
## Run the short (quick) tests

From the build directory, invoke ctest specifying only tests
including "short" to run that are known to be stable.
```
ctest -j 16 -R short -LE unstable --output-on-failure
```

These tests currently take up to around one hour. On average, all
tests should pass at a three sigma level of reliability. Any
initially failing test should pass when rerun.

## Run individual tests

Individual tests can be run by specifying their name
```
ctest -R name-of-test-to-run
```

# Contributing

Contributions of any size are very welcome. Guidance for contributing to QMCPACK is included in the manual
https://qmcpack.readthedocs.io/en/develop/introduction.html#contributing-to-qmcpack. We use a git flow model including pull
request reviews. A continuous integration system runs on pull requests. See https://github.com/QMCPACK/qmcpack/wiki for details.
For an extensive contribution, it can be helpful to discuss on the [Google QMCPACK
group](https://groups.google.com/forum/#!forum/qmcpack), to create a GitHub issue, or to talk directly with a developer in
advance.

Contributions are made under the same UIUC/NCSA open source license that covers QMCPACK. Please contact us if this is problematic.


