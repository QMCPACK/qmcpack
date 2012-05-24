/* src/ohmms-config.h.in.  Generated from configure.in by autoheader.  */
// -*- c++  -*-
//
//Ohmms Configuration Header. Automatically Generated
//
//See the LICENSE file in the top-level directory for copyright notices
//
#ifndef QMCPLUSPLUS_CONFIGURATION_H
#define QMCPLUSPLUS_CONFIGURATION_H

/* define the major version */
#define QMCPLUSPLUS_VERSION_MAJOR  0

/* define the minor version */
#define QMCPLUSPLUS_VERSION_MINOR  6

/* define the patch version */
#define QMCPLUSPLUS_VERSION_PATCH  1

/* define the release version */
/* #undef QMCPLUSPLUS_RELEASE */

/* define the linearscale version */
/* #undef QMCPLUSPLUS_LINEARSCALE */

/* define the subversion branch */
#define QMCPLUSPLUS_BRANCH  5487

/* define the subversion last changed date */
#define QMCPLUSPLUS_LAST_CHANGED_DATE  "2012-04-22 10:24:14 -0400 (Sun, 22 Apr 2012)"

/* define QMC_BUILD_LEVEL */
#define QMC_BUILD_LEVEL 3

/* define PRINT_DEBUG */
/* #undef PRINT_DEBUG */

/* Enable OpenMP parallelization. */
#define ENABLE_OPENMP 1

/* Define to 1 if you have the `hdf5' library (-lhdf5). */
#define HAVE_LIBHDF5 1

/* Define to 1 if you want to use parallel hdf5 for frequent output */
/* #undef ENABLE_PHDF5 */

/* Define to 1 if you have the `boost' library */
#define HAVE_LIBBOOST 1

/* Define to 1 if you have the `sprng' library (-lsprng). */
/* #undef HAVE_LIBSPRNG */

/* Define to 1 if you have the `blitz' library */
/* #undef HAVE_LIBBLITZ */

/* Define to 1 if you have libxml2 */
#define HAVE_LIBXML2 1

/* Define to 1 if you have fftw */
#define HAVE_LIBFFTW 1

/* Define to 1 if you have libxml++ */
/* #undef HAVE_LIBXMLPP */

/* Define to 1 if you have gsl */
/* #undef HAVE_LIBGSL */

/* Define to 1 if you have MPI library */
/* #undef HAVE_MPI */

/* Define to 1 if you have OOMPI library */
/* #undef HAVE_OOMPI */

/* Define the base precision: float, double */
#define APP_PRECISION double

/* Define the physical dimension of appliation. */
#define OHMMS_DIM 3

/* Define the index type: int, long */
#define OHMMS_INDEXTYPE int

/* Define the base precision: float, double */
#define OHMMS_PRECISION double

/* Define to 1 if complex wavefunctions are used */
/* #undef QMC_COMPLEX */

/* Define to 1 if using AYSNC comm for estimator */
/* #undef QMC_ASYNC_COLLECT */

/* Define to 1 if using recursive SK evaluation */
/* #undef QMC_SK_USE_RECURSIVE */

/* Define if the code is specialized for orthorhombic supercell */
#define OHMMS_ORTHO 0

/* Define the index of the walker iterator. NOT USED */
#define QMC_FASTWALKER 1

/* Define if sincos function exists */
/* #undef HAVE_SINCOS */

/* Define if std::round function exists */
#define HAVE_STD_ROUND 1

/* Define if floor function exists */
#define HAVE_FLOOR 1

/* Define if einspline lib exists */
#define HAVE_EINSPLINE 1

/* Define if external einspline is found */
/* #undef HAVE_EINSPLINE_EXT */

#ifndef HAVE_EINSPLINE_EXT

/* Define if posix_memalign function exists */
#define HAVE_POSIX_MEMALIGN 1

/* Define if pow function exists */
#define HAVE_POW 1

/* Define if sqrt function exists */
#define HAVE_SQRT 1

/* Define if dlfcn.h exists */
#define HAVE_DLFCN_H 1

/* Define if inttypes.h exists */
#define HAVE_INTTYPES_H 1

/* Define if memory.h exists */
#define HAVE_MEMORY_H 1

/* Define if pmmintrin.h exists */
#define HAVE_PMMINTRIN_H 1

/* Define if emmintrin.h exists */
#define HAVE_EMMINTRIN_H 1

/* Define if sys/stat.h exists */
#define HAVE_SYS_STAT_H 1

/* Define if sys/time.h exists */
#define HAVE_SYS_TIME_H 1

/* Define if sys/types.h exists */
#define HAVE_SYS_TYPES_H 1

/* Define if unistd.h exists */
#define HAVE_UNISTD_H 1

/* Define if mmx support exists */
/* #undef HAVE_MMX */

/* Define if sse support exists */
#define HAVE_SSE 1

/* Define if sse2 support exists */
#define HAVE_SSE2 1

/* Define if sse3 support exists */
#define HAVE_SSE3 1

/* Define if ssse3 support exists */
#define HAVE_SSSE3 1

/* Define if c variable array support exists */
/* #undef HAVE_C_VARARRAYS */

/* Prefetch loop lead distance  */
#define PREFETCH_AHEAD 10

/* Use SSE prefetch  */
#define USE_PREFETCH 1

#endif /* HAVE_EINSPLINE_EXT */

/* Find mkl library */
/* #undef HAVE_MKL */

/* Find mkl/vml library */
/* #undef HAVE_MKL_VML */

/* Find essl library */
/* #undef HAVE_ESSL */

/* Fund acml library */
/* #undef HAVE_ACML */

/* Using CUDA for GPU execution */
/* #undef QMC_CUDA */

/* Setting precision for CUDA core kernels */
#define CUDA_PRECISION float

/* #undef DEBUG_PSIBUFFER_ON */

#endif // QMCPLUSPLUS_CONFIGURATION_H

