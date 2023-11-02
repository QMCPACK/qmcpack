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
#define QMCPACK_VERSION_MAJOR  3

/* define the minor version */
#define QMCPACK_VERSION_MINOR  17

/* define the patch version */
#define QMCPACK_VERSION_PATCH  9

/* define the release version */
/* #undef QMCPACK_RELEASE */

/* define the git last commit date */
/* #undef QMCPLUSPLUS_LAST_CHANGED_DATE */

/* building from Git repository or not */
#define IS_GIT_PROJECT  1

/* define QMC_BUILD_SANDBOX_ONLY */
/* #undef QMC_BUILD_SANDBOX_ONLY */

/* define PRINT_DEBUG */
/* #undef PRINT_DEBUG */

/* Enable OpenMP offload. */
/* #undef ENABLE_OFFLOAD */

#ifdef ENABLE_OFFLOAD
  #define PRAGMA_OFFLOAD(x) _Pragma(x)
#else
  #define PRAGMA_OFFLOAD(x)
#endif

/* Enable OpenMP taskloop. */
#define ENABLE_OMP_TASKLOOP ON

#ifdef ENABLE_OMP_TASKLOOP
  #define PRAGMA_OMP_TASKLOOP(x) _Pragma(x)
#else
  #define PRAGMA_OMP_TASKLOOP(x) _Pragma("omp taskgroup")
#endif

/* Manage memory allocations via vendor APIs and associate them with the OpenMP runtime */
/* #undef QMC_OFFLOAD_MEM_ASSOCIATED */

/* Define to 1 if you have MPI library */
#define HAVE_MPI 1

/* Define the physical dimension of appliation. */
#define OHMMS_DIM 3

/* Define the index type: int, long */
#define OHMMS_INDEXTYPE int

/* Define the base precision: float, double */
#define OHMMS_PRECISION double

/* Define the full precision: double, long double */
#define OHMMS_PRECISION_FULL double

/* Define Cache/SIMD alignment in bytes */
#define QMC_SIMD_ALIGNMENT 64

/* Define to 1 if precision is mixed, only for the CPU code */
/* #undef MIXED_PRECISION */

/* Define to 1 if complex wavefunctions are used */
#define QMC_COMPLEX 1

/* Define if sincos function exists */
#define HAVE_SINCOS 1

/* Define if external einspline is found */
/* #undef HAVE_EINSPLINE_EXT */

#ifndef HAVE_EINSPLINE_EXT

/* Define if posix_memalign function exists */
#define HAVE_POSIX_MEMALIGN 1

/* Define if pmmintrin.h exists */
/* #undef HAVE_PMMINTRIN_H */

/* Define if emmintrin.h exists */
/* #undef HAVE_EMMINTRIN_H */

#endif /* HAVE_EINSPLINE_EXT */

/* Find ESSL library */
/* #undef HAVE_ESSL */

/* Using translation of CUDA to HIP for GPU execution */
/* #undef QMC_CUDA2HIP */

/* Disable hipHostRegister/hipHostUnregister */
/* #undef QMC_DISABLE_HIP_HOST_REGISTER */

/* Using CUDA for GPU execution, next generation */
/* #undef ENABLE_CUDA */

/* Using SYCL for GPU execution */
/* #undef ENABLE_SYCL */

/* Using boost::stacktrace */
/* #undef ENABLE_STACKTRACE */

/* For AFQMC compilation  */
/* #undef BUILD_AFQMC */

/* For AFQMC compilation  */
/* #undef BUILD_AFQMC_WITH_NCCL */

/* For FCIQMC compilation  */
/* #undef BUILD_FCIQMC */

/* #undef DEBUG_PSIBUFFER_ON */

/* Disable trace manager and associated features */
/* #undef DISABLE_TRACEMANAGER */

/* Fully remove trace manager and associated features */
/* #undef REMOVE_TRACEMANAGER */

/* Fixed Size Walker Properties */
#define WALKER_MAX_PROPERTIES 2048

/* Internal timers */
#define ENABLE_TIMERS 1

/* Use VTune API */
/* #undef USE_VTUNE_API */

/* Use VTune Task API with timers */
/* #undef USE_VTUNE_TASKS */

/* Enable NVTX regions in CUDA code. */
/* #undef USE_NVTX_API */

#endif // QMCPLUSPLUS_CONFIGURATION_H

