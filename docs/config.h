////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**\mainpage QMCPACK Documentation
 * \section intro_sec Summary
 *
 *
 * QMCPACK, framework for Quantum Monte Carlo simulations, implements advanced
 * QMC algorithms.  Generic programming enabled by templates in C++ is
 * extensively utilized to achieve high efficiency. It is designed for high
 * performance systems using MPI and OpenMP.
 *
 * The code development is led by J. Kim and the main
 * contributors are the members of the electron structure group of Profs.
 * Martin and Ceperley at University of Illinois at Urbana-Champaign. 
 * - Active developers: J. Kim, K. Esler, J. McMinis, B. Clark, J. Gergely, C. Yang
 * - Past contributors: S. Chiesa, K. Delaney, J. Vincent 
 *
 * The development of QMCPACK is supported by Materials Computation Center and
 * National Center for Supercomputing Applications at University of Illinois at
 * Urbana-Champaign, and funded by the National Science Foundation and
 * Department of Energy.
 *
 * \section change_secs Major changes
 * \htmlonly
 * <h2>2009-06-xx</h2>
 * <h2>2008-07-22</h2>
  <ul> 
  <li>Numerous bug fixes.
  <li>Support TrialWaveFunction cloning for multi-threaded applications
  <li>EinsplineOrbitalSet uses Einspline library 
  <li>BsplinFunctor for One- and Two-Body Jastrow
  </ul>
  <h2> 2005-07-25</h2>
  <ul> 
    <li> updated doxygen documentations
    <li> docs directory is added to cvs repository. 
    <li> To generate doxygen documentation on a local host,
      <ul>
      <li> cd docs; doxygen Doxyfile
      <li> doxygen/html/index.html is the main page
      </ul>
  <li> Introduced groups that follow the directory structure
  <ul>
   <li>Orbital group
   <li>Orbital builder group
   <li>Many-body wave function group
   <li>Hamiltonian group
   <li>QMC Driver group
   <li>QMC Drivers using walker-by-walker update
   <li>QMC Drivers using particle-by-particle update
   <li>QMC Drivers for energy differences
   <li>QMC Application group
  </ul>
  </ul>
 \endhtmlonly

 * \page license University of Illinois/NCSA Open Source License

Copyright (c) 2003, University of Illinois Board of Trustees.
All rights reserved.

Developed by:   
  Jeongnim Kim
  Condensed Matter Physics,
  National Center for Supercomputing Applications, University of Illinois
  Materials computation Center, University of Illinois
  http://www.mcc.uiuc.edu/qmc/

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
``Software''), to deal with the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

        * Redistributions of source code must retain the above copyright 
          notice, this list of conditions and the following disclaimers.
        * Redistributions in binary form must reproduce the above copyright 
          notice, this list of conditions and the following disclaimers in 
          the documentation and/or other materials provided with the 
          distribution.
        * Neither the names of the NCSA, the MCC, the University of Illinois, 
          nor the names of its contributors may be used to endorse or promote 
          products derived from this Software without specific prior written 
          permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS WITH THE SOFTWARE.
 */
#ifndef QMCPLUSPLUS_CONFIGURATION_H
#define QMCPLUSPLUS_CONFIGURATION_H

/* define the major version */
#define QMCPLUSPLUS_VERSION_MAJOR  0

/* define the minor version */
#define QMCPLUSPLUS_VERSION_MINOR  5

/* define the patch version */
#define QMCPLUSPLUS_VERSION_PATCH  0

/* define the release version */
/* #undef QMCPLUSPLUS_RELEASE */

/* define the linearscale version */
/* #undef QMCPLUSPLUS_LINEARSCALE */

/* define the subversion branch */
#define QMCPLUSPLUS_BRANCH  3680

/* define the subversion last changed date */
#define QMCPLUSPLUS_LAST_CHANGED_DATE  "2009-03-13 14:48:19 -0500 (Fri, 13 Mar 2009)"

/* define QMC_BUILD_COMPLETE */
#define QMC_BUILD_COMPLETE 1

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
/* #undef HAVE_EINSPLINE */

/* Define if external einspline is found */
/* #undef HAVE_EINSPLINE_EXT */

#ifndef HAVE_EINSPLINE_EXT

/* Define if posix_memalign function exists */
/* #undef HAVE_POSIX_MEMALIGN */

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
/* #undef HAVE_SSE */

/* Define if sse2 support exists */
/* #undef HAVE_SSE2 */

/* Define if sse3 support exists */
/* #undef HAVE_SSE3 */

/* Define if ssse3 support exists */
/* #undef HAVE_SSSE3 */

/* Define if c variable array support exists */
#define HAVE_C_VARARRAYS 1

/* Prefetch loop lead distance  */
#define PREFETCH_AHEAD 12

/* Use SSE prefetch  */
/* #undef USE_PREFETCH */

#endif /* HAVE_EINSPLINE_EXT */

/* Fund mkl library */
/* #undef HAVE_MKL */

/* Fund mkl/vml library */
/* #undef HAVE_MKL_VML */

/* Fund acml library */
/* #undef HAVE_ACML */

#endif // QMCPLUSPLUS_CONFIGURATION_H

