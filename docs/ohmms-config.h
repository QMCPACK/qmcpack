////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**\mainpage qmcPlusPlus main page
 * \section intro_sec Introduction
 *
 * \section license
 *
 */
#ifndef OHMMS_QMC_CONFIGURATION_H
#define OHMMS_QMC_CONFIGURATION_H

/* Enable OpenMP parallelization. */
/* #undef ENABLE_OPENMP 0 */

/* Define to 1 if you have the `hdf5' library (-lhdf5). */
#define HAVE_LIBHDF5 1

/* Define to 1 if you have the `boost' library */
#define HAVE_LIBBOOST 1

/* Define to 1 if you have the `sprng' library (-lsprng). */
/* #undef HAVE_LIBSPRNG 0 */

/* Define to 1 if you have the `blitz' library */
#define HAVE_LIBBLITZ 1

/* Define to 1 if you have libxml2 */
#define HAVE_LIBXML2 1

/* Define to 1 if you have libxml++ */
/* #undef HAVE_LIBXMLPP 0 */

/* Define to 1 if you have gsl */
#define HAVE_LIBGSL 1

/* Define to 1 if you have MPI library */
/* #undef HAVE_MPI 0 */

/* Define to 1 if you have OOMPI library */
/* #undef HAVE_OOMPI 0 */

/* Define the physical dimension of appliation. */
#define OHMMS_DIM 3

/* Define the index type: int, long */
#define OHMMS_INDEXTYPE int

/* Define the base precision: float, double */
#define OHMMS_PRECISION double

/* Define if the code is specialized for orthorhombic supercell */
#define OHMMS_ORTHO 0

/* Define the index of the walker iterator. NOT USED */
#define QMC_FASTWALKER 1

#endif // OHMMS_QMC_CONFIGURATION_H
