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
