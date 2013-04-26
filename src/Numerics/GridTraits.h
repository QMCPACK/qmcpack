//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
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
/** @file GridTraits.h
 *
 * Define data types for any GridType
 */
#ifndef QMCPLUSPLUS_ONEDIMENSIONALGRID_TRAITS_H
#define QMCPLUSPLUS_ONEDIMENSIONALGRID_TRAITS_H
#include <vector>
#include <complex>
#include <limits>

/** enumeration of one-dimensional grid type
 */
enum {LINEAR_1DGRID, LOG_1DGRID, LOGZERO_1DGRID, CUSTOM_1DGRID};

/** enumeration of boundary conditions
 */
enum {PBC_CONSTRAINTS, FIRSTDERIV_CONSTRAINTS, NATURAL_CONSTRAINTS};

template <class T> struct GridTraits {};

template<>
struct GridTraits<double>
{
  typedef double          point_type;
  typedef double          value_type;
};

template<>
struct GridTraits<std::complex<double> >
{
  typedef double               point_type;
  typedef std::complex<double> value_type;
};
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1580 $   $Date: 2007-01-04 10:00:43 -0600 (Thu, 04 Jan 2007) $
 * $Id: TricubicBsplineSet.h 1580 2007-01-04 16:00:43Z jnkim $
 ***************************************************************************/
