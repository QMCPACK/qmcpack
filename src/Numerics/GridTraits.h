//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
