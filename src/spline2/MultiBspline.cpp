//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file MultiBspline.cpp
 *
 * Initialize the static data for float and double
 */
#include <spline2/MultiBspline.hpp>

namespace qmcplusplus
{
  /** initialization of static data for MultiBspline<float> */
  template<>
  QMC_ALIGNAS const float MultiBspline<float>::A44[16] = 
  { -1.0f/6.0f,  3.0f/6.0f, -3.0f/6.0f, 1.0f/6.0f,
     3.0f/6.0f, -6.0f/6.0f,  0.0f/6.0f, 4.0f/6.0f,
    -3.0f/6.0f,  3.0f/6.0f,  3.0f/6.0f, 1.0f/6.0f,
     1.0f/6.0f,  0.0f/6.0f,  0.0f/6.0f, 0.0f/6.0f };

  template<>
  QMC_ALIGNAS const float MultiBspline<float>::dA44[16] = 
  {  0.0f, -0.5f,  1.0f, -0.5f,
     0.0f,  1.5f, -2.0f,  0.0f,
     0.0f, -1.5f,  1.0f,  0.5f,
     0.0f,  0.5f,  0.0f,  0.0f };

  template<>
  QMC_ALIGNAS const float MultiBspline<float>::d2A44[16] =
  {  0.0f, 0.0f, -1.0f,  1.0f,
     0.0f, 0.0f,  3.0f, -2.0f,
     0.0f, 0.0f, -3.0f,  1.0f,
     0.0f, 0.0f,  1.0f,  0.0f };

  /** initialization of static data for MultiBspline<double> */
  template<>
  QMC_ALIGNAS const double MultiBspline<double>::A44[16] = 
  { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 };

  template<>
  QMC_ALIGNAS const double MultiBspline<double>::dA44[16] = 
  {  0.0, -0.5,  1.0, -0.5,
     0.0,  1.5, -2.0,  0.0,
     0.0, -1.5,  1.0,  0.5,
     0.0,  0.5,  0.0,  0.0 };

  template<>
  QMC_ALIGNAS const double MultiBspline<double>::d2A44[16] =
  {  0.0, 0.0, -1.0,  1.0,
     0.0, 0.0,  3.0, -2.0,
     0.0, 0.0, -3.0,  1.0,
     0.0, 0.0,  1.0,  0.0 };
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/
