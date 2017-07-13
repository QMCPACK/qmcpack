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
  /** initialization of static data for MultiBsplineData<float> */
  template<>
  QMC_ALIGNAS const float MultiBsplineData<float>::A44[16] =
  { -1.0f/6.0f,  3.0f/6.0f, -3.0f/6.0f, 1.0f/6.0f,
     3.0f/6.0f, -6.0f/6.0f,  0.0f/6.0f, 4.0f/6.0f,
    -3.0f/6.0f,  3.0f/6.0f,  3.0f/6.0f, 1.0f/6.0f,
     1.0f/6.0f,  0.0f/6.0f,  0.0f/6.0f, 0.0f/6.0f };

  template<>
  QMC_ALIGNAS const float MultiBsplineData<float>::dA44[16] =
  {  0.0f, -0.5f,  1.0f, -0.5f,
     0.0f,  1.5f, -2.0f,  0.0f,
     0.0f, -1.5f,  1.0f,  0.5f,
     0.0f,  0.5f,  0.0f,  0.0f };

  template<>
  QMC_ALIGNAS const float MultiBsplineData<float>::d2A44[16] =
  {  0.0f, 0.0f, -1.0f,  1.0f,
     0.0f, 0.0f,  3.0f, -2.0f,
     0.0f, 0.0f, -3.0f,  1.0f,
     0.0f, 0.0f,  1.0f,  0.0f };

  /** initialization of static data for MultiBsplineData<double> */
  template<>
  QMC_ALIGNAS const double MultiBsplineData<double>::A44[16] =
  { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 };

  template<>
  QMC_ALIGNAS const double MultiBsplineData<double>::dA44[16] =
  {  0.0, -0.5,  1.0, -0.5,
     0.0,  1.5, -2.0,  0.0,
     0.0, -1.5,  1.0,  0.5,
     0.0,  0.5,  0.0,  0.0 };

  template<>
  QMC_ALIGNAS const double MultiBsplineData<double>::d2A44[16] =
  {  0.0, 0.0, -1.0,  1.0,
     0.0, 0.0,  3.0, -2.0,
     0.0, 0.0, -3.0,  1.0,
     0.0, 0.0,  1.0,  0.0 };
}

