//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Lab
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMC_MACROS_H
#define QMCPLUSPLUS_QMC_MACROS_H


#define QMC_FOREACH_REAL_TYPE(MACRO) \
  MACRO(float)                       \
  MACRO(double)

#define QMC_FOREACH_COMPLEX_TYPE(MACRO) \
  MACRO(std::complex<float>)            \
  MACRO(std::complex<double>)

#define QMC_FOREACH_TYPE(MACRO) \
  QMC_FOREACH_REAL_TYPE(MACRO)  \
  QMC_FOREACH_COMPLEX_TYPE(MACRO)

#endif
