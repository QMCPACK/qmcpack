//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COULOMB_TYPE_HPP
#define QMCPLUSPLUS_COULOMB_TYPE_HPP
#include <config.h>
//model to implement the mixed precision

#define DECLARE_COULOMB_TYPES                              \
  using pRealType    = OHMMS_PRECISION;                    \
  using mRealType    = OHMMS_PRECISION_FULL;               \
  using PosType      = TinyVector<pRealType, OHMMS_DIM>;
#endif
