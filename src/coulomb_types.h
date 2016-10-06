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
#if (__cplusplus >= 201103L)

#define DECLARE_COULOMB_TYPES \
  using pRealType=OHMMS_PRECISION;              \
  using mRealType=OHMMS_PRECISION_FULL;                       \
  using pComplexType=std::complex<pRealType>;   \
  using mComplexType=std::complex<OHMMS_PRECISION_FULL>;      \
  using PosType=TinyVector<pRealType,OHMMS_DIM>; 

#else

#define DECLARE_COULOMB_TYPES \
  typedef OHMMS_PRECISION pRealType;               \
  typedef OHMMS_PRECISION_FULL mRealType;                        \
  typedef std::complex<pRealType> pComplexType;    \
  typedef std::complex<OHMMS_PRECISION_FULL> mComplexType;       \
  typedef TinyVector<pRealType,OHMMS_DIM> PosType;

#endif

#endif
