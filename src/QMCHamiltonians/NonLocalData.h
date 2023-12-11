//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NONLOCALDATA_H
#define QMCPLUSPLUS_NONLOCALDATA_H

#include "Configuration.h"

namespace qmcplusplus
{
struct NonLocalData : public QMCTraits
{
  IndexType PID;
  RealType Weight;
  PosType Delta;
  inline NonLocalData() : PID(-1), Weight(1.0) {}
  inline NonLocalData(IndexType id, RealType w, const PosType& d) : PID(id), Weight(w), Delta(d) {}
};
} // namespace qmcplusplus
#endif
