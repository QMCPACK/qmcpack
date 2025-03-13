//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TMOVEKIND_H
#define QMCPLUSPLUS_TMOVEKIND_H

namespace qmcplusplus
{
/// Tmove options
enum class TmoveKind
{
  OFF = 0, // no Tmove
  V0,      // M. Casula, PRB 74, 161102(R) (2006)
  V1,      // version 1, M. Casula et al., JCP 132, 154113 (2010)
  V3,      // an approximation to version 1 but much faster.
};
} // namespace qmcplusplus
#endif
