//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DATALOCALITY_H
#define QMCPLUSPLUS_DATALOCALITY_H

namespace qmcplusplus
{

///data locality with respect to walker buffer
enum class DataLocality
{
  process = 0,
  rank,
  crowd,
  queue,
  walker
};

}

#endif
