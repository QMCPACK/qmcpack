//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from VMC.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CROWD_H
#define QMCPLUSPLUS_CROWD_H

namespace qmcplusplus
{
class Crowd
{
public:
  /** This is the data structure for walkers within a crowd
   */
  struct Walkers
  {
  public:
    int live;
  };

private:
  Walkers walkers_;

public:
};
} // namespace qmcplusplus
#endif
