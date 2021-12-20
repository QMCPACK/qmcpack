//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_INPUTNODE_HPP
#define QMCPLUSPLUS_INPUTNODE_HPP

#include <memory>

namespace qmcplusplus
{

/** Interface to allow Input nodes to be type erased from the standship of memory ownership.
 *  Inputs can greatly vary in size so we don't want to try to do this with variant.
 */
struct InputNode {
  virtual ~InputNode() {};
};

}
#endif
