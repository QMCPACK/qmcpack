//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
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

/** Interface to allow input nodes to be type erased for the purposes of ownership
 */
struct InputNode
{
  virtual ~InputNode(){};
};

} // namespace qmcplusplus
#endif
