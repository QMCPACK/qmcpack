//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ONEDIMGRIDFACTORY_H
#define QMCPLUSPLUS_ONEDIMGRIDFACTORY_H
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/LibxmlNumericIO.h"

namespace qmcplusplus
{
/** Factory class using Singleton pattern
 */
template <typename T>
struct OneDimGridFactory
{
  using RealType = T;
  ///typedef of the one-dimensional grid
  using GridType = OneDimGridBase<RealType>;

  /** return a GridType*
   * @param cur xmlnode for the grid definition
   */
  static std::unique_ptr<GridType> createGrid(xmlNodePtr cur);
};
} // namespace qmcplusplus
#endif
