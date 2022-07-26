//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HEGGRID_H
#define QMCPLUSPLUS_HEGGRID_H

#include "Lattice/CrystalLattice.h"
#include <array>

namespace qmcplusplus
{
//three-d specialization
template<class T>
struct HEGGrid
{
  using PL_t = CrystalLattice<T, OHMMS_DIM>;

  const PL_t& Lattice;
  static constexpr std::array<int, 31> n_within_shell{1,   7,   19,  27,  33,  57,  81,  93,  123, 147, 171,
                                                      179, 203, 251, 257, 305, 341, 365, 389, 437, 461, 485,
                                                      515, 587, 619, 691, 739, 751, 799, 847, 895};

  HEGGrid(const PL_t& lat) : Lattice(lat) {}
  ~HEGGrid() = default;


  /** return the estimated number of grid in each direction */
  inline int getNC(int nup) const { return static_cast<int>(std::pow(static_cast<T>(nup), 1.0 / 3.0)) / 2 + 1; }

  //return the number of k-points upto nsh-shell
  inline int getNumberOfKpoints(int nsh) const
  {
    if (nsh < n_within_shell.size())
      return n_within_shell[nsh];
    else
      return -1;
  }

  //return the shell index for nkpt k-points
  inline int getShellIndex(int nkpt) const
  {
    auto loc = std::upper_bound(n_within_shell.begin(), n_within_shell.end(), nkpt);
    if (loc < n_within_shell.end())
      return loc - n_within_shell.begin() - 1;
    else
      return getNC(nkpt);
  }

  /** return the cell size  for the number of particles and rs
   * @param nptcl number of particles
   * @param rs_in rs
   */
  inline T getCellLength(int nptcl, T rs_in) const { return std::pow(4.0 * M_PI * nptcl / 3.0, 1.0 / 3.0) * rs_in; }
};

} // namespace qmcplusplus
#endif
