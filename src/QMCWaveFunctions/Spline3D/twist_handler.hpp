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
    
    
#ifndef QMCPLUSPLUS_TWIST_HANDLER_H
#define QMCPLUSPLUS_TWIST_HANDLER_H

#include <Lattice/CrystalLattice.h>

namespace qmcplusplus
{

  template<typename T, unsigned D>
    struct twist_handler
    {
      typedef TinyVector<T, D> point_type;
      typedef CrystalLattice<T, D> unitcell_type;
      /** k-vector in the reduced unit */
      point_type TwistVector;
      /** k-vector in the Cartesian unit */
      point_type kVector;
      /** Lattice for the super cell
      */
      unitcell_type SuperLattice;
      /** Lattice for the primitive cell
      */
      unitcell_type PrimLattice;
      /** type matrix */
      Tensor<int, D> TileMatrix;
      /** boolean to unpack the copies */
      std::vector<bool> MakeTwoCopies;
      /** k-points for each orbital */
      std::vector<point_type> kPoints;
      /** \f$k^2\$ */
      std::vector<T> ksq;
      /** character of G-vectors with real orbitals
       *
       * HalfG[i]=0, if twist is zero
       */
      TinyVector<int, D> HalfG;
      /** virtual destructor */
      virtual ~twist_handler() { }
    };
}
#endif
