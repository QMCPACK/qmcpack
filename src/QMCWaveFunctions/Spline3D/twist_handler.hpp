//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: twist_handler.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
