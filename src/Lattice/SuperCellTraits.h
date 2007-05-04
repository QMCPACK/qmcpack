//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SUPERCELL_TRAITS_H
#define QMCPLUSPLUS_SUPERCELL_TRAITS_H
#include "OhmmsPETE/TinyVector.h"
namespace qmcplusplus
{

  /** enumeration to clssify a CrystalLattice
  */
  enum {SUPERCELL_OPEN=0, SUPERCELL_WIRE=1, SUPERCELL_SLAB=3, SUPERCELL_BULK=7};

  /** dummy class to detemine the supercell type **/
  template<unsigned D>
    struct SuperCellType {};

  /** specialization of SuperCellType for 3-dimensional cell
  */
  template<>
    struct SuperCellType<3> {
      /** convert box to an integer
       * @param box 3-dimensional boolean vector
       *
       * When all the directions are open, returns 0.
       * - 7 (1,1,1) bulk system
       * - 3 (1,1,0) slab system
       * - 1 (1,0,0) wire system
       */
      inline static int apply(const TinyVector<int,3>& box) {
        return box[0]+2*(box[1]+box[2]*2);
      }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1189 $   $Date: 2006-07-17 10:08:11 -0500 (Mon, 17 Jul 2006) $
 * $Id: DistanceTable.cpp 1189 2006-07-17 15:08:11Z jnkim $ 
 ***************************************************************************/
