//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMSHF_ZOVERRFUNCTOR_H
#define OHMMSHF_ZOVERRFUNCTOR_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf {

  /**
     @ingroup RadialPotential
     @class ZOverRPotential
     @brief implements the Nuclear potential of 
     \f[ 
     V_{Nuclear}(r) = -\frac{Z}{r} 
     \f]
  */
  struct ZOverRPotential: public RadialPotentialBase {
    value_type Z;
    ZOverRPotential(value_type z);
    value_type evaluate(const BasisSetType& psi, 
			RadialOrbitalSet_t& V, int norb);
    ///return \f$ n-l-1 \f$
    int getNumOfNodes(int n, int l);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

  
