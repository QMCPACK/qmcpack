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
#ifndef OHMMSHF_STEPPOTENTIAL_H
#define OHMMSHF_STEPPOTENTIAL_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf {

  /**
     @ingroup RadialPotential
     @class StepPotential
     @brief implements the multiple step potential profile.
     *
     *\image html steppot.png
   */
  struct StepPotential: public RadialPotentialBase {
    vector<value_type> Rseg;
    vector<value_type> Vseg;
    StepPotential();
    value_type evaluate(const BasisSetType& psi, 
			RadialOrbitalSet_t& V, int norb);
    ///return \f$ n \f$
    int getNumOfNodes(int n, int l);

    bool put(xmlNodePtr cur);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
