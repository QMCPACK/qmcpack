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
#ifndef OHMMSHF_SPHERICAL_HEG_POTENTIAL_H
#define OHMMSHF_SPHERICAL_HEG_POTENTIAL_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf
{

/**
   @ingroup RadialPotential
   @class SHEGPotential
   @brief implements spherical jellium
 */
struct SHEGPotential: public RadialPotentialBase
{
  ///number of particles
  value_type Ntot;
  ///external background charge
  value_type Zext;
  value_type Zgauss;
  value_type Sgauss;
  ///density
  value_type Rs;
  ///maximum radius of the background charge
  value_type Rmax;

  SHEGPotential(int nel, value_type rs=1);
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  ///return \f$ n \f$
  int getNumOfNodes(int n, int l);

  bool put(xmlNodePtr cur);
};


}
#endif
/***************************************************************************
 * $RCSfile$   $Author: qmc $
 * $Revision: 42 $   $Date: 2004-09-15 18:05:20 -0500 (Wed, 15 Sep 2004) $
 * $Id: SHEGPotential.h 42 2004-09-15 23:05:20Z qmc $
 ***************************************************************************/
