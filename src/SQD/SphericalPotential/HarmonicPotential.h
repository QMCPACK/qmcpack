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
#ifndef OHMMSHF_HARMONICPOTENTIAL_H
#define OHMMSHF_HARMONICPOTENTIAL_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf
{

/**
   @ingroup RadialPotential
   @class HarmonicPotential
   @brief implements the Harmonic potential of
   \f[
   V_{Harmonic}(r) = \frac{1}{2}\omega^2 r^2
   \f]
 */
struct HarmonicPotential: public RadialPotentialBase
{
  value_type Omega;
  HarmonicPotential(value_type omega);
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  ///return \f$ n \f$
  int getNumOfNodes(int n, int l);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
