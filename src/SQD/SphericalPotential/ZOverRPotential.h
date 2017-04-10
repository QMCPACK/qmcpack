//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMSHF_ZOVERRFUNCTOR_H
#define OHMMSHF_ZOVERRFUNCTOR_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf
{

/**
   @ingroup RadialPotential
   @class ZOverRPotential
   @brief implements the Nuclear potential of
   \f[
   V_{Nuclear}(r) = -\frac{Z}{r}
   \f]
*/
struct ZOverRPotential: public RadialPotentialBase
{
  value_type Z;
  ZOverRPotential(value_type z);
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  ///return \f$ n-l-1 \f$
  int getNumOfNodes(int n, int l);
};
}
#endif


