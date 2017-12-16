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
    
    


#ifndef OHMMSHF_SJPSEUDOPOTENTIAL_H
#define OHMMSHF_SJPSEUDOPOTENTIAL_H
#include "SQD/SphericalPotential/RadialPotential.h"
#include "Optimize/VarList.h"
namespace ohmmshf
{

/**
   @ingroup RadialPotential
   @class SJPseudoPotential
   @brief implements the Starkloff-Joannopoulos pseudopotential
   *
   \f[
   V_{SJ}(r) = -\frac{Z_{Eff}}{r}\frac{1-e^{-\lambda r}}
   {1+e^{-\lambda (r-r_c)}}
   \f]
   *
   See Th. Starkloff and J.D. Joannopoulos, Phys. Rev. B, \textbf{16},
   5212, (1977).
*/
struct SJPseudoPotential: public RadialPotentialBase
{
  value_type Zeff, SJ_lambda, rc;
  SJPseudoPotential(VarRegistry<value_type>&,value_type,
                    value_type,value_type);
  SJPseudoPotential(value_type,value_type,value_type);
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t&V, int norb);
  ///must always return 1
  inline int getNumOfNodes(int n, int l);
};
}
#endif
