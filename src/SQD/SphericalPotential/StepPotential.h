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
    
    


#ifndef OHMMSHF_STEPPOTENTIAL_H
#define OHMMSHF_STEPPOTENTIAL_H

#include "SQD/SphericalPotential/RadialPotential.h"

namespace ohmmshf
{

/**
   @ingroup RadialPotential
   @class StepPotential
   @brief implements the multiple step potential profile.
   *
   *\image html steppot.png
 */
struct StepPotential: public RadialPotentialBase
{
  std::vector<value_type> Rseg;
  std::vector<value_type> Vseg;
  StepPotential();
  value_type evaluate(const BasisSetType& psi,
                      RadialOrbitalSet_t& V, int norb);
  ///return \f$ n \f$
  int getNumOfNodes(int n, int l);

  bool put(xmlNodePtr cur);
};
}
#endif
