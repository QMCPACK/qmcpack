//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file ACForce.cpp
 *@brief Implementation of ACForce, Assaraf-Caffarel ZVZB style force estimation.
 */
#include "QMCHamiltonians/ACForce.h"
#include <sstream>

namespace qmcplusplus
{

ACForce::ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi_in, QMCHamiltonian& H):
FirstForceIndex(-1), ions(source), elns(target), psi(psi_in), Nions(0)
{
  prefix="ACForce";
  Nions=ions.getTotalNum();
};

ACForce::Return_t ACForce::evaluate(ParticleSet& P)
{
  return 0.0;
};  

void ACForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if(FirstForceIndex<0)
    FirstForceIndex=plist.size();
  for(int iat=0; iat<Nions; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      std::ostringstream obsName;
      obsName << prefix << "_" << iat << "_" << x;
      plist.add(obsName.str());
    }
  }
};
void ACForce::setObservables(PropertySetType& plist)
{
};
void ACForce::setParticlePropertyList(PropertySetType& plist, int offset)
{
};

}
