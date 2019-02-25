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

namespace qmcplusplus
{

ACForce::ACForce(ParticleSet& source, ParticleSet& target, TrialWaveFunction& psi, QMCHamiltonian& H)
{};

ACForce::Return_t ACForce::evaluate(ParticleSet& P)
{
  return 0.0;
};  

void ACForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
};
void ACForce::setObservables(PropertySetType& plist)
{
};
void ACForce::setParticlePropertyList(PropertySetType& plist, int offset)
{
};

}
