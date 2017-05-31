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
    
    
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/Jastrow/JAAPBCBuilder.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus
{

/** constructor
 * @param p target ParticleSet whose wave function is to be initialized
 *@param psi the wavefunction
 *@param psets a vector containing the ParticleSets
 *
 *Jastrow wavefunctions chose the needed distance tables and the
 *DistanceTableData objects are initialized based on the source
 *and target particle sets.
 */
JAAPBCBuilder::JAAPBCBuilder(ParticleSet& p, TrialWaveFunction& psi):
  OrbitalBuilderBase(p,psi)
{ }

bool JAAPBCBuilder::put(xmlNodePtr cur)
{
  typedef LRJastrowSingleton::LRHandlerType HandlerType;
  HandlerType* handler = LRJastrowSingleton::getHandler(targetPtcl);
  LRTwoBodyJastrow *J2 = new LRTwoBodyJastrow(targetPtcl, handler);
  bool success = J2->put(cur, targetPsi.VarList);
  J2->setOptimizable(true);
  targetPsi.addOrbital(J2);
  return success;
}
}
