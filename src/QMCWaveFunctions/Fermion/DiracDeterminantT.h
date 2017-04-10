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
    
    
/** @file DiracDeterminantT.h
 * @brief declaration of DiracDeterminantT
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_TEMPLATE_H
#define QMCPLUSPLUS_DIRACDETERMINANT_TEMPLATE_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
namespace qmcplusplus
{
template <class SPOSet>
class DiracDeterminantT: public DiracDeterminantBase
{
  SPOSet& Phi;

public:

  DiracDeterminantT(SPOSet& spos, int first =0): DiracDeterminantBase(first), Phi(spos) { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    Phi.resetTargetParticleSet(P);
  }
  void reset()
  {
    Phi.reset();
  }

  void evaluateSingle(ParticleSet& P, int iat)
  {
    Phi.evaluate(P,iat,psiV);
  }
  void evaluateSingleAll(ParticleSet& P, int iat)
  {
    Phi.evaluate(P, iat, psiV, dpsiV, d2psiV);
  }
  void evaluateAll(ParticleSet& P)
  {
    Phi.evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  }
};
}
#endif
