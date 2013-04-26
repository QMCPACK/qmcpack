//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
