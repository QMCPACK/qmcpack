//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_GENERALIZED_JASTROWBUILDER_H
#define QMCPLUSPLUS_GENERALIZED_JASTROWBUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

class OrbitalConstraintsBase;

/** Jastrow Jastrow Builder with constraints
 */
class JastrowBuilder: public OrbitalBuilderBase
{

public:

  JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets);

  bool put(xmlNodePtr cur);

private:
  ///particleset pool to get ParticleSet other than the target
  PtclPoolType& ptclPool;
  ///index for the jastrow type: 1, 2, 3
  int JastrowType;
  ///jastrow/@name
  string nameOpt;
  ///jastrow/@type
  string typeOpt;
  ///jastrow/@function
  string funcOpt;
  ///jastrow/@spin
  string spinOpt;
  ///jastrow/@transform
  string transformOpt;
  ///jastrow/@source
  string sourceOpt;
  ///reset the options
  void resetOptions();
  ///add one-body term
  bool addOneBody(xmlNodePtr cur);
  ///add two-body term
  bool addTwoBody(xmlNodePtr cur);
  ///add three-body term
  bool addThreeBody(xmlNodePtr cur);
  /// add electron-electron ion term
  bool add_eeI (xmlNodePtr cur);
  ///add k-Space term
  bool addkSpace(xmlNodePtr cur);
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1666 $   $Date: 2007-01-30 13:00:15 -0600 (Tue, 30 Jan 2007) $
 * $Id: JastrowBuilder.h 1666 2007-01-30 19:00:15Z jnkim $
 ***************************************************************************/
