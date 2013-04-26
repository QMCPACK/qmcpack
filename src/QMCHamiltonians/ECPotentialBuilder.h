//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ECPPOTENTIAL_BUILDER_H
#define QMCPLUSPLUS_ECPPOTENTIAL_BUILDER_H
#include "Configuration.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
namespace qmcplusplus
{

class QMCHamiltonian;
class ParticleSet;
class TrialWaveFunction;

struct ECPotentialBuilder: public MPIObjectBase, public QMCTraits
{

  typedef LocalECPotential::RadialPotentialType RadialPotentialType;
  typedef LocalECPotential::GridType GridType;
  bool hasLocalPot;
  bool hasNonLocalPot;

  QMCHamiltonian&  targetH;
  ParticleSet& IonConfig;
  ParticleSet& targetPtcl;
  TrialWaveFunction& targetPsi;

  vector<RealType>  localZeff;
  vector<RadialPotentialType*>  localPot;
  vector<NonLocalECPComponent*>  nonLocalPot;

  ECPotentialBuilder(QMCHamiltonian& h,
                     ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi,
                     Communicate* c);

  bool put(xmlNodePtr cur);

  void useSimpleTableFormat();
  void useXmlFormat(xmlNodePtr cur);
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

