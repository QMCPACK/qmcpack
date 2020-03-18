//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ECPPOTENTIAL_BUILDER_H
#define QMCPLUSPLUS_ECPPOTENTIAL_BUILDER_H
#include "Configuration.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "QMCHamiltonians/SOECPotential.h"
#include "QMCHamiltonians/L2Potential.h"
namespace qmcplusplus
{
class QMCHamiltonian;
class ParticleSet;
class TrialWaveFunction;

struct ECPotentialBuilder : public MPIObjectBase, public QMCTraits
{
  typedef LocalECPotential::RadialPotentialType RadialPotentialType;
  typedef LocalECPotential::GridType GridType;
  bool hasLocalPot;
  bool hasNonLocalPot;
  bool hasSOPot;
  bool hasL2Pot;

  QMCHamiltonian& targetH;
  ParticleSet& IonConfig;
  ParticleSet& targetPtcl;
  TrialWaveFunction& targetPsi;

  std::vector<RealType> localZeff;
  std::vector<RadialPotentialType*> localPot;
  std::vector<NonLocalECPComponent*> nonLocalPot;
  std::vector<SOECPComponent*> soPot;
  std::vector<L2RadialPotential*> L2Pot;

  ECPotentialBuilder(QMCHamiltonian& h, ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi, Communicate* c);

  bool put(xmlNodePtr cur);

  void useSimpleTableFormat();
  void useXmlFormat(xmlNodePtr cur);
};
} // namespace qmcplusplus
#endif
