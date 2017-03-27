//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  std::string nameOpt;
  ///jastrow/@type
  std::string typeOpt;
  ///jastrow/@function
  std::string funcOpt;
  ///jastrow/@spin
  std::string spinOpt;
  ///jastrow/@transform
  std::string transformOpt;
  ///jastrow/@source
  std::string sourceOpt;
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
