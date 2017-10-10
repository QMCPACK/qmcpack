//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ION_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_ION_ORBITAL_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

class OrbitalConstraintsBase;

/** IonOrbital IonOrbital Builder with constraints
 */
class IonOrbitalBuilder: public OrbitalBuilderBase
{

public:

  IonOrbitalBuilder(ParticleSet& p, TrialWaveFunction& psi,
                    PtclPoolType& psets);

  bool put(xmlNodePtr cur);

private:
  ///particleset pool to get ParticleSet other than the target
  PtclPoolType& ptclPool;
  ///index for the jastrow type: 1, 2, 3
  int IonOrbitalType;
  ///name
  std::string nameOpt;
  ///type
  std::string typeOpt;
  ///function
  Vector<RealType> widthOpt;
  ///spin
  std::string spinOpt;
  ///transform
  std::string transformOpt;
  ///source
  std::string sourceOpt;
  ///reset the options
  void resetOptions();
};

}
#endif
