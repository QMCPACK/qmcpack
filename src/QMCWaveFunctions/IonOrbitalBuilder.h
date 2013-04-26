//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Ken Esler
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
  string nameOpt;
  ///type
  string typeOpt;
  ///function
  Vector<RealType> widthOpt;
  ///spin
  string spinOpt;
  ///transform
  string transformOpt;
  ///source
  string sourceOpt;
  ///reset the options
  void resetOptions();
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: esler $
 * $Revision: 1666 $   $Date: 2007-01-30 13:00:15 -0600 (Tue, 30 Jan 2007) $
 * $Id: IonOrbitalBuilder.h 1666 2007-01-30 19:00:15Z jnkim $
 ***************************************************************************/
