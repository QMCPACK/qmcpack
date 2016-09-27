//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_IOXML_H
#define QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_IOXML_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.h"

namespace qmcplusplus
{

/** OrbitalBuilder with HFAtomicSTOSet
 */
class HFAtomicSTOSetBuilder: public OrbitalBuilderBase
{

private:

  typedef HFAtomicSTOSet::RadialOrbital_t RadialOrbital_t;
  typedef HFAtomicSTOSet::SPO_t           SPO_t;
  DistanceTableData* d_table;

  /// Maximum Angular Momentum
  int Lmax;

  /// mapping function for Rnl[ RnlID[name] ]
  std::map<std::string,int> RnlID;

  ///temporary storage of radial functions
  std::vector<RadialOrbital_t*> Rnl;

  ///single-particle wave functions
  std::map<std::string, SPO_t*>   OrbSet;

  bool getBasis(xmlNodePtr cur);
  HFAtomicSTOSet* getOrbital(xmlNodePtr cur);

public:

  HFAtomicSTOSetBuilder(ParticleSet& els, TrialWaveFunction& psi, ParticleSet& ions);
  bool put(xmlNodePtr);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
