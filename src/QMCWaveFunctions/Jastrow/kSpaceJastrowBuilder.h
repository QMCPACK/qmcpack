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


#ifndef QMCPLUSPLUS_KSPACE_JASTROW_BUILDER_H
#define QMCPLUSPLUS_KSPACE_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"

namespace qmcplusplus
{
//forward declaration
class ParticleSet;

class kSpaceJastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  const ParticleSet& sourcePtcl;
  std::map<std::string, kSpaceJastrow::SymmetryType> SymmMap;
  // One-body constructor
  kSpaceJastrowBuilder(Communicate* comm, ParticleSet& target, const ParticleSet& source)
      : WaveFunctionComponentBuilder(comm, target), sourcePtcl(source)
  {
    // nothing for now
    SymmMap["crystal"]   = kSpaceJastrow::CRYSTAL;
    SymmMap["isotropic"] = kSpaceJastrow::ISOTROPIC;
    SymmMap["none"]      = kSpaceJastrow::NOSYMM;
  }

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;
  void outputJastrow(kSpaceJastrow* jastrow);
};

} // namespace qmcplusplus
#endif
