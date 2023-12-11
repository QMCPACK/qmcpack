//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EEI_JASTROW_BUILDER_H
#define QMCPLUSPLUS_EEI_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
//forward declaration
class ParticleSet;

class eeI_JastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  ParticleSet* sourcePtcl;
  // Two-body constructor
  eeI_JastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  template<typename J3type>
  bool putkids(xmlNodePtr kids, J3type& J3);
};

} // namespace qmcplusplus
#endif
