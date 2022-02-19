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
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
class OrbitalConstraintsBase;

/** Jastrow Jastrow Builder with constraints
 */
class JastrowBuilder : public WaveFunctionComponentBuilder
{
public:
  JastrowBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  ///particleset pool to get ParticleSet other than the target
  const PtclPoolType& ptclPool;
  ///index for the jastrow type: 1, 2, 3
  int JastrowType;
  /// \xmla{jastrow,name}
  std::string nameOpt;
  /// \xmla{jastrow,type}
  std::string typeOpt;
  /// \xmla{jastrow,function}
  std::string funcOpt;
  /// \xmla{jastrow,spin}
  std::string spinOpt;
  /// \xmla{jastrow,transform}
  std::string transformOpt;
  /// \xmla{jastrow,source}
  std::string sourceOpt;
  ///reset the options
  void resetOptions();
  /// build one-body term
  std::unique_ptr<WaveFunctionComponent> buildOneBody(xmlNodePtr cur);
  /// build two-body term
  std::unique_ptr<WaveFunctionComponent> buildTwoBody(xmlNodePtr cur);
  /// build electron-electron ion term
  std::unique_ptr<WaveFunctionComponent> build_eeI(xmlNodePtr cur);
  /// build k-Space term
  std::unique_ptr<WaveFunctionComponent> buildkSpace(xmlNodePtr cur);
  /// build number-counting term
  std::unique_ptr<WaveFunctionComponent> buildCounting(xmlNodePtr cur);
};

} // namespace qmcplusplus
#endif
