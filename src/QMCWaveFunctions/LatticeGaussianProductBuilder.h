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


#ifndef QMCPLUSPLUS_LATTICE_GAUSSIAN_BUILDER_H
#define QMCPLUSPLUS_LATTICE_GAUSSIAN_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
/** LatticeGaussianProduct LatticeGaussianProduct Builder with constraints
 */
class LatticeGaussianProductBuilder : public WaveFunctionComponentBuilder
{
public:
  LatticeGaussianProductBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets);

  std::unique_ptr<WaveFunctionComponent> buildComponent(xmlNodePtr cur) override;

private:
  ///particleset pool to get ParticleSet other than the target
  const PtclPoolType& ptclPool;
  ///index for the jastrow type: 1, 2, 3
  int LatticeGaussianProductType;
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

} // namespace qmcplusplus
#endif
