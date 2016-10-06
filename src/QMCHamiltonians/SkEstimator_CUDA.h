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
    
    
#ifndef QMCPLUSPLUS_SK_ESTIMATOR_CUDA_H
#define QMCPLUSPLUS_SK_ESTIMATOR_CUDA_H
#include "QMCHamiltonians/SkEstimator.h"

namespace qmcplusplus
{
class SkEstimator_CUDA : public SkEstimator
{
public:
  SkEstimator_CUDA(ParticleSet& elns) : SkEstimator(elns) {}
  void addEnergy(MCWalkerConfiguration &W,  std::vector<RealType> &LocalEnergy);
};
}

#endif
