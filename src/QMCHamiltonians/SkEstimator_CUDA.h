//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SK_ESTIMATOR_CUDA_H
#define QMCPLUSPLUS_SK_ESTIMATOR_CUDA_H
#include "QMCHamiltonians/SkEstimator.h"

namespace qmcplusplus
{
class SkEstimator_CUDA : public SkEstimator
{
public:
  SkEstimator_CUDA(ParticleSet& elns) : SkEstimator(elns) {}
  void addEnergy(MCWalkerConfiguration &W,  vector<RealType> &LocalEnergy);
};
}

#endif
