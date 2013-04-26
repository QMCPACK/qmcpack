#ifndef LOCAL_ECPOTENTIAL_CUDA_H
#define LOCAL_ECPOTENTIAL_CUDA_H

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

#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/CudaCoulomb.h"

class MCWalkerConfiguration;

namespace qmcplusplus
{
struct LocalECPotential_CUDA : public LocalECPotential
{
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  //// Short-range part
  int NumIons, NumElecs, NumElecGroups, NumIonSpecies;
  ParticleSet &ElecRef, &IonRef;
  vector<int> IonFirst, IonLast;
  // This is indexed by the ion species
  vector<TextureSpline*> SRSplines;
  TextureSpline *V0Spline;
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::host_vector<CUDA_PRECISION>    SumHost;
  gpu::device_vector<CUDA_PRECISION>  IGPU;
  gpu::device_vector<CUDA_PRECISION>  ZionGPU;

  vector<PosType> SortedIons;
  void add(int groupID, RadialPotentialType* ppot, RealType zion);

  void addEnergy(MCWalkerConfiguration &W,
                 vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  LocalECPotential_CUDA(ParticleSet& ions, ParticleSet& elns);

};
}

#endif
