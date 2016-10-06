//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_NONLOCALECPOTENTIAL_CUDA_H
#define QMCPLUSPLUS_NONLOCALECPOTENTIAL_CUDA_H

#include "NonLocalECPotential.h"

namespace qmcplusplus
{

class NonLocalECPotential_CUDA: public NonLocalECPotential
{
protected:
  //////////////////////////////////
  // Vectorized evaluation on GPU //
  //////////////////////////////////
  bool UsePBC;
  int NumIonGroups;
  std::vector<int> IonFirst, IonLast;
  gpu::device_vector<CUDA_PRECISION> Ions_GPU, L, Linv;
  gpu::device_vector<int> Elecs_GPU;
  gpu::host_vector<int> Elecs_host;
  gpu::device_vector<CUDA_PRECISION> Dist_GPU;
  gpu::host_vector<CUDA_PRECISION> Dist_host;
  gpu::device_vector<int*> Eleclist_GPU;
  gpu::device_vector<CUDA_PRECISION*> Distlist_GPU;
  gpu::device_vector<int> NumPairs_GPU;
  gpu::host_vector<int> NumPairs_host;
  gpu::host_vector<int*> Eleclist_host;
  gpu::host_vector<CUDA_PRECISION*> Distlist_host;
  gpu::host_vector<CUDA_PRECISION*> RatioPoslist_host;
  gpu::host_vector<CUDA_PRECISION*> Ratiolist_host;
  gpu::host_vector<CUDA_PRECISION*> CosThetalist_host;

  int NumElecs;
  // The maximum number of quadrature points over all the ions species
  int MaxKnots, MaxPairs, RatiosPerWalker;
  // These are the positions at which we have to evalate the WF ratios
  // It has size OHMMS_DIM * MaxPairs * MaxKnots * NumWalkers
  gpu::device_vector<CUDA_PRECISION> RatioPos_GPU, CosTheta_GPU;
  gpu::host_vector<CUDA_PRECISION> RatioPos_host, CosTheta_host;
  gpu::device_vector<CUDA_PRECISION*> RatioPoslist_GPU, CosThetalist_GPU;

  // Quadrature points
  std::vector<gpu::device_vector<CUDA_PRECISION> > QuadPoints_GPU;
  std::vector<std::vector<CUDA_PRECISION> > QuadPoints_host;
  int CurrentNumWalkers;

  // These are used in calling Psi->NLratios
  std::vector<NLjob> JobList;
  std::vector<PosType> QuadPosList;
  std::vector<ValueType> RatioList;


  std::vector<PosType> SortedIons;

  void setupCUDA(ParticleSet &elecs);
  void resizeCUDA(int nw);

public:
  NonLocalECPotential_CUDA(ParticleSet& ions, ParticleSet& els,
                           TrialWaveFunction& psi, bool usePBC,
                           bool doForces=false);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy);
  void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy,
                 std::vector<std::vector<NonLocalData> > &Txy);
};


}

#endif
