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

#include "QMCHamiltonians/CoulombPotential_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{
CoulombPotentialAA_CUDA::CoulombPotentialAA_CUDA(ParticleSet* s, bool quantum):
  CoulombPotential<OHMMS_PRECISION>(s,0,quantum),
  SumGPU("CoulombPotentialAA_CUDA::SumGPU")
{
  NumElecs = s->getTotalNum();
}

void
CoulombPotentialAA_CUDA::addEnergy(MCWalkerConfiguration &W,
                                   vector<RealType> &LocalEnergy)
{
  vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  if (SumGPU.size() < nw)
    SumGPU.resize(nw);
  // Evaluate sum on GPU
  CoulombAA_Sum(W.RList_GPU.data(), NumElecs, SumGPU.data(), nw);
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] =
      SumHost[iw];
    LocalEnergy[iw] += SumHost[iw];
  }
}

QMCHamiltonianBase*
CoulombPotentialAA_CUDA::makeClone (ParticleSet& qp, TrialWaveFunction& psi)
{
  return new CoulombPotentialAA_CUDA(&qp,true);
}


CoulombPotentialAB_CUDA::CoulombPotentialAB_CUDA(ParticleSet* s, ParticleSet* t)
  : CoulombPotential<OHMMS_PRECISION>(s,t,true),
    SumGPU("CoulombPotentialAB_CUDA::SumGPU"),
    IGPU("CoulombPotentialAB_CUDA::IGPU"),
    ZionGPU("CoulombPotentialAB_CUDA::ZionGPOU")
{
  SpeciesSet &sSet = s->getSpeciesSet();
  NumIonSpecies = sSet.getTotalNum();
  NumIons  = s->getTotalNum();
  NumElecs = t->getTotalNum();
  // Copy center positions to GPU, sorting by GroupID
  gpu::host_vector<CUDA_PRECISION> I_host(OHMMS_DIM*NumIons);
  int index=0;
  for (int cgroup=0; cgroup<NumIonSpecies; cgroup++)
  {
    IonFirst.push_back(index);
    for (int i=0; i<NumIons; i++)
    {
      if (s->GroupID[i] == cgroup)
      {
        for (int dim=0; dim<OHMMS_DIM; dim++)
          I_host[OHMMS_DIM*index+dim] = s->R[i][dim];
        SortedIons.push_back(s->R[i]);
        index++;
      }
    }
    IonLast.push_back(index-1);
  }
  IGPU = I_host;
}

void
CoulombPotentialAB_CUDA::addEnergy(MCWalkerConfiguration &W,
                                   vector<RealType> &LocalEnergy)
{
  vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  if (SumGPU.size() < nw)
    SumGPU.resize(nw);
  // Evaluate sum on GPU
  CoulombAB_Sum (W.RList_GPU.data(), NumElecs,
                 IGPU.data(), ZionGPU.data(), NumIons, SumGPU.data(), nw);
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = SumHost[iw];
    LocalEnergy[iw] += SumHost[iw];
  }
}

QMCHamiltonianBase*
CoulombPotentialAB_CUDA::makeClone (ParticleSet& qp, TrialWaveFunction& psi)
{
  return new CoulombPotentialAB_CUDA(Pa,&qp);
}

}
