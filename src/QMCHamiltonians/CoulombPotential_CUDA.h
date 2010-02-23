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
#ifndef QMCPLUSPLUS_COULOMB_POTENTIAL_CUDA_H
#define QMCPLUSPLUS_COULOMB_POTENTIAL_CUDA_H
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/CudaCoulomb.h"

class MCWalkerConfiguration;

namespace qmcplusplus {
  struct CoulombPotentialAA_CUDA : public CoulombPotentialAA
  {
    int NumElecs;
    //////////////////////////////////
    // Vectorized evaluation on GPU //
    //////////////////////////////////
    CoulombPotentialAA_CUDA(ParticleSet& ref);

    ParticleSet &PtclRef;
    gpu::device_vector<CUDA_PRECISION> SumGPU;
    gpu::host_vector<CUDA_PRECISION>  SumHost;
    void addEnergy(MCWalkerConfiguration &W, 
		   vector<RealType> &LocalEnergy);

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  };

  struct CoulombPotentialAB_CUDA : public CoulombPotentialAB
  {
    ParticleSet &PtclRef;
    int NumElecs, NumIons, NumIonSpecies;
    vector<int> IonFirst, IonLast;
    vector<PosType> SortedIons;

    //////////////////////////////////
    // Vectorized evaluation on GPU //
    //////////////////////////////////
    CoulombPotentialAB_CUDA(ParticleSet& ref, ParticleSet &ions);

    gpu::host_vector<CUDA_PRECISION>   SumHost;
    gpu::device_vector<CUDA_PRECISION>  SumGPU;
    gpu::device_vector<CUDA_PRECISION> ZionGPU;
    gpu::device_vector<CUDA_PRECISION>    IGPU;

    void addEnergy(MCWalkerConfiguration &W, 
		   vector<RealType> &LocalEnergy);

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  };

}
#endif
