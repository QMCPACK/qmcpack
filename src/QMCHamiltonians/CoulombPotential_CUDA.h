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


namespace qmcplusplus
{

class MCWalkerConfiguration;

/** CoulombPotential for el-el
 *
 * Derived from CoulombPotential<T>
 */
struct CoulombPotentialAA_CUDA : public CoulombPotential<OHMMS_PRECISION>
{
  int NumElecs;
  CoulombPotentialAA_CUDA(ParticleSet* s, bool quantum);

  gpu::device_vector<CUDA_PRECISION> SumGPU;
  gpu::host_vector<CUDA_PRECISION>  SumHost;
  void addEnergy(MCWalkerConfiguration &W,
                 vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

/** CoulombPotential for ion-el
 *
 * Derived from CoulombPotential<T>
 */
struct CoulombPotentialAB_CUDA : public CoulombPotential<OHMMS_PRECISION>
{
  int NumElecs, NumIons, NumIonSpecies;
  vector<int> IonFirst, IonLast;
  vector<PosType> SortedIons;

  CoulombPotentialAB_CUDA(ParticleSet* s, ParticleSet* t);

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
