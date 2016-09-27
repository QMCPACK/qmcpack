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
                 std::vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

/** CoulombPotential for ion-el
 *
 * Derived from CoulombPotential<T>
 */
struct CoulombPotentialAB_CUDA : public CoulombPotential<OHMMS_PRECISION>
{
  int NumElecs, NumIons, NumIonSpecies;
  std::vector<int> IonFirst, IonLast;
  std::vector<PosType> SortedIons;

  CoulombPotentialAB_CUDA(ParticleSet* s, ParticleSet* t);

  gpu::host_vector<CUDA_PRECISION>   SumHost;
  gpu::device_vector<CUDA_PRECISION>  SumGPU;
  gpu::device_vector<CUDA_PRECISION> ZionGPU;
  gpu::device_vector<CUDA_PRECISION>    IGPU;

  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
};

}
#endif
