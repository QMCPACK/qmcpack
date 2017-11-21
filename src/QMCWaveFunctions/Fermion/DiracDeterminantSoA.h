//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantSoA.h
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_SOA_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_SOA_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"

namespace qmcplusplus
{

  /** a DiracDeterminantBase which uses SoA for G & L and implements both memory and compute optimization*/
class DiracDeterminantSoA: public DiracDeterminantBase
{

  public:
  
  typedef VectorSoaContainer<ValueType,DIM+1> GLVector_t;
  typedef VectorSoaContainer<ValueType,DIM+2> VGLVector_t;

  bool ComputeDeterminant;
  /** full GL container to compute gradient and laplacian
   * 
   * mGL[i][j] returns (gx,gy,gz,lap) for Phi(i,j) SPO
   */
  aligned_vector<GLVector_t> mGL;
  /** a vector of (gx,gy,gz,lap) for the active particle */
  VGLVector_t vVGL; 
  /** memory management for this object */
  aligned_vector<ValueType> memoryPool;
  /** size of Norb for alignment */
  size_t NorbPad;
  /** NorbPad*(OHMMS_DIM+1) for the internal size */
  size_t BlockSize;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantSoA(SPOSetBasePtr const &spos, int first=0);

  ///final destructure
  ~DiracDeterminantSoA();

  /**copy constructor
   * @param s existing DiracDeterminantSoA
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantSoA(const DiracDeterminantSoA& s);

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);
  ///cloning with spo
  DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;

  RealType evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L) ;

  /** recompute evaluate everything clean */
  void recompute(ParticleSet& P);

  void updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L);

  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratio(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  void acceptMove(ParticleSet& P, int iat);

  void registerData(ParticleSet& P, WFBufferType& buf);
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

#if 0
  GradType evalGradSource(ParticleSet &P, ParticleSet &source, int iat);

  GradType evalGradSource
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  GradType evalGradSourcep
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

#endif


};

}
#endif
