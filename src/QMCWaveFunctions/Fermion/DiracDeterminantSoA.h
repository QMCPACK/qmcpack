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
  
  typedef OrbitalSetTraits<ValueType>::GLVector_t       GLVector_t;

  bool ComputeDeterminant;
  /** full GL container to compute gradient and laplacian
   * 
   * mGL[i][j] returns (gx,gy,gz,lap) for Phi(i,j) SPO
   */
  std::vector<GLVector_t> mGL;
  /** a vector of (gx,gy,gz,lap) for the active particle */
  GLVector_t vGL; 
  /** memory management for this object */
  aligned_vector<ValueType> memoryPool;
  /** size of Norb for alignment */
  size_t NorbPad;
  /** NorbPad*(OHMMS_DIM+1) for the internal size */
  size_t BlockSize;

  /** handle inverse, update */
  DiracMatrix<ValueType> detEng;
#ifdef MIXED_PRECISION
  DiracMatrix<mValueType> detEng_hp;
#endif

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

  void recompute(ParticleSet& P);


  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);

  void updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L);
 
  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratio(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat);

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat);


#if 0
  void update(ParticleSet& P,
      ParticleSet::ParticleGradient_t& dG,
      ParticleSet::ParticleLaplacian_t& dL,
      int iat);

  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  GradType evalGrad(ParticleSet& P, int iat);
  GradType evalGradSource(ParticleSet &P, ParticleSet &source, int iat);

  GradType evalGradSource
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  GradType evalGradSourcep
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  virtual ValueType logRatio(ParticleSet& P, int iat,
                             ParticleSet::ParticleGradient_t& dG,
                             ParticleSet::ParticleLaplacian_t& dL);
#endif


};

}
#endif
