//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file RNDiracDeterminantBaseAlternateBase.h
 * @brief Declaration of RNDiracDeterminantBaseAlternate with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_RNDIRACDETERMINANTWITHBASEALTERNATE_H
#define QMCPLUSPLUS_RNDIRACDETERMINANTWITHBASEALTERNATE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
// #include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "DiracDeterminantBase.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

class RNDiracDeterminantBaseAlternate: public DiracDeterminantBase
{
public:
  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  RNDiracDeterminantBaseAlternate(SPOSetBasePtr const &spos, int first=0);

  ///default destructor
  ~RNDiracDeterminantBaseAlternate();

  /**copy constructor
   * @param s existing RNDiracDeterminantBaseAlternate
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  RNDiracDeterminantBaseAlternate(const RNDiracDeterminantBaseAlternate& s);

  RNDiracDeterminantBaseAlternate& operator=(const RNDiracDeterminantBaseAlternate& s);

  void registerData(ParticleSet& P, WFBufferType& buf);

  void resize(int nel, int morb);
  void restore(int iat);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  ValueType ratio(ParticleSet& P, int iat);

  ValueType alternateRatio(ParticleSet& P);

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  GradType evalGrad(ParticleSet& P, int iat);
  GradType alternateEvalGrad(ParticleSet& P, int iat);
  GradType evalGradSource(ParticleSet &P, ParticleSet &source,
                          int iat);

  GradType evalGradSource
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  GradType evalGradSourcep
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat);

  RealType getAlternatePhaseDiff()
  {
    return evaluatePhase(curRatio);
  }
  RealType getAlternatePhaseDiff(int iat)
  {
    return evaluatePhase(curRatio);
  }

  ///evaluate log of determinant for a particle set: should not be called
  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L) ;

  DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;


  ValueType alternateCurRatio;
  //    double ComputeExtraTerms(int ptcl_gradient, int elDim,int ionDim);
  ParticleSet::ParticleGradient_t myG_alternate;
  ParticleSet::ParticleLaplacian_t myL_alternate;
  ValueType logepsilon;
  ValueType alternateLogValue;
  ValueType alternatePhaseValue;
//
  inline void setLogEpsilon(ValueType x)
  {
    logepsilon=x;
  }

};



}
#endif
