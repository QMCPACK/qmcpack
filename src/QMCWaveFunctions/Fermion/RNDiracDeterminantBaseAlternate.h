//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);

  void resize(int nel, int morb);
  void restore(int iat);

  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  ValueType ratio(ParticleSet& P, int iat);

  ValueType alternateRatio(ParticleSet& P);

  /** return the ratio
   * @param P current configuration
   * @param iat particle whose position is moved
   * @param dG differential Gradients
   * @param dL differential Laplacians
   *
   * Data member *_temp contain the data assuming that the move is accepted
   * and are used to evaluate differential Gradients and Laplacians.
   */
  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL);

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

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);

  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);
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
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4473 $   $Date: 2009-12-08 11:38:31 -0600 (Tue, 08 Dec 2009) $
 * $Id: RNDiracDeterminantBaseAlternate.h 4473 2009-12-08 17:38:31Z jmcminis $
 ***************************************************************************/
