//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBACKFLOW_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBACKFLOW_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include<cmath>

namespace qmcplusplus
{

class SlaterDetWithBackflow: public SlaterDet
{
public:
  BackflowTransformation *BFTrans;

  /**  constructor
   * @param targetPtcl target Particleset
   * @param rn release node
   */
  SlaterDetWithBackflow(ParticleSet& targetPtcl,BackflowTransformation *BF);

  ///destructor
  ~SlaterDetWithBackflow();

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    BFTrans = bf;
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->setBF(BFTrans);
  }

  void resetTargetParticleSet(ParticleSet& P);
  void checkInVariables(opt_variables_type& active)
  {
    //if(Optimizable) {
    if(BFTrans->isOptimizable())
    {
      BFTrans->checkInVariables(active);
      for(int i=0; i<Dets.size(); i++)
        Dets[i]->checkInVariables(active);
    }
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    //if(Optimizable) {
    if(BFTrans->isOptimizable())
    {
      BFTrans->checkOutVariables(active);
      for(int i=0; i<Dets.size(); i++)
        Dets[i]->checkOutVariables(active);
    }
  }

  ///reset all the Dirac determinants, Optimizable is true
  void resetParameters(const opt_variables_type& active)
  {
    //if(Optimizable) {
    if(BFTrans->isOptimizable())
    {
      BFTrans->resetParameters(active);
      for(int i=0; i<Dets.size(); i++)
        Dets[i]->resetParameters(active);
    }
  }

  ValueType evaluate(ParticleSet& P
                     ,ParticleSet::ParticleGradient_t& G
                     ,ParticleSet::ParticleLaplacian_t& L);

  RealType evaluateLog(ParticleSet& P
                       ,ParticleSet::ParticleGradient_t& G
                       ,ParticleSet::ParticleLaplacian_t& L);

  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);
  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
  void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);
  void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

  inline ValueType ratio(ParticleSet& P, int iat,
                         ParticleSet::ParticleGradient_t& dG,
                         ParticleSet::ParticleLaplacian_t& dL)
  {
    BFTrans->evaluatePbyPAll(P,iat);
    //BFTrans->evaluate(P);
    ValueType psi=1.0;
    for(int i=0; i<Dets.size(); ++i)
      psi *= Dets[i]->ratio(P,iat,dG,dL);
    return psi;
  }

  inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    BFTrans->evaluatePbyPWithGrad(P,iat);
    //BFTrans->evaluate(P);
    ValueType psi=1.0;
    for(int i=0; i<Dets.size(); ++i)
      psi *= Dets[i]->ratioGrad(P,iat,grad_iat);
    return psi;
  }

  inline ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::ratioGrad() \n");
    return ValueType();
  }

  GradType evalGrad(ParticleSet& P, int iat)
  {
    QMCTraits::GradType g;
    for(int i=0; i<Dets.size(); ++i)
      g += Dets[i]->evalGrad(P,iat);
    return g;
  }

  GradType alternateEvalGrad(ParticleSet& P, int iat)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::alternateEvalGrad() \n");
    return GradType();
  }

  GradType evalGradSource(ParticleSet& P, ParticleSet &src, int iat)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
    return ValueType();
  }

  GradType evalGradSource (ParticleSet& P, ParticleSet& src, int iat,
                           TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
                           TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
    return GradType();
  }

  inline ValueType logRatio(ParticleSet& P, int iat,
                            ParticleSet::ParticleGradient_t& dG,
                            ParticleSet::ParticleLaplacian_t& dL)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::logRatio() \n");
    return ValueType();
  }

  inline void acceptMove(ParticleSet& P, int iat)
  {
    BFTrans->acceptMove(P,iat);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->acceptMove(P,iat);
  }

  inline void restore(int iat)
  {
    BFTrans->restore(iat);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->restore(iat);
  }


  inline ValueType ratio(ParticleSet& P, int iat)
  {
    BFTrans->evaluatePbyP(P,iat);
    //BFTrans->evaluate(P);
    ValueType ratio=1.0;
    for(int i=0; i<Dets.size(); ++i)
      ratio*=Dets[i]->ratio(P,iat);
    return ratio;
  }

  inline ValueType alternateRatio(ParticleSet& P)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::alternateRatio() \n");
    return ValueType();
  }

  inline void alternateGrad(ParticleSet::ParticleGradient_t& G)
  {
    APP_ABORT("Need to implement SlaterDetWithBackflow::alternateRatio() \n");
  }

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat)
  {
    for(int i=0; i<Dets.size(); ++i)
      Dets[i]->update(P,dG,dL,iat);
  }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  SPOSetBasePtr getPhi(int i=0)
  {
    return Dets[i]->getPhi();
  }

  void get_ratios(ParticleSet& P, vector<ValueType>& ratios);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           vector<RealType>& dlogpsi,
                           vector<RealType>& dhpsioverpsi);

  void testDerivGL(ParticleSet& P);

  //private:
  //SlaterDetWithBackflow() {}
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 4711 $   $Date: 2010-03-07 18:06:54 -0600 (Sun, 07 Mar 2010) $
 * $Id: SlaterDetWithBackflow.h 4711 2010-03-08 00:06:54Z jeongnim.kim $
 ***************************************************************************/
