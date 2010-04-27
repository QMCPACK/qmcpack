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

    void checkInVariables(opt_variables_type& active)
    {
      if(Optimizable) {
        BFTrans->checkInVariables(active);
        for(int i=0; i<Dets.size(); i++) Dets[i]->checkInVariables(active);
      }
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      if(Optimizable) {
        BFTrans->checkOutVariables(active);
        for(int i=0; i<Dets.size(); i++) Dets[i]->checkOutVariables(active);
      }
    }

    ///reset all the Dirac determinants, Optimizable is true
    void resetParameters(const opt_variables_type& active)
    {
      if(Optimizable) {
        BFTrans->resetParameters(active);
        for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
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
      dG=0;
      dL=0;   
      ValueType det = evaluate(P,dG,dL);
      dG-=P.G; 
      dL-=P.L; 
      return det/std::exp(P.Properties(LOGPSI));
    }

    inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::ratioGrad() \n");
      return 0.0;
    }

    inline ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::ratioGrad() \n");
      return 0.0;
    }

    GradType evalGrad(ParticleSet& P, int iat)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::evalGrad() \n");
      return 0.0;
    }

    GradType alternateEvalGrad(ParticleSet& P, int iat)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::alternateEvalGrad() \n");
      return 0.0;
    }

    GradType evalGradSource(ParticleSet& P, ParticleSet &src, int iat)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
      return 0.0;
    }

    GradType evalGradSource (ParticleSet& P, ParticleSet& src, int iat,
        TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
        TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad) 
    {      
      APP_ABORT("Need to implement SlaterDetWithBackflow::evalGradSource() \n");
      return 0.0;
    }

    inline ValueType logRatio(ParticleSet& P, int iat,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::logRatio() \n");
      return 0.0;
    }

    inline void acceptMove(ParticleSet& P, int iat)
    {
      for(int i=0; i<Dets.size(); i++)
        Dets[i]->acceptMove(P,iat);
    }

    inline ValueType ratio(ParticleSet& P, int iat)
    {
      BFTrans->evaluate(P);

      RealType ratio=1.0;
      for(int i=0; i<Dets.size(); ++i)
        ratio*=Dets[i]->ratio(P,0);
      return ratio; 

    }

    inline ValueType alternateRatio(ParticleSet& P)
    {
      APP_ABORT("Need to implement SlaterDetWithBackflow::alternateRatio() \n");
      return 0.0;
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
      return Dets[DetID[iat]]->update(P,dG,dL,iat);
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
