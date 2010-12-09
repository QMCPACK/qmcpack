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
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#ifdef QMC_CUDA
  #include "QMCWaveFunctions/Fermion/DiracDeterminantCUDA.h"
#else
  #include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#endif
#include <map>

namespace qmcplusplus
{
// NOTE NOTE NOTE
// template<bool backflow>
//  class SlaterDet: public OrbitalBase {}
//     then change SlaterDet to SlaterDet<false>
//     and SlaterDeterminantWithBackflow to SlaterDet<true>
//     and remove all virtuals and inline them 

  class SlaterDet: public OrbitalBase
  {
    public:
    typedef DiracDeterminantBase Determinant_t;
    ///container for the DiracDeterminants
    vector<Determinant_t*>  Dets;
    vector<int> M;
    vector<int> DetID;
    map<string,SPOSetBasePtr> mySPOSet;

    /**  constructor
     * @param targetPtcl target Particleset
     * @param rn release node
     */
    SlaterDet(ParticleSet& targetPtcl);

    ///destructor
    ~SlaterDet();

    ///add a SPOSet 
    void add(SPOSetBasePtr sposet, const string& aname);

    ///add a new DiracDeterminant to the list of determinants
    virtual
    void add(Determinant_t* det, int ispin);

    virtual void checkInVariables(opt_variables_type& active);

    virtual void checkOutVariables(const opt_variables_type& active);

    ///reset all the Dirac determinants, Optimizable is true
    virtual void resetParameters(const opt_variables_type& optVariables);

    void reportStatus(ostream& os);

    void resetTargetParticleSet(ParticleSet& P);

    virtual
    ValueType evaluate(ParticleSet& P
          ,ParticleSet::ParticleGradient_t& G
          ,ParticleSet::ParticleLaplacian_t& L);

    virtual
    RealType evaluateLog(ParticleSet& P
          ,ParticleSet::ParticleGradient_t& G
          ,ParticleSet::ParticleLaplacian_t& L);
    
    virtual    
    RealType evaluateLog(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G, 
      ParticleSet::ParticleLaplacian_t& L,
      PooledData<RealType>& buf,
      bool fillBuffer);
      
    void registerDataForDerivatives(ParticleSet& P, BufferType& buf);

    ///return the total number of Dirac determinants
    inline int size() const
    {
      return Dets.size();
    }

    ///return the column dimension of the i-th Dirac determinant
    inline int size(int i) const
    {
      return Dets[i]->cols();
    }

    virtual
    RealType registerData(ParticleSet& P, PooledData<RealType>& buf);

    virtual
    RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

    virtual
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    virtual
    void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);

    virtual
    void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

    virtual
    RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

    virtual
    inline ValueType ratio(ParticleSet& P, int iat,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL)
    {
      return Dets[DetID[iat]]->ratio(P,iat,dG,dL);
    }

    virtual
    inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      return Dets[DetID[iat]]->ratioGrad(P,iat,grad_iat);
    }

    virtual
    inline ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      return Dets[DetID[iat]]->alternateRatioGrad(P,iat,grad_iat);
    }

    virtual
    GradType evalGrad(ParticleSet& P, int iat)
    {
      return Dets[DetID[iat]]->evalGrad(P,iat);
    }

    virtual
    GradType alternateEvalGrad(ParticleSet& P, int iat)
    {
      return Dets[DetID[iat]]->alternateEvalGrad(P,iat);
    }

    virtual
    GradType evalGradSource(ParticleSet& P, ParticleSet &src, int iat)
    {
      GradType G = GradType();
      for (int iz=0; iz < size(); iz++)
        G += Dets[iz]->evalGradSource(P, src, iat);
      return G;
    }

    virtual
    GradType evalGradSource (ParticleSet& P, ParticleSet& src, int iat,
        TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
        TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad) 
    {
      GradType G = GradType(); 
      for (int iz=0; iz < size(); iz++) 
        G += Dets[iz]->evalGradSource(P, src, iat, grad_grad, lapl_grad); 
      return G; 
    }

    virtual
    inline ValueType logRatio(ParticleSet& P, int iat,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL)
    {
      ValueType r = Dets[DetID[iat]]->ratio(P,iat,dG,dL);
      return evaluateLogAndPhase(r,PhaseValue);
    }

    virtual
    inline void restore(int iat)
    {
      return Dets[DetID[iat]]->restore(iat);
    }

    RealType getAlternatePhaseDiff()
    {
      RealType ap(0.0);
      for (int iz=0; iz < size(); iz++) ap += Dets[iz]->getAlternatePhaseDiff();
      return ap;
    }

    RealType getAlternatePhaseDiff(int iat)
    {
      return Dets[DetID[iat]]->getAlternatePhaseDiff();
    }

    virtual
    inline void acceptMove(ParticleSet& P, int iat)
    {
      Dets[DetID[iat]]->acceptMove(P,iat);
    }

    virtual
    inline ValueType ratio(ParticleSet& P, int iat)
    {
      return Dets[DetID[iat]]->ratio(P,iat);
    }

    virtual
    inline ValueType alternateRatio(ParticleSet& P)
    {
      ValueType v(1.0);
      for(int i=0; i<Dets.size(); ++i) v *= Dets[i]->alternateRatio(P);
      return v;
    }

    virtual
    inline void alternateGrad(ParticleSet::ParticleGradient_t& G)
    {
      for(int i=0; i<Dets.size(); ++i) Dets[i]->alternateGrad(G);
    }

    virtual
    void update(ParticleSet& P,
        ParticleSet::ParticleGradient_t& dG,
        ParticleSet::ParticleLaplacian_t& dL,
        int iat)
    {
      return Dets[DetID[iat]]->update(P,dG,dL,iat);
    }

    virtual
    OrbitalBasePtr makeClone(ParticleSet& tqp) const;

    virtual
    SPOSetBasePtr getPhi(int i=0)
    {
      return Dets[i]->getPhi();
    }

    virtual
    void get_ratios(ParticleSet& P, vector<ValueType>& ratios);

    void evaluateDerivatives(ParticleSet& P,
			     const opt_variables_type& active,
			     vector<RealType>& dlogpsi,
			     vector<RealType>& dhpsioverpsi)
    {
      // First zero out values, since each determinant only adds on
      // its contribution (i.e. +=) , rather than setting the value
      // (i.e. =)
      for (int k=0; k<myVars.size(); ++k) {
	int kk=myVars.where(k);
	if (kk >= 0)
	  dlogpsi[kk] = dhpsioverpsi[kk] = 0.0;
      }

      
      // Now add on contribution from each determinant to the derivatives
      for (int i=0; i<Dets.size(); i++)
	Dets[i]->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
    }

#ifdef QMC_CUDA
    /////////////////////////////////////////////////////
    // Functions for vectorized evaluation and updates //
    /////////////////////////////////////////////////////
    void recompute(MCWalkerConfiguration &W, bool firstTime)
    {
      for (int id=0; id<Dets.size(); id++)
        Dets[id]->recompute(W, firstTime);
    }

    void reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool)
    {
      for (int id=0; id<Dets.size(); id++)
        Dets[id]->reserve(pool);
    }

    void addLog (MCWalkerConfiguration &W, vector<RealType> &logPsi)
    {
      for (int id=0; id<Dets.size(); id++)
        Dets[id]->addLog(W, logPsi);
    }

    void 
      ratio (MCWalkerConfiguration &W, int iat
          , vector<ValueType> &psi_ratios,vector<GradType>  &grad, vector<ValueType> &lapl)
      {
        Dets[DetID[iat]]->ratio(W, iat, psi_ratios, grad, lapl);
      }

    void ratio (vector<Walker_t*> &walkers,    vector<int> &iatList,
        vector<PosType> &rNew, vector<ValueType> &psi_ratios, 
        vector<GradType>  &grad, vector<ValueType> &lapl);

    void addGradient(MCWalkerConfiguration &W, int iat, vector<GradType> &grad)
    {
      Dets[DetID[iat]]->addGradient(W, iat, grad);
    }

    void update (vector<Walker_t*> &walkers, int iat)
    {
      Dets[DetID[iat]]->update(walkers, iat);
    }

    void update (const vector<Walker_t*> &walkers, const vector<int> &iatList);

    void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads, ValueMatrix_t &lapl) 
    {
      for (int id=0; id<Dets.size(); id++) Dets[id]->gradLapl(W, grads, lapl);
    }

    void NLratios (MCWalkerConfiguration &W,  vector<NLjob> &jobList
        , vector<PosType> &quadPoints, vector<ValueType> &psi_ratios)
    {
      for (int id=0; id<Dets.size(); id++) 
        Dets[id]->NLratios(W, jobList, quadPoints, psi_ratios);
    }
#endif  

    private:
    SlaterDet() {}
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
