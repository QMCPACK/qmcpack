//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#define QMCPLUSPLUS_SLATERDETERMINANT_WITHBASE_H
#include "QMCWaveFunctions/FermionBase.h"
#ifdef QMC_CUDA
#include "QMCWaveFunctions/Fermion/DiracDeterminantCUDA.h"
#else
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#endif
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include <map>

namespace qmcplusplus
{
// NOTE NOTE NOTE
// template<bool backflow>
//  class SlaterDet: public OrbitalBase {}
//     then change SlaterDet to SlaterDet<false>
//     and SlaterDeterminantWithBackflow to SlaterDet<true>
//     and remove all virtuals and inline them

class SlaterDet: public OrbitalBase, public FermionBase
{
public:
  typedef DiracDeterminantBase Determinant_t;
  ///container for the DiracDeterminants
  std::vector<Determinant_t*>  Dets;
  ///the last particle of each group
  std::vector<int> Last;
  std::map<std::string,SPOSetBasePtr> mySPOSet;

  /**  constructor
   * @param targetPtcl target Particleset
   * @param rn release node
   */
  SlaterDet(ParticleSet& targetPtcl);

  ///destructor
  ~SlaterDet();

  //get Det ID
  inline int getDetID(const int iat)
  {
    int id=0;
    while(iat>Last[id])
      id++;
    return id;
  }

  ///add a SPOSet
  void add(SPOSetBasePtr sposet, const std::string& aname);

  ///add a new DiracDeterminant to the list of determinants
  virtual
  void add(Determinant_t* det, int ispin);

  ///set BF pointers
  virtual
  void setBF(BackflowTransformation* BFTrans) {}

  virtual void checkInVariables(opt_variables_type& active);

  virtual void checkOutVariables(const opt_variables_type& active);

  ///reset all the Dirac determinants, Optimizable is true
  virtual void resetParameters(const opt_variables_type& optVariables);

  void reportStatus(std::ostream& os);

  virtual void resetTargetParticleSet(ParticleSet& P);

  virtual
  ValueType evaluate(ParticleSet& P
                     ,ParticleSet::ParticleGradient_t& G
                     ,ParticleSet::ParticleLaplacian_t& L);

  virtual
  RealType evaluateLog(ParticleSet& P
                       ,ParticleSet::ParticleGradient_t& G
                       ,ParticleSet::ParticleLaplacian_t& L);

  virtual void recompute(ParticleSet& P);

  virtual
  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);

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

  virtual void registerData(ParticleSet& P, WFBufferType& buf);

  virtual void updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L);

  virtual
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  virtual
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  virtual
  inline void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    return Dets[VP.activeGroup]->evaluateRatios(VP,ratios);
  }

  virtual
  inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    return Dets[getDetID(iat)]->ratioGrad(P,iat,grad_iat);
  }

  virtual
  inline ValueType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    return Dets[getDetID(iat)]->alternateRatioGrad(P,iat,grad_iat);
  }

  virtual
  GradType evalGrad(ParticleSet& P, int iat)
  {
    return Dets[getDetID(iat)]->evalGrad(P,iat);
  }

  virtual
  GradType alternateEvalGrad(ParticleSet& P, int iat)
  {
    return Dets[getDetID(iat)]->alternateEvalGrad(P,iat);
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
  inline void restore(int iat)
  {
    return Dets[getDetID(iat)]->restore(iat);
  }

  RealType getAlternatePhaseDiff()
  {
    RealType ap(0.0);
    for (int iz=0; iz < size(); iz++)
      ap += Dets[iz]->getAlternatePhaseDiff();
    return ap;
  }

  RealType getAlternatePhaseDiff(int iat)
  {
    return Dets[getDetID(iat)]->getAlternatePhaseDiff();
  }

  virtual
  inline void acceptMove(ParticleSet& P, int iat)
  {
    Dets[getDetID(iat)]->acceptMove(P,iat);

    LogValue=0.0;
    PhaseValue=0.0;
    for(int i=0; i<Dets.size(); ++i) {
      LogValue+= Dets[i]->LogValue;
      PhaseValue+= Dets[i]->PhaseValue;
    }
  }

  virtual
  inline ValueType ratio(ParticleSet& P, int iat)
  {
    return Dets[getDetID(iat)]->ratio(P,iat);
  }

  virtual
  inline ValueType alternateRatio(ParticleSet& P)
  {
    ValueType v(1.0);
    for(int i=0; i<Dets.size(); ++i)
      v *= Dets[i]->alternateRatio(P);
    return v;
  }

  virtual
  inline void alternateGrad(ParticleSet::ParticleGradient_t& G)
  {
    for(int i=0; i<Dets.size(); ++i)
      Dets[i]->alternateGrad(G);
  }

  virtual
  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  virtual
  SPOSetBasePtr getPhi(int i=0)
  {
    return Dets[i]->getPhi();
  }

  virtual
  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi)
  {
    // First zero out values, since each determinant only adds on
    // its contribution (i.e. +=) , rather than setting the value
    // (i.e. =)
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk >= 0)
        dlogpsi[kk] = dhpsioverpsi[kk] = 0.0;
    }
    // Now add on contribution from each determinant to the derivatives
    for (int i=0; i<Dets.size(); i++)
      Dets[i]->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
  }

  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                               std::vector<RealType>& dgradlogpsi)
  {
    for (int i=0; i<Dets.size(); i++)
      Dets[i]->evaluateGradDerivatives(G_in, dgradlogpsi);
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

  void reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->reserve(pool);
  }

  void addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->addLog(W, logPsi);
  }

  void
  ratio (MCWalkerConfiguration &W, int iat
         , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->ratio(W, iat, psi_ratios, grad, lapl);
  }

  void
  calcRatio (MCWalkerConfiguration &W, int iat
             , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->calcRatio(W, iat, psi_ratios, grad, lapl);
  }

  void
  addRatio (MCWalkerConfiguration &W, int iat
            , std::vector<ValueType> &psi_ratios,std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    Dets[getDetID(iat)]->addRatio(W, iat, psi_ratios, grad, lapl);
  }


  void ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
              std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
              std::vector<GradType>  &grad, std::vector<ValueType> &lapl);

  void calcGradient(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
  {
    Dets[getDetID(iat)]->calcGradient(W, iat, grad);
  }

  void addGradient(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
  {
    Dets[getDetID(iat)]->addGradient(W, iat, grad);
  }

  void update (std::vector<Walker_t*> &walkers, int iat)
  {
    Dets[getDetID(iat)]->update(walkers, iat);
  }

  void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads, ValueMatrix_t &lapl)
  {
    for (int id=0; id<Dets.size(); id++)
      Dets[id]->gradLapl(W, grads, lapl);
  }

  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList
                 , std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
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
