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
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file TrialWaveFunction.h
 *@brief Declaration of a TrialWaveFunction
 */
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_H

#include "Message/MPIObjectBase.h"
#include "Particle/VirtualParticleSet.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/FermionBase.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Utilities/NewTimer.h"
/**@defgroup MBWfs Many-body wave function group
 * @brief Classes to handle many-body trial wave functions
 */

namespace qmcplusplus
{

struct CoefficientHolder
{
  typedef OrbitalBase::RealType           RealType;

  std::vector<opt_variables_type> oldVars;
  std::vector<RealType> energies;
  std::vector<RealType> variances;

  void addParams(opt_variables_type& var, RealType e, RealType v)
  {
    oldVars.push_back(var);
    energies.push_back(e);
    variances.push_back(v);
  }

  //this returns the "best" parameters stored in the history. we=energy weight and wv=weightvariance
  opt_variables_type getBestCoefficients(RealType we, RealType wv, bool print=0)
  {
    int best(0);
    RealType bestVal(0);
    bestVal=energies[0]*we + variances[0]*wv;
    if (print)
      app_log()<<0<<" "<<energies[0]<<" "<<variances[0]<<" "<<bestVal<< std::endl;
    for(int ix=1; ix<energies.size(); ix++)
    {
      if (print)
        app_log()<<ix<<" "<<energies[ix]<<" "<<variances[ix]<<" "<<energies[ix]*we+variances[ix]*wv<< std::endl;
      if (energies[ix]*we+variances[ix]*wv < bestVal)
      {
        bestVal=energies[ix]*we+variances[ix]*wv;
        best=ix;
      }
    }
    return oldVars[best];
  }
  opt_variables_type getAvgCoefficients(int lastN)
  {
    opt_variables_type return_params(oldVars[0]);
    for (int i=0; i<return_params.size(); i++)
      return_params[i]=0;
    int start(std::max((int)(oldVars.size()-lastN),0));
    for (int x=start; x<oldVars.size(); x++)
      for (int i=0; i<return_params.size(); i++)
        return_params[i]+=oldVars[x][i];
    RealType nrm(1.0/(oldVars.size()-start));
    for (int i=0; i<return_params.size(); i++)
      return_params[i]*=nrm;
    return return_params;
  }
};

/** @ingroup MBWfs
 * @brief Class to represent a many-body trial wave function
 *
 *A many-body trial wave function is represented by
 *\f$\Psi({\bf R}) = \prod_i \psi_i({\bf R})\f$,
 *where each function \f$\psi_i({\bf R})\f$ is an OrbitalBase
 (see OrbitalComponent).
 *A Composite Pattern is used to handle \f$\prod\f$ operations.
 *Each OrbitalBase should provide proper evaluate functions
 *for the value, gradient and laplacian values.
 */
class TrialWaveFunction: public MPIObjectBase
{

public:

  typedef OrbitalBase::RealType           RealType;
  typedef OrbitalBase::ValueType          ValueType;
  typedef OrbitalBase::PosType            PosType;
  typedef OrbitalBase::GradType           GradType;
  typedef OrbitalBase::BufferType         BufferType;
  typedef OrbitalBase::WFBufferType       WFBufferType;
  typedef OrbitalBase::HessType           HessType;
  typedef OrbitalBase::HessVector_t       HessVector_t;
#ifdef QMC_CUDA
  typedef OrbitalBase::CudaValueType   CudaValueType;
  typedef OrbitalBase::CudaGradType    CudaGradType;
  typedef OrbitalBase::CudaRealType    CudaRealType;
  typedef OrbitalBase::RealMatrix_t    RealMatrix_t;
  typedef OrbitalBase::ValueMatrix_t   ValueMatrix_t;
  typedef OrbitalBase::GradMatrix_t    GradMatrix_t;
  typedef ParticleSet::Walker_t        Walker_t;
#endif

  ///differential gradients
  ParticleSet::ParticleGradient_t G;
  ///differential laplacians
  ParticleSet::ParticleLaplacian_t L;

  TrialWaveFunction(Communicate* c);

  ~TrialWaveFunction();

  inline int size() const
  {
    return Z.size();
  }
  inline RealType getPhase() const
  {
    return PhaseValue;
  }

  inline void setPhase(RealType PhaseValue_new)
  {
    PhaseValue = PhaseValue_new;
  }
  void getLogs(std::vector<RealType>& lvals);
  void getPhases(std::vector<RealType>& pvals);

  inline RealType getAlternatePhaseDiff()
  {
    RealType apd=0.0;
    for (int i=0; i<Z.size(); i++)
    {
      apd += Z[i]->getAlternatePhaseDiff();
    }
    return apd;
  }
  inline RealType getAlternatePhaseDiff(int iat)
  {
    RealType apd=0.0;
    for (int i=0; i<Z.size(); i++)
    {
      apd += Z[i]->getAlternatePhaseDiff(iat);
    }
    return apd;
  }
  inline void alternateGrad(ParticleSet::ParticleGradient_t& G)
  {
    for (int i=0; i<Z.size(); i++)
      Z[i]->alternateGrad(G);
  }

  inline RealType getPhaseDiff() const
  {
    return PhaseDiff;
  }
  inline void resetPhaseDiff()
  {
    PhaseDiff=0.0;
    for (int i=0; i<Z.size(); i++)
      Z[i]->resetPhaseDiff();
  }
  inline RealType getLogPsi() const
  {
    return LogValue;
  }
  inline void setLogPsi(RealType LogPsi_new)
  {
    LogValue = LogPsi_new;
  }

  ///Add an OrbitalBase
  //void addOrbital(OrbitalBase* aterm);
  void addOrbital(OrbitalBase* aterm, const std::string& aname, bool fermion=false);

  ///read from xmlNode
  bool put(xmlNodePtr cur);
  ///implement the virtual function
  void reset();
  /** set OrbitalBase::IsOptimizing to true */
  void startOptimization();
  /** set OrbitalBase::IsOptimizing to flase */
  void stopOptimization();
  /** check in an optimizable parameter
   * * @param o a super set of optimizable variables
   *
   * Update myOptIndex if o is found among the "active" paramemters.
   */
  void checkInVariables(opt_variables_type& o);
  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& o);
  ///reset member data
  void resetParameters(const opt_variables_type& active);
  /** print out state of the trial wavefunction
   */
  void reportStatus(std::ostream& os);

  /** recursively change the ParticleSet whose G and L are evaluated */
  void resetTargetParticleSet(ParticleSet& P);

  /////Check if aname-ed Single-Particle-Orbital set exists
  //bool hasSPOSet(const std::string& aname);
  /////add a Single-Particle-Orbital set
  //void addSPOSet(OhmmsElementBase* spo);
  /////return the aname-ed Single-Particle-Orbital set.
  //OhmmsElementBase* getSPOSet(const std::string& aname);

  /** evalaute the values of the wavefunction, gradient and laplacian  for a walkers */
  ValueType evaluate(ParticleSet& P);

  /** evalaute the values of the wavefunction, gradient and laplacian  for a walkers */
  RealType evaluateLogOnly(ParticleSet& P);

  /** evalaute the log of the trial wave function */
  RealType evaluateLog(ParticleSet& P);
  
  /** recompute the value of the orbitals which require critical accuracy */
  void recompute(ParticleSet& P);

  RealType evaluateDeltaLog(ParticleSet& P, bool recompute=false);

  void evaluateDeltaLog(ParticleSet& P,
                        RealType& logpsi_fixed,
                        RealType& logpsi_opt,
                        ParticleSet::ParticleGradient_t& fixedG,
                        ParticleSet::ParticleLaplacian_t& fixedL);

  /** functions to handle particle-by-particle update */
  RealType ratio(ParticleSet& P, int iat);
  /** evaluate ratios for EE */
  RealType ratioVector(ParticleSet& P, int iat, std::vector<RealType>& ratios);
  /** evaluate ratio for RMC */
  RealType alternateRatio(ParticleSet& P);
  ValueType full_ratio(ParticleSet& P, int iat);

  /** compulte multiple ratios to handle non-local moves and other virtual moves
   */
  void evaluateRatios(VirtualParticleSet& P, std::vector<RealType>& ratios);
  /** compute both ratios and deriatives of ratio with respect to the optimizables*/
  void evaluateDerivRatios(VirtualParticleSet& P, const opt_variables_type& optvars,
      std::vector<RealType>& ratios, Matrix<RealType>& dratio);

  void printGL(ParticleSet::ParticleGradient_t& G,
               ParticleSet::ParticleLaplacian_t& L, std::string tag = "GL");

  /** Returns the logarithmic gradient of the trial wave function
   *  with respect to the iat^th atom of the source ParticleSet. */
  GradType evalGradSource(ParticleSet& P, ParticleSet &source, int iat);
  /** Returns the logarithmic gradient of the w.r.t. the iat^th atom
   * of the source ParticleSet of the sum of laplacians w.r.t. the
   * electrons (target ParticleSet) of the trial wave function. */
  GradType evalGradSource
  (ParticleSet& P, ParticleSet &source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  RealType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  RealType alternateRatioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  GradType evalGrad(ParticleSet& P, int iat);
  GradType alternateEvalGrad(ParticleSet& P, int iat);

  void rejectMove(int iat);
  void acceptMove(ParticleSet& P, int iat);

  /** register all the wavefunction components in buffer.
   *  See OrbitalBase::registerData for more detail */
  void registerData(ParticleSet& P, WFBufferType& buf);
  /** update all the wavefunction components in buffer.
   *  See OrbitalBase::updateBuffer for more detail */
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);
  /** copy all the wavefunction components from buffer.
   *  See OrbitalBase::updateBuffer for more detail */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  RealType KECorrection() const;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi,
                           bool project=false);

  void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
                               std::vector<RealType>& dgradlogpsi);

  /** evalaute the values of the wavefunction, gradient and laplacian  for all the walkers */
  //void evaluate(WalkerSetRef& W, OrbitalBase::ValueVectorType& psi);
  /** evaluate the hessian w.r.t. electronic coordinates of particle iat **/
 // void evaluateHessian(ParticleSet & P, int iat, HessType& grad_grad_psi);
  /** evaluate the hessian hessian w.r.t. electronic coordinates of particle iat **/
  void evaluateHessian(ParticleSet & P, HessVector_t& all_grad_grad_psi);
  
  void reverse();

  inline void resizeTempP(ParticleSet& P)
  {
    tempP = new ParticleSet(P);
  }

  TrialWaveFunction* makeClone(ParticleSet& tqp) const;

  std::vector<OrbitalBase*>& getOrbitals()
  {
    return Z;
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  void setTwist(std::vector<RealType> t)
  {
    myTwist=t;
  }
  const std::vector<RealType> twist()
  {
    return myTwist;
  }

  CoefficientHolder coefficientHistory;

  inline void setMassTerm(ParticleSet& P)
  {
    OneOverM=1.0/P.Mass[0];
    //SpeciesSet tspecies(P.getSpeciesSet());
    //int massind=tspecies.addAttribute("mass");
    //RealType mass = tspecies(massind,0);
    //OneOverM = 1.0/mass;
  }

  FermionBase* getFermionWF()
  {
    return FermionWF;
  }

private:

  ///control how ratio is calculated
  bool Ordered;

  ///the size of ParticleSet
  int NumPtcls;

  ///the size of gradient component (QMCTraits::DIM)*the number of particles
  int TotalDim;

  ///index of the active particle
  int WorkingPtcl;

  ///starting index of the buffer
  size_t BufferCursor;

  ///starting index of the scalar buffer
  size_t BufferCursor_scalar;

  ///sign of the trial wave function
  RealType PhaseValue;

  ///diff of the phase of the trial wave function during ratio calls
  RealType PhaseDiff;

  ///log of the trial wave function
  RealType LogValue;

  ///One over mass of target particleset, needed for Local Energy Derivatives
  RealType OneOverM;

  ///a list of OrbitalBases constituting many-body wave functions
  std::vector<OrbitalBase*> Z;

  ///fermionic wavefunction
  FermionBase* FermionWF;

  ///fake particleset
  ParticleSet* tempP;

  TrialWaveFunction();

  std::vector<NewTimer*> myTimers;
  std::vector<RealType> myTwist;

  ///////////////////////////////////////////
  // Vectorized version for GPU evaluation //
  ///////////////////////////////////////////
#ifdef QMC_CUDA
private:
  gpu::device_host_vector<CudaValueType>   GPUratios;
  gpu::device_host_vector<CudaGradType>    GPUgrads;
  gpu::device_host_vector<CudaValueType>   GPUlapls;

public:
  void freeGPUmem();

  void recompute (MCWalkerConfiguration &W, bool firstTime=true);

  void reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool,
                bool onlyOptimizable=false);
  void getGradient (MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);
  void calcGradient (MCWalkerConfiguration &W, int iat,
                     std::vector<GradType> &grad);
  void addGradient (MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);
  void evaluateLog (MCWalkerConfiguration &W,
                    std::vector<RealType> &logPsi);
  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios);
  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG);
  void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG,
              std::vector<ValueType> &newL);
  void calcRatio (MCWalkerConfiguration &W, int iat,
                  std::vector<ValueType> &psi_ratios,
                  std::vector<GradType> &newG,
                  std::vector<ValueType> &newL);
  void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<ValueType> &psi_ratios,
                 std::vector<GradType> &newG,
                 std::vector<ValueType> &newL);
#ifdef QMC_COMPLEX
  void convertRatiosFromComplexToReal (std::vector<ValueType> &psi_ratios,
                                       std::vector<RealType> &psi_ratios_real);
#endif
  void ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
              std::vector<PosType> &rNew,
              std::vector<ValueType> &psi_ratios,
              std::vector<GradType> &newG,
              std::vector<ValueType> &newL);

  void NLratios (MCWalkerConfiguration &W,
                 gpu::device_vector<CUDA_PRECISION*> &Rlist,
                 gpu::device_vector<int*>            &ElecList,
                 gpu::device_vector<int>             &NumCoreElecs,
                 gpu::device_vector<CUDA_PRECISION*> &QuadPosList,
                 gpu::device_vector<CUDA_PRECISION*> &RatioList,
                 int numQuadPoints);

  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);

  void update (std::vector<Walker_t*> &walkers, int iat);
  void update (const std::vector<Walker_t*> &walkers,
               const std::vector<int> &iatList);

  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);


  void evaluateDeltaLog(MCWalkerConfiguration &W,
                        std::vector<RealType>& logpsi_opt);

  void evaluateDeltaLog(MCWalkerConfiguration &W,
                        std::vector<RealType>& logpsi_fixed,
                        std::vector<RealType>& logpsi_opt,
                        GradMatrix_t&  fixedG,
                        ValueMatrix_t& fixedL);

  void evaluateOptimizableLog (MCWalkerConfiguration &W,
                               std::vector<RealType>& logpsi_opt,
                               GradMatrix_t&  optG,
                               ValueMatrix_t& optL);

  void evaluateDerivatives (MCWalkerConfiguration &W,
                            const opt_variables_type& optvars,
                            RealMatrix_t &dlogpsi,
                            RealMatrix_t &dhpsioverpsi);

#endif


};
/**@}*/
}
#endif
