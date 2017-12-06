//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantBaseBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"

namespace qmcplusplus
{

class DiracDeterminantBase: public OrbitalBase
{
protected:
  ParticleSet *targetPtcl;
public:
  bool Optimizable;
  void registerTimers();
  NewTimer UpdateTimer, RatioTimer, InverseTimer, BufferTimer, SPOVTimer, SPOVGLTimer;
  // Optimizable parameters
  opt_variables_type myVars;


  typedef SPOSetBase::IndexVector_t IndexVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
  typedef SPOSetBase::HessMatrix_t  HessMatrix_t;
  typedef SPOSetBase::HessVector_t  HessVector_t;
  typedef SPOSetBase::HessType      HessType;

#ifdef MIXED_PRECISION
  typedef ParticleSet::ParticleValue_t mValueType;
  typedef OrbitalSetTraits<mValueType>::ValueMatrix_t ValueMatrix_hp_t;
#else
  typedef ValueType mValueType;
#endif
  typedef TinyVector<mValueType,DIM> mGradType;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantBase(SPOSetBasePtr const &spos, int first=0);

  ///default destructor
  virtual ~DiracDeterminantBase();

  /**copy constructor
   * @param s existing DiracDeterminantBase
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantBase(const DiracDeterminantBase& s);

  DiracDeterminantBase& operator=(const DiracDeterminantBase& s);

  ///** return a clone of Phi
  // */
  //SPOSetBasePtr clonePhi() const;

  SPOSetBasePtr getPhi()
  {
    return Phi;
  };

  inline IndexType rows() const
  {
    return NumPtcls;
  }

  inline IndexType cols() const
  {
    return NumOrbitals;
  }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  virtual void set(int first, int nel);
  virtual RealType getAlternatePhaseDiff()
  {
    return 0.0;
  }
  virtual RealType getAlternatePhaseDiff(int iat)
  {
    return 0.0;
  }

  ///set BF pointers
  virtual
  void setBF(BackflowTransformation* BFTrans) {}

  ///optimizations  are disabled
  virtual inline void checkInVariables(opt_variables_type& active)
  {
    Phi->checkInVariables(active);
    Phi->checkInVariables(myVars);
  }

  virtual inline void checkOutVariables(const opt_variables_type& active)
  {
    Phi->checkOutVariables(active);
    myVars.clear();
    myVars.insertFrom(Phi->myVars);
    myVars.getIndex(active);
  }

  virtual void resetParameters(const opt_variables_type& active)
  {
    Phi->resetParameters(active);
    for(int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if(ii>=0)
        myVars[i]= active[ii];
    }
  }

  ///invert psiM or its copies
  void invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat);

  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   std::vector<RealType>& dlogpsi,
                                   std::vector<RealType>& dhpsioverpsi);

  // used by DiracDeterminantWithBackflow
  virtual void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   int offset,
                                   Matrix<RealType>& dlogpsi,
                                   Array<GradType,3>& dG,
                                   Matrix<RealType>& dL) {}

  //virtual void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
  //                                     std::vector<RealType>& dgradlogpsi);

  inline void reportStatus(std::ostream& os)
  {
  }
  virtual void resetTargetParticleSet(ParticleSet& P)
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  ///reset the size: with the number of particles and number of orbtials
  virtual void resize(int nel, int morb);

  virtual void registerData(ParticleSet& P, WFBufferType& buf);

  virtual void updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L);

  virtual RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  virtual ValueType ratio(ParticleSet& P, int iat);

  /** compute multiple ratios for a particle move
   */
  virtual void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  virtual ValueType alternateRatio(ParticleSet& P)
  {
    return 1.0;
  }

  virtual ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  virtual GradType evalGrad(ParticleSet& P, int iat);
  virtual GradType evalGradSource(ParticleSet &P, ParticleSet &source,
                                  int iat);

  virtual GradType evalGradSource
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  virtual GradType evalGradSourcep
  (ParticleSet& P, ParticleSet& source, int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  /** move was accepted, update the real container
   */
  virtual void acceptMove(ParticleSet& P, int iat);

  /** move was rejected. copy the real container to the temporary to move on
   */
  virtual void restore(int iat);

  ///evaluate log of determinant for a particle set: should not be called
  virtual RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L) ;

  virtual void recompute(ParticleSet& P);

  virtual ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L);
           
  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);

  virtual OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  virtual DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;
//       virtual DiracDeterminantBase* makeCopy(ParticleSet& tqp, SPOSetBase* spo) const {return makeCopy(spo); };

  virtual void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);
  ///total number of particles
  int NP;
  ///number of single-particle orbitals which belong to this Dirac determinant
  int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///index of the particle (or row)
  int WorkingIndex;
  ///a set of single-particle orbitals used to fill in the  values of the matrix
  SPOSetBasePtr Phi;

  /////Current determinant value
  //ValueType CurrentDet;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM, psiM_temp;

  /// temporary container for testing
  ValueMatrix_t psiMinv;

  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix_t  dpsiM, dpsiM_temp;

  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  ValueMatrix_t d2psiM, d2psiM_temp;

  /// Used for force computations
  GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
  HessMatrix_t  grad_grad_source_psiM;
  
  GradMatrix_t phi_alpha_Minv, grad_phi_Minv;
  ValueMatrix_t lapl_phi_Minv;
  HessMatrix_t grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector_t psiV;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;
  ValueVector_t workV1, workV2;
  GradVector_t workG;

#ifdef MIXED_PRECISION
  /// temporal matrix and workspace in higher precision for the accurate inversion.
  ValueMatrix_hp_t psiM_hp;
  Vector<ParticleSet::ParticleValue_t> WorkSpace_hp;
  DiracMatrix<mValueType> detEng_hp;
#endif
  DiracMatrix<ValueType> detEng;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueType curRatio,cumRatio;
  ParticleSet::ParticleValue_t *FirstAddressOfG;
  ParticleSet::ParticleValue_t *LastAddressOfG;
  ValueType *FirstAddressOfdV;
  ValueType *LastAddressOfdV;
  //    double ComputeExtraTerms(int ptcl_gradient, int elDim,int ionDim);
  ParticleSet::ParticleGradient_t myG, myG_temp;
  ParticleSet::ParticleLaplacian_t myL, myL_temp;
//
  virtual inline void setLogEpsilon(ValueType x) { }

#ifdef QMC_CUDA
  /////////////////////////////////////////////////////
  // Functions for vectorized evaluation and updates //
  /////////////////////////////////////////////////////
  virtual void recompute(MCWalkerConfiguration &W, bool firstTime)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::recompute.\n";
    abort();
  }

  virtual void
  reserve (PointerPool<gpu::device_vector<CudaValueType> > &pool)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::reserve.\n";
    abort();
  }

  virtual void
  addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::addLog.\n";
    abort();
  }

  virtual void
  ratio (MCWalkerConfiguration &W, int iat,
         std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
         std::vector<ValueType> &lapl)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::ratio.\n";
    abort();
  }
  virtual void
  addRatio (MCWalkerConfiguration &W, int iat,
            std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
            std::vector<ValueType> &lapl)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::addRatio.\n";
    abort();
  }
  virtual void
  calcRatio (MCWalkerConfiguration &W, int iat,
             std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
             std::vector<ValueType> &lapl)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::calcRatio.\n";
    abort();
  }



  virtual void
  ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
         std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
         std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::ratio.\n";
    abort();
  }

  virtual void
  addGradient(MCWalkerConfiguration &W, int iat,
              std::vector<GradType> &grad)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::addGradient.\n";
    abort();
  }
  virtual void
  calcGradient(MCWalkerConfiguration &W, int iat,
               std::vector<GradType> &grad)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::calcGradient.\n";
    abort();
  }


  virtual void update (std::vector<Walker_t*> &walkers, int iat)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::update.\n";
    abort();
  }

  virtual void
  update (const std::vector<Walker_t*> &walkers,
          const std::vector<int> &iatList)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::update.\n";
    abort();
  }

  void
  gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
            ValueMatrix_t &lapl)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::gradLapl.\n";
    abort();
  }

  virtual void
  NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
            std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
  {
    std::cerr << "Need specialization of DiracDetermiantBase::NLratios.\n";
    abort();
  }
#endif
};



}
#endif
