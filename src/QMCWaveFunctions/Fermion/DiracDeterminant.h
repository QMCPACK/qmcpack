//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
///  Refactored from DiracDeterminantBase.h by:
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminant.h
 * @brief Declaration of DiracDeterminant template class with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_H
#define QMCPLUSPLUS_DIRACDETERMINANT_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"

namespace qmcplusplus
{

template<Batching batching = Batching::SINGLE>
class DiracDeterminant;

template<>
class DiracDeterminant<Batching::SINGLE> : public WaveFunctionComponent
{
protected:
  ParticleSet* targetPtcl;
  static constexpr Batching B = Batching::SINGLE;
public:
  bool Optimizable;
  void registerTimers();
  NewTimer UpdateTimer, RatioTimer, InverseTimer, BufferTimer, SPOVTimer, SPOVGLTimer;
  // Optimizable parameters
  opt_variables_type myVars;

  using SSTA = SPOSetTypeAliases;
  using QMCT = QMCTraits;
  using WFCA = WaveFunctionComponentTypeAliases;
  using Walker_t =  WFCA::Walker_t;

  typedef SSTA::IndexVector_t IndexVector_t;
  typedef SSTA::ValueVector_t ValueVector_t;
  typedef SSTA::ValueMatrix_t ValueMatrix_t;
  typedef SSTA::GradVector_t GradVector_t;
  typedef SSTA::GradMatrix_t GradMatrix_t;
  typedef SSTA::HessMatrix_t HessMatrix_t;
  typedef SSTA::HessVector_t HessVector_t;
  typedef SSTA::HessType HessType;

  SPOSet<Batching::SINGLE>* Phi;
#ifdef MIXED_PRECISION
  typedef ParticleSet::SingleParticleValue_t mValueType;
  typedef OrbitalSetTraits<mValueType>::ValueMatrix_t ValueMatrix_hp_t;
#else
  typedef ValueType mValueType;
#endif
  typedef TinyVector<mValueType, DIM> mGradType;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminant(SPOSet<B>* const &spos, int first=0);
  DiracDeterminant(int first = 0);

  ///default destructor
  virtual ~DiracDeterminant();

  /**copy constructor
   * @param s existing DiracDeterminantBase
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminant(const DiracDeterminant& s);

  DiracDeterminant& operator=(const DiracDeterminant& s);

  ///** return a clone of Phi
  // */
  //SPOSetPtr clonePhi() const;
  virtual SPOSet<B>* getPhi() { return dynamic_cast<SPOSet<B>*>(Phi); }
  
  inline IndexType rows() const { return NumPtcls; }

  inline IndexType cols() const { return NumOrbitals; }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  virtual void set(int first, int nel);

  ///set BF pointers
  virtual void setBF(BackflowTransformation* BFTrans) {}


  ///invert psiM or its copies
  virtual void invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat);

  inline void reportStatus(std::ostream& os) {}

  ///reset the size: with the number of particles and number of orbtials
  virtual void resize(int nel, int morb);

  virtual void registerData(ParticleSet& P, WFBufferType& buf);

  virtual void copyFromBuffer(ParticleSet& P, WFBufferType& buf);


  /** move was accepted, update the real container
   */
  virtual void acceptMove(ParticleSet& P, int iat);

  /** move was rejected. copy the real container to the temporary to move on
   */
  virtual void restore(int iat);

  ///evaluate log of determinant for a particle set: should not be called
  virtual RealType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  virtual WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  virtual DiracDeterminant* makeCopy(SPOSet<B>* spo) const;
  //       virtual DiracDeterminantBase* makeCopy(ParticleSet& tqp, SPOSet* spo) const {return makeCopy(spo); };


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

  ValueType
  ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void updateAfterSweep(ParticleSet& P,
			ParticleSet::ParticleGradient_t& G,
			ParticleSet::ParticleLaplacian_t& L);

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

  virtual void resetTargetParticleSet(ParticleSet& P)
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  /** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
  ValueType ratio(ParticleSet& P, int iat);
  
public:
  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  GradType evalGrad(ParticleSet& P, int iat);

  void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  GradType evalGradSource(ParticleSet& P, ParticleSet& source,
			  int iat);

  GradType evalGradSourcep(ParticleSet& P, ParticleSet& source,int iat,
			   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
			   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  GradType evalGradSource(ParticleSet& P, ParticleSet& source,int iat,
			  TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
			  TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);
    
  void recompute(ParticleSet& P);
    
  void NLratios_CPU (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                     std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);

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
  //SPOSetPtr Phi;

  /////Current determinant value
  //ValueType CurrentDet;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM, psiM_temp;

  /// memory pool for temporal data
  aligned_vector<ValueType> memoryPool;

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
  Vector<ParticleSet::SingleParticleValue_t> WorkSpace_hp;
  DiracMatrix<mValueType> detEng_hp;
#endif
  DiracMatrix<ValueType> detEng;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueType curRatio, cumRatio;
  ParticleSet::SingleParticleValue_t* FirstAddressOfG;
  ParticleSet::SingleParticleValue_t* LastAddressOfG;
  ValueType* FirstAddressOfdV;
  ValueType* LastAddressOfdV;
};


} // namespace qmcplusplus
#endif
