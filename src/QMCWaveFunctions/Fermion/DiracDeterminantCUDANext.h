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
    
    
/**@file DiracDeterminantCUDANext.h
 * @brief Declaration of DiracDeterminantCUDANext with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANT_CUDANEXT_H
#define QMCPLUSPLUS_DIRACDETERMINANT_CUDANEXT_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "QMCWaveFunctions/Fermion/DelayedUpdateCUDA4.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"

namespace qmcplusplus
{

class DiracDeterminantCUDANext: public DiracDeterminant
{
protected:
  ParticleSet *targetPtcl;
  int ndelay;
public:
  bool Optimizable;
  void registerTimers();
  NewTimer UpdateTimer, RatioTimer, InverseTimer, BufferTimer, SPOVTimer, SPOVGLTimer;
  // Optimizable parameters
  opt_variables_type myVars;


  typedef SPOSet::IndexVector_t IndexVector_t;
  typedef SPOSet::ValueVector_t ValueVector_t;
  typedef SPOSet::ValueMatrix_t ValueMatrix_t;
  typedef SPOSet::GradVector_t  GradVector_t;
  typedef SPOSet::GradMatrix_t  GradMatrix_t;
  typedef SPOSet::HessMatrix_t  HessMatrix_t;
  typedef SPOSet::HessVector_t  HessVector_t;
  typedef SPOSet::HessType      HessType;

  typedef ParticleSet::SingleParticleValue_t mValueType;
  typedef OrbitalSetTraits<mValueType>::ValueMatrix_t ValueMatrix_hp_t;
  typedef TinyVector<mValueType,DIM> mGradType;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantCUDANext(SPOSetPtr const &spos, int first=0);

  ///default destructor
   ~DiracDeterminantCUDANext();

  /**copy constructor
   * @param s existing DiracDeterminantCUDANext
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantCUDANext(const DiracDeterminantCUDANext& s);

  DiracDeterminantCUDANext& operator=(const DiracDeterminantCUDANext& s);

  ///** return a clone of Phi
  // */
  //SPOSetPtr clonePhi() const;

  SPOSetPtr getPhi()
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
   void set(int first, int nel, int delay=1);

  ///set BF pointers
  
  void setBF(BackflowTransformation* BFTrans) {}

  ///optimizations  are disabled
   inline void checkInVariables(opt_variables_type& active)
  {
    Phi->checkInVariables(active);
    Phi->checkInVariables(myVars);
  }

   inline void checkOutVariables(const opt_variables_type& active)
  {
    Phi->checkOutVariables(active);
    myVars.clear();
    myVars.insertFrom(Phi->myVars);
    myVars.getIndex(active);
  }

   void resetParameters(const opt_variables_type& active)
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

   void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   std::vector<RealType>& dlogpsi,
                                   std::vector<RealType>& dhpsioverpsi);

  // used by DiracDeterminantCUDANextWithBackflow
   void evaluateDerivatives(ParticleSet& P,
                                   const opt_variables_type& active,
                                   int offset,
                                   Matrix<RealType>& dlogpsi,
                                   Array<GradType,3>& dG,
                                   Matrix<RealType>& dL) {}

  // void evaluateGradDerivatives(const ParticleSet::ParticleGradient_t& G_in,
  //                                     std::vector<RealType>& dgradlogpsi);

  inline void reportStatus(std::ostream& os)
  {
  }
   void resetTargetParticleSet(ParticleSet& P)
  {
    Phi->resetTargetParticleSet(P);
    targetPtcl = &P;
  }

  ///reset the size: with the number of particles and number of orbtials
   void resize(int nel, int morb);

   void registerData(ParticleSet& P, WFBufferType& buf);

   void updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L);

   RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

   void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
   ValueType ratio(ParticleSet& P, int iat);

  /** compute multiple ratios for a particle move
   */
   void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios);

   ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
   GradType evalGrad(ParticleSet& P, int iat);
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
   void completeUpdates(ParticleSet& P);

  /** move was rejected. copy the real container to the temporary to move on
   */
   void restore(int iat);

  ///evaluate log of determinant for a particle set: should not be called
   RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L) ;

   void recompute(ParticleSet& P);

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);

   WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
   DiracDeterminantCUDANext* makeCopy(SPOSet* spo) const;
//        DiracDeterminantCUDANext* makeCopy(ParticleSet& tqp, SPOSet* spo) const {return makeCopy(spo); };

   void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);
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
  SPOSetPtr Phi;

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

  /// temporal matrix in higher precision for the accurate inversion.
  ValueMatrix_hp_t psiM_hp;
  DiracMatrix<mValueType> detEng;
  DelayedUpdateCUDA<ValueType> updateEng;

  ValueType curRatio,cumRatio;
  ParticleSet::SingleParticleValue_t *FirstAddressOfG;
  ParticleSet::SingleParticleValue_t *LastAddressOfG;
  ValueType *FirstAddressOfdV;
  ValueType *LastAddressOfdV;

#ifdef QMC_CUDA
  using WaveFunctionComponent::recompute;
  using WaveFunctionComponent::reserve;
  using WaveFunctionComponent::addLog;
  using WaveFunctionComponent::ratio;
  using WaveFunctionComponent::addRatio;
  using WaveFunctionComponent::calcRatio;
  using WaveFunctionComponent::addGradient;
  using WaveFunctionComponent::calcGradient;
  using WaveFunctionComponent::update;
  using WaveFunctionComponent::gradLapl;
  using WaveFunctionComponent::NLratios;
#endif

};



}
#endif
