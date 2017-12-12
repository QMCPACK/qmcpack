//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantWithBackflowBase.h
 * @brief Declaration of DiracDeterminantWithBackflow with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBACKFLOW_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBACKFLOW_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{

/** class to handle determinants with backflow
 */
class DiracDeterminantWithBackflow: public DiracDeterminantBase
{
public:

  typedef SPOSetBase::IndexVector_t IndexVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
  typedef SPOSetBase::HessMatrix_t  HessMatrix_t;
  typedef SPOSetBase::HessVector_t  HessVector_t;
  typedef SPOSetBase::HessType      HessType;
  typedef SPOSetBase::GGGType       GGGType;
  typedef SPOSetBase::GGGVector_t   GGGVector_t;
  typedef SPOSetBase::GGGMatrix_t   GGGMatrix_t;
  typedef SPOSetBase::HessArray_t HessArray_t;
  //typedef Array<GradType,3>       GradArray_t;
  //typedef Array<PosType,3>        PosArray_t;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantWithBackflow(ParticleSet &ptcl, SPOSetBasePtr const &spos, BackflowTransformation * BF, int first=0);

  ///default destructor
  ~DiracDeterminantWithBackflow();

  /**copy constructor
   * @param s existing DiracDeterminantWithBackflow
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantWithBackflow(const DiracDeterminantWithBackflow& s);

  DiracDeterminantWithBackflow& operator=(const DiracDeterminantWithBackflow& s);

  ///** return a clone of Phi
  // */
  //SPOSetBasePtr clonePhi() const;

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    BFTrans = bf;
  }

  // in general, assume that P is the quasiparticle set
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi,
                           ParticleSet::ParticleGradient_t* G0,
                           ParticleSet::ParticleLaplacian_t* L0,
                           int k);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           int offset,
                           Matrix<RealType>& dlogpsi,
                           Array<GradType,OHMMS_DIM>& dG,
                           Matrix<RealType>& dL);

  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  void registerData(ParticleSet& P, WFBufferType& buf);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  ValueType ratio(ParticleSet& P, int iat);

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  ValueType alternateRatio(ParticleSet& P)
  {
    return 1.0;
  }

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

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat);

  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L) ;

  ValueType
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  DiracDeterminantWithBackflow* makeCopy(SPOSetBase* spo) const;

  inline void setLogEpsilon(ValueType x) { }
  inline ValueType rcdot(TinyVector<RealType,OHMMS_DIM>& lhs, TinyVector<ValueType,OHMMS_DIM>& rhs)
  {
    ValueType ret(0);
    for (int i(0); i<OHMMS_DIM; i++)
      ret+=lhs[i]*rhs[i];
    return ret;
  };
#ifdef QMC_COMPLEX
  inline ValueType rcdot(TinyVector<ValueType,OHMMS_DIM>& lhs, TinyVector<RealType,OHMMS_DIM>& rhs)
  {
    ValueType ret(0);
    for (int i(0); i<OHMMS_DIM; i++)
      ret+=lhs[i]*rhs[i];
    return ret;
  };
#endif
  int NumParticles;
  GradMatrix_t dFa;
  HessMatrix_t grad_grad_psiM;
  HessVector_t grad_gradV;
  HessMatrix_t grad_grad_psiM_temp;
  GGGMatrix_t  grad_grad_grad_psiM;
  BackflowTransformation *BFTrans;
  ParticleSet::ParticleGradient_t Gtemp;
  ValueType La1,La2,La3;
  HessMatrix_t Ajk_sum,Qmat;
  GradMatrix_t Fmat;
  GradVector_t Fmatdiag;
  GradVector_t Fmatdiag_temp;

  ValueMatrix_t psiMinv_temp;
  ValueType *FirstAddressOfGGG;
  ValueType *LastAddressOfGGG;
  ValueType *FirstAddressOfFm;
  ValueType *LastAddressOfFm;

  void testDerivFjj(ParticleSet& P, int pa);
  void testGGG(ParticleSet& P);
  void testGG(ParticleSet& P);
  void testDerivLi(ParticleSet& P, int pa);
  void testL(ParticleSet& P);
  void dummyEvalLi(ValueType& L1, ValueType& L2, ValueType& L3);


  void evaluate_SPO(ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);
  void evaluate_SPO(ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

};



}
#endif
