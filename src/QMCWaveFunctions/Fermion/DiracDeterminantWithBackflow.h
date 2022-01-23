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


/**@file
 * @brief Declaration of DiracDeterminantWithBackflow with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBACKFLOW_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBACKFLOW_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "Utilities/TimerManager.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{
class BackflowTransformation;

/** class to handle determinants with backflow
 */
class DiracDeterminantWithBackflow : public DiracDeterminantBase
{
public:
  typedef SPOSet::IndexVector IndexVector;
  typedef SPOSet::ValueVector ValueVector;
  typedef SPOSet::ValueMatrix ValueMatrix;
  typedef SPOSet::GradVector GradVector;
  typedef SPOSet::GradMatrix GradMatrix;
  typedef SPOSet::HessMatrix HessMatrix;
  typedef SPOSet::HessVector HessVector;
  typedef SPOSet::HessType HessType;
  typedef SPOSet::GGGType GGGType;
  typedef SPOSet::GGGVector GGGVector;
  typedef SPOSet::GGGMatrix GGGMatrix;
  typedef SPOSet::HessArray HessArray;
  //typedef Array<GradType,3>       GradArray_t;
  //typedef Array<PosType,3>        PosArray_t;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantWithBackflow(std::shared_ptr<SPOSet>&& spos, BackflowTransformation& BF, int first, int last);

  ///default destructor
  ~DiracDeterminantWithBackflow() override;

  // copy constructor and assign operator disabled
  DiracDeterminantWithBackflow(const DiracDeterminantWithBackflow& s) = delete;
  DiracDeterminantWithBackflow& operator=(const DiracDeterminantWithBackflow& s) = delete;

  // in general, assume that P is the quasiparticle set
  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi) override;

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
                           Array<GradType, OHMMS_DIM>& dG,
                           Matrix<RealType>& dL) override;

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  PsiValueType ratio(ParticleSet& P, int iat) override;

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;
  GradType evalGrad(ParticleSet& P, int iat) override;
  GradType evalGradSource(ParticleSet& P, ParticleSet& source, int iat) override;

  GradType evalGradSource(ParticleSet& P,
                          ParticleSet& source,
                          int iat,
                          TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                          TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad) override;

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient_t& G,
                           ParticleSet::ParticleLaplacian_t& L) override;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  std::unique_ptr<DiracDeterminantWithBackflow> makeCopyWithBF(std::shared_ptr<SPOSet>&& spo, BackflowTransformation& BF) const;
  std::unique_ptr<DiracDeterminantBase> makeCopy(std::shared_ptr<SPOSet>&& spo) const override
  {
    throw std::runtime_error("makeCopy spo should not be called.");
    return nullptr;
  }

  void testDerivFjj(ParticleSet& P, int pa);
  void testGGG(ParticleSet& P);
  void testGG(ParticleSet& P);
  void testDerivLi(ParticleSet& P, int pa);
  void testL(ParticleSet& P);

  BackflowTransformation& BFTrans_;

private:
  ///reset the size: with the number of particles and number of orbtials
  void resize(int nel, int morb);

  inline ValueType rcdot(TinyVector<RealType, OHMMS_DIM>& lhs, TinyVector<ValueType, OHMMS_DIM>& rhs)
  {
    ValueType ret(0);
    for (int i(0); i < OHMMS_DIM; i++)
      ret += lhs[i] * rhs[i];
    return ret;
  };

#ifdef QMC_COMPLEX
  inline ValueType rcdot(TinyVector<ValueType, OHMMS_DIM>& lhs, TinyVector<RealType, OHMMS_DIM>& rhs)
  {
    ValueType ret(0);
    for (int i(0); i < OHMMS_DIM; i++)
      ret += lhs[i] * rhs[i];
    return ret;
  };
#endif

  ///total number of particles. Ye: used to track first time allocation but I still feel it very strange.
  int NP;
  int NumParticles;
  GradMatrix dFa;
  HessMatrix grad_grad_psiM;
  HessVector grad_gradV;
  HessMatrix grad_grad_psiM_temp;
  GGGMatrix grad_grad_grad_psiM;
  ParticleSet::ParticleGradient_t Gtemp;
  ValueType La1, La2, La3;
  HessMatrix Ajk_sum, Qmat;
  GradMatrix Fmat;
  GradVector Fmatdiag;
  GradVector Fmatdiag_temp;

  /////Current determinant value
  //ValueType CurrentDet;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix psiM, psiM_temp;

  /// temporary container for testing
  ValueMatrix psiMinv;

  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix dpsiM, dpsiM_temp;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector psiV;
  GradVector dpsiV;
  ValueVector d2psiV;

  PsiValueType curRatio;
  ParticleSet::SingleParticleValue_t* FirstAddressOfG;
  ParticleSet::SingleParticleValue_t* LastAddressOfG;
  ValueType* FirstAddressOfdV;
  ValueType* LastAddressOfdV;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueMatrix psiMinv_temp;
  ValueType* FirstAddressOfGGG;
  ValueType* LastAddressOfGGG;
  ValueType* FirstAddressOfFm;
  ValueType* LastAddressOfFm;

  ParticleSet::ParticleGradient_t myG, myG_temp;
  ParticleSet::ParticleLaplacian_t myL, myL_temp;

  void dummyEvalLi(ValueType& L1, ValueType& L2, ValueType& L3);

  void evaluate_SPO(ValueMatrix& logdet, GradMatrix& dlogdet, HessMatrix& grad_grad_logdet);
  void evaluate_SPO(ValueMatrix& logdet,
                    GradMatrix& dlogdet,
                    HessMatrix& grad_grad_logdet,
                    GGGMatrix& grad_grad_grad_logdet);
};


} // namespace qmcplusplus
#endif
