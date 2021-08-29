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
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{
/** class to handle determinants with backflow
 */
class DiracDeterminantWithBackflow : public DiracDeterminantBase
{
public:
  typedef SPOSet::IndexVector_t IndexVector_t;
  typedef SPOSet::ValueVector_t ValueVector_t;
  typedef SPOSet::ValueMatrix_t ValueMatrix_t;
  typedef SPOSet::GradVector_t GradVector_t;
  typedef SPOSet::GradMatrix_t GradMatrix_t;
  typedef SPOSet::HessMatrix_t HessMatrix_t;
  typedef SPOSet::HessVector_t HessVector_t;
  typedef SPOSet::HessType HessType;
  typedef SPOSet::GGGType GGGType;
  typedef SPOSet::GGGVector_t GGGVector_t;
  typedef SPOSet::GGGMatrix_t GGGMatrix_t;
  typedef SPOSet::HessArray_t HessArray_t;
  //typedef Array<GradType,3>       GradArray_t;
  //typedef Array<PosType,3>        PosArray_t;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantWithBackflow(std::shared_ptr<SPOSet>&& spos,
                               std::shared_ptr<BackflowTransformation> BF,
                               int first, int last);

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

  LogValueType evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) override;

  /** cloning function
   * @param tqp target particleset
   * @param spo spo set
   *
   * This interface is exposed only to SlaterDet and its derived classes
   * can overwrite to clone itself correctly.
   */
  DiracDeterminantWithBackflow* makeCopy(std::shared_ptr<SPOSet>&& spo, std::shared_ptr<BackflowTransformation> BF) const;
  DiracDeterminantWithBackflow* makeCopy(std::shared_ptr<SPOSet>&& spo) const override
  {
    throw std::runtime_error("makeCopy spo should not be called.");
    return nullptr;
  }

  void testDerivFjj(ParticleSet& P, int pa);
  void testGGG(ParticleSet& P);
  void testGG(ParticleSet& P);
  void testDerivLi(ParticleSet& P, int pa);
  void testL(ParticleSet& P);

  std::shared_ptr<BackflowTransformation> BFTrans;

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
  GradMatrix_t dFa;
  HessMatrix_t grad_grad_psiM;
  HessVector_t grad_gradV;
  HessMatrix_t grad_grad_psiM_temp;
  GGGMatrix_t grad_grad_grad_psiM;
  ParticleSet::ParticleGradient_t Gtemp;
  ValueType La1, La2, La3;
  HessMatrix_t Ajk_sum, Qmat;
  GradMatrix_t Fmat;
  GradVector_t Fmatdiag;
  GradVector_t Fmatdiag_temp;

  /////Current determinant value
  //ValueType CurrentDet;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM, psiM_temp;

  /// temporary container for testing
  ValueMatrix_t psiMinv;

  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix_t dpsiM, dpsiM_temp;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector_t psiV;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;

  PsiValueType curRatio;
  ParticleSet::SingleParticleValue_t* FirstAddressOfG;
  ParticleSet::SingleParticleValue_t* LastAddressOfG;
  ValueType* FirstAddressOfdV;
  ValueType* LastAddressOfdV;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueMatrix_t psiMinv_temp;
  ValueType* FirstAddressOfGGG;
  ValueType* LastAddressOfGGG;
  ValueType* FirstAddressOfFm;
  ValueType* LastAddressOfFm;

  ParticleSet::ParticleGradient_t myG, myG_temp;
  ParticleSet::ParticleLaplacian_t myL, myL_temp;

  void dummyEvalLi(ValueType& L1, ValueType& L2, ValueType& L3);

  void evaluate_SPO(ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);
  void evaluate_SPO(ValueMatrix_t& logdet,
                    GradMatrix_t& dlogdet,
                    HessMatrix_t& grad_grad_logdet,
                    GGGMatrix_t& grad_grad_grad_logdet);
};


} // namespace qmcplusplus
#endif
