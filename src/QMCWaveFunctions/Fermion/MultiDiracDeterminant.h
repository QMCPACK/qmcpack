//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////


/**@file DiracDeterminant.h
 * @brief Declaration of DiracDeterminant with a S(ingle)P(article)O(rbital)Set
 */
#ifndef QMCPLUSPLUS_MULTIDIRACDETERMINANT_H
#define QMCPLUSPLUS_MULTIDIRACDETERMINANT_H
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantCalculator.h"
#include "Message/Communicate.h"
#include "Numerics/DeterminantOperators.h"
//#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
class MultiDiracDeterminant : public WaveFunctionComponent
{
public:
  bool Optimizable;
  void registerTimers();
  NewTimer &UpdateTimer, &RatioTimer, &MWRatioTimer, &InverseTimer, &buildTableTimer, &readMatTimer, &evalWTimer,
      &evalOrbTimer, &evalOrb1Timer;
  NewTimer &readMatGradTimer, &buildTableGradTimer, &ExtraStuffTimer;
  // Optimizable parameter
  opt_variables_type myVars;

  using IndexVector = SPOSet::IndexVector;
  using ValueVector = SPOSet::ValueVector;
  using ValueMatrix = SPOSet::ValueMatrix;
  using GradVector  = SPOSet::GradVector;
  using GradMatrix  = SPOSet::GradMatrix;
  using HessMatrix  = SPOSet::HessMatrix;
  using HessType    = SPOSet::HessType;

  //lookup table mapping the unique determinants to their element position in C2_node vector
  std::vector<std::vector<int>> lookup_tbl;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  MultiDiracDeterminant(std::unique_ptr<SPOSet>&& spos, bool spinor);

  ///default destructor
  ~MultiDiracDeterminant() override;

  /**copy constructor
   * @param s existing DiracDeterminant
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  MultiDiracDeterminant(const MultiDiracDeterminant& s);

  MultiDiracDeterminant& operator=(const MultiDiracDeterminant& s) = delete;

  /** return a clone of Phi
   */
  std::unique_ptr<SPOSet> clonePhi() const;

  SPOSetPtr getPhi() { return Phi.get(); };

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   *@param ref_det_id id of the reference determinant
   *
   * Note: ciConfigList should have been populated when calling this function
   */
  void set(int first, int nel, int ref_det_id);

  ///optimizations  are disabled
  inline void checkInVariables(opt_variables_type& active) override { Phi->checkInVariables(active); }

  inline void checkOutVariables(const opt_variables_type& active) override { Phi->checkOutVariables(active); }

  /// create optimizable orbital rotation parameters
  void buildOptVariables(std::vector<size_t>& C2node);
  ///helper function to buildOptVariables
  int build_occ_vec(const std::vector<int>& data, const size_t nel, const size_t nmo, std::vector<int>& occ_vec);

  void resetParameters(const opt_variables_type& active) override { Phi->resetParameters(active); }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<ValueType>& dlogpsi,
                           std::vector<ValueType>& dhpsioverpsi,
                           const MultiDiracDeterminant& pseudo_dn,
                           const ValueType& psiCurrent,
                           const std::vector<ValueType>& Coeff,
                           const std::vector<size_t>& C2node_up,
                           const std::vector<size_t>& C2node_dn)
  {
    if (!Optimizable)
      return;

    const ValueVector& detValues_up = getRatiosToRefDet();
    const ValueVector& detValues_dn = pseudo_dn.getRatiosToRefDet();
    const GradMatrix& grads_up      = grads;
    const GradMatrix& grads_dn      = pseudo_dn.grads;
    const ValueMatrix& lapls_up     = lapls;
    const ValueMatrix& lapls_dn     = pseudo_dn.lapls;
    const ValueMatrix& M_up         = psiM;
    const ValueMatrix& M_dn         = pseudo_dn.psiM;
    const ValueMatrix& Minv_up      = psiMinv;
    const ValueMatrix& Minv_dn      = pseudo_dn.psiMinv;
    const GradMatrix& B_grad        = dpsiM;
    const ValueMatrix& B_lapl       = d2psiM;

    const size_t N1  = FirstIndex;
    const size_t N2  = pseudo_dn.FirstIndex;
    const size_t NP1 = NumPtcls;
    const size_t NP2 = pseudo_dn.NumPtcls;

    Phi->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, psiCurrent, Coeff, C2node_up, C2node_dn, detValues_up,
                             detValues_dn, grads_up, grads_dn, lapls_up, lapls_dn, M_up, M_dn, Minv_up, Minv_dn, B_grad,
                             B_lapl, *detData, N1, N2, NP1, NP2, lookup_tbl);
  }

  void evaluateDerivativesWF(ParticleSet& P,
                             const opt_variables_type& optvars,
                             std::vector<ValueType>& dlogpsi,
                             const MultiDiracDeterminant& pseudo_dn,
                             const PsiValueType& psiCurrent,
                             const std::vector<ValueType>& Coeff,
                             const std::vector<size_t>& C2node_up,
                             const std::vector<size_t>& C2node_dn)
  {
    if (!Optimizable)
      return;

    const ValueVector& detValues_up = getRatiosToRefDet();
    const ValueVector& detValues_dn = pseudo_dn.getRatiosToRefDet();
    const ValueMatrix& M_up         = psiM;
    const ValueMatrix& M_dn         = pseudo_dn.psiM;
    const ValueMatrix& Minv_up      = psiMinv;
    const ValueMatrix& Minv_dn      = pseudo_dn.psiMinv;

    Phi->evaluateDerivativesWF(P, optvars, dlogpsi, psiCurrent, Coeff, C2node_up, C2node_dn, detValues_up, detValues_dn,
                               M_up, M_dn, Minv_up, Minv_dn, *detData, lookup_tbl);
  }


  inline void reportStatus(std::ostream& os) override {}

  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  /****************************************************************************
   * These functions should not be called.
   ***************************************************************************/

  PsiValueType ratio(ParticleSet& P, int iat) override
  {
    APP_ABORT("  MultiDiracDeterminant: This should not be called. \n");
    return PsiValueType();
  }

  GradType evalGrad(ParticleSet& P, int iat) override
  {
    APP_ABORT("  MultiDiracDeterminant: This should not be called. \n");
    return GradType();
  }

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override
  {
    APP_ABORT("  MultiDiracDeterminant: This should not be called. \n");
    return PsiValueType();
  }

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override
  {
    APP_ABORT("  MultiDiracDeterminant: This should not be called. \n");
    return 0.0;
  }

  ValueType evaluate(const ParticleSet& P, ParticleSet::ParticleGradient& G, ParticleSet::ParticleLaplacian& L)
  {
    APP_ABORT("  MultiDiracDeterminant: This should not be called. \n");
    return ValueType();
  }


  /****************************************************************************
   * END END END
   ***************************************************************************/

  // create necessary structures used in the evaluation of the determinants
  // this works with confgList, which shouldn't change during a simulation
  void createDetData(const ci_configuration2& ref,
                     std::vector<int>& data,
                     std::vector<std::pair<int, int>>& pairs,
                     std::vector<RealType>& sign);

  template<typename ITER>
  inline ValueType CalculateRatioFromMatrixElements(int n, ValueMatrix& dotProducts, ITER it)
  {
    switch (n)
    {
    case 0:
      return 1.0;
    case 1:
      return dotProducts(*it, *(it + 1));
    case 2: {
      const int i = *it;
      const int j = *(it + 1);
      const int a = *(it + 2);
      const int b = *(it + 3);
      return dotProducts(i, a) * dotProducts(j, b) - dotProducts(i, b) * dotProducts(j, a);
    }
    case 3: {
      const int i1 = *it;
      const int i2 = *(it + 1);
      const int i3 = *(it + 2);
      const int a1 = *(it + 3);
      const int a2 = *(it + 4);
      const int a3 = *(it + 5);
      return DetCalculator.evaluate(dotProducts(i1, a1), dotProducts(i1, a2), dotProducts(i1, a3), dotProducts(i2, a1),
                                    dotProducts(i2, a2), dotProducts(i2, a3), dotProducts(i3, a1), dotProducts(i3, a2),
                                    dotProducts(i3, a3));
    }
    case 4: {
      const int i1 = *it;
      const int i2 = *(it + 1);
      const int i3 = *(it + 2);
      const int i4 = *(it + 3);
      const int a1 = *(it + 4);
      const int a2 = *(it + 5);
      const int a3 = *(it + 6);
      const int a4 = *(it + 7);
      return DetCalculator.evaluate(dotProducts(i1, a1), dotProducts(i1, a2), dotProducts(i1, a3), dotProducts(i1, a4),
                                    dotProducts(i2, a1), dotProducts(i2, a2), dotProducts(i2, a3), dotProducts(i2, a4),
                                    dotProducts(i3, a1), dotProducts(i3, a2), dotProducts(i3, a3), dotProducts(i3, a4),
                                    dotProducts(i4, a1), dotProducts(i4, a2), dotProducts(i4, a3), dotProducts(i4, a4));
    }
    case 5: {
      const int i1 = *it;
      const int i2 = *(it + 1);
      const int i3 = *(it + 2);
      const int i4 = *(it + 3);
      const int i5 = *(it + 4);
      const int a1 = *(it + 5);
      const int a2 = *(it + 6);
      const int a3 = *(it + 7);
      const int a4 = *(it + 8);
      const int a5 = *(it + 9);
      return DetCalculator.evaluate(dotProducts(i1, a1), dotProducts(i1, a2), dotProducts(i1, a3), dotProducts(i1, a4),
                                    dotProducts(i1, a5), dotProducts(i2, a1), dotProducts(i2, a2), dotProducts(i2, a3),
                                    dotProducts(i2, a4), dotProducts(i2, a5), dotProducts(i3, a1), dotProducts(i3, a2),
                                    dotProducts(i3, a3), dotProducts(i3, a4), dotProducts(i3, a5), dotProducts(i4, a1),
                                    dotProducts(i4, a2), dotProducts(i4, a3), dotProducts(i4, a4), dotProducts(i4, a5),
                                    dotProducts(i5, a1), dotProducts(i5, a2), dotProducts(i5, a3), dotProducts(i5, a4),
                                    dotProducts(i5, a5));
    }
    default:
      return DetCalculator.evaluate(dotProducts, it, n);
    }
    return 0.0;
  }

  void BuildDotProductsAndCalculateRatios_impl(int ref,
                                               ValueType det0,
                                               ValueType* restrict ratios,
                                               const ValueMatrix& psiinv,
                                               const ValueMatrix& psi,
                                               ValueMatrix& dotProducts,
                                               const std::vector<int>& data,
                                               const std::vector<std::pair<int, int>>& pairs,
                                               const std::vector<RealType>& sign);

  void BuildDotProductsAndCalculateRatios(int ref,
                                          ValueVector& ratios,
                                          const ValueMatrix& psiinv,
                                          const ValueMatrix& psi,
                                          ValueMatrix& dotProducts,
                                          const std::vector<int>& data,
                                          const std::vector<std::pair<int, int>>& pairs,
                                          const std::vector<RealType>& sign);

  void BuildDotProductsAndCalculateRatios(int ref,
                                          int iat,
                                          GradMatrix& ratios,
                                          ValueMatrix& psiinv,
                                          ValueMatrix& psi,
                                          ValueMatrix& dotProducts,
                                          std::vector<int>& data,
                                          std::vector<std::pair<int, int>>& pairs,
                                          std::vector<RealType>& sign,
                                          int dx);

  void BuildDotProductsAndCalculateRatios(int ref,
                                          int iat,
                                          ValueMatrix& ratios,
                                          ValueMatrix& psiinv,
                                          ValueMatrix& psi,
                                          ValueMatrix& dotProducts,
                                          std::vector<int>& data,
                                          std::vector<std::pair<int, int>>& pairs,
                                          std::vector<RealType>& sign);

  //   Finish this at some point
  inline void InverseUpdateByColumn_GRAD(ValueMatrix& Minv,
                                         GradVector& newcol,
                                         ValueVector& rvec,
                                         ValueVector& rvecinv,
                                         int colchanged,
                                         ValueType c_ratio,
                                         int dx)
  {
    c_ratio = (RealType)1.0 / c_ratio;
    int m   = Minv.rows();
    BLAS::gemv('N', m, m, c_ratio, Minv.data(), m, newcol[0].data() + dx, 3, 0.0, rvec.data(), 1);
    rvec[colchanged] = (RealType)1.0 - c_ratio;
    BLAS::copy(m, Minv.data() + colchanged, m, rvecinv.data(), 1);
    BLAS::ger(m, m, -1.0, rvec.data(), 1, rvecinv.data(), 1, Minv.data(), m);
  }
  /*
      inline void InverseUpdateByColumn(ValueMatrix& Minv
            , GradVector& dM, ValueVector& rvec
            , ValueVector& rvecinv, int colchanged
            , ValueType c_ratio, std::vector<int>::iterator& it)
      {
            ValueType c_ratio=1.0/ratioLapl;
            BLAS::gemv('N', NumPtcls, NumPtcls, c_ratio, dpsiMinv.data(), NumPtcls, tv, 1, T(), workV1, 1);
            workV1[colchanged]=1.0-c_ratio;
            BLAS::copy(m,pinv+colchanged,m,workV2,1);
            BLAS::ger(m,m,-1.0,workV1,1,workV2,1,dpsiMinv.data(),m);
      }
  */

  /** evaluate the value of all the unique determinants with one electron moved. Used by the table method
   *@param P particle set which provides the positions
   *@param iat the index of the moved electron
   *@param refPtcl if given, the id of the reference particle in virtual moves
   */
  void evaluateDetsForPtclMove(const ParticleSet& P, int iat, int refPtcl = -1);
  /// multi walker version of evaluateDetsForPtclMove
  void static mw_evaluateDetsForPtclMove(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                                         const RefVectorWithLeader<ParticleSet>& P_list,
                                         int iat);

  /// evaluate the value and gradients of all the unique determinants with one electron moved. Used by the table method
  void evaluateDetsAndGradsForPtclMove(const ParticleSet& P, int iat);
  /// multi walker version of mw_evaluateDetsAndGradsForPtclMove
  void static mw_evaluateDetsAndGradsForPtclMove(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                                                 const RefVectorWithLeader<ParticleSet>& P_list,
                                                 int iat);
  /// evaluate the value and gradients of all the unique determinants with one electron moved. Used by the table method. Includes Spin Gradient data
  void evaluateDetsAndGradsForPtclMoveWithSpin(const ParticleSet& P, int iat);


  /// evaluate the gradients of all the unique determinants with one electron moved. Used by the table method
  void evaluateGrads(ParticleSet& P, int iat);
  /// multi walker version of mw_evaluateGrads
  void static mw_evaluateGrads(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                               const RefVectorWithLeader<ParticleSet>& P_list,
                               int iat);
  /// evaluate the gradients of all the unique determinants with one electron moved. Used by the table method. Includes Spin Gradient data
  void evaluateGradsWithSpin(ParticleSet& P, int iat);

  // full evaluation of all the structures from scratch, used in evaluateLog for example
  void evaluateForWalkerMove(const ParticleSet& P, bool fromScratch = true);
  // full evaluation of all the structures from scratch, used in evaluateLog for example. Includes spin gradients for spin moves
  void evaluateForWalkerMoveWithSpin(const ParticleSet& P, bool fromScratch = true);

  // accessors
  inline int getNumDets() const { return ciConfigList->size(); }
  inline int getNumPtcls() const { return NumPtcls; }
  inline int getFirstIndex() const { return FirstIndex; }
  inline std::vector<ci_configuration2>& getCIConfigList() { return *ciConfigList; }

  const ValueVector& getRatiosToRefDet() const { return ratios_to_ref_; }
  const ValueVector& getNewRatiosToRefDet() const { return new_ratios_to_ref_; }
  const GradMatrix& getGrads() const { return grads; }
  const GradMatrix& getNewGrads() const { return new_grads; }
  const ValueMatrix& getLapls() const { return lapls; }
  const ValueMatrix& getNewLapls() const { return new_lapls; }
  const ValueMatrix& getSpinGrads() const { return spingrads; }
  const ValueMatrix& getNewSpinGrads() const { return new_spingrads; }

  PsiValueType getRefDetRatio() const { return static_cast<PsiValueType>(curRatio); }
  LogValueType getLogValueRefDet() const { return log_value_ref_det_; }

private:
  ///reset the size: with the number of particles
  void resize(int nel);

  ///a set of single-particle orbitals used to fill in the  values of the matrix
  const std::unique_ptr<SPOSet> Phi;
  ///number of single-particle orbitals which belong to this Dirac determinant
  const int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///use shared_ptr
  std::shared_ptr<std::vector<ci_configuration2>> ciConfigList;
  // the reference determinant never changes, so there is no need to store it.
  // if its value is zero, then use a data from backup, but always use this one
  // by default
  int ReferenceDeterminant;
  // flag to determine if spin arrays need to be resized and used. Set by ParticleSet::is_spinor_ in SlaterDetBuilder
  const bool is_spinor_;

  /// psiM(i,j) \f$= \psi_j({\bf r}_i)\f$
  /// TpsiM(i,j) \f$= psiM(j,i) \f$
  ValueMatrix psiM, TpsiM;
  /// inverse Dirac determinant matrix of the reference det
  ValueMatrix psiMinv, psiMinv_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix dpsiM;
  // temporaty storage
  ValueMatrix dpsiMinv;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  ValueMatrix d2psiM;
  /* dspin_psiM(i,j) \f$= \partial_{s_i} \psi_j({\bf r}_i,s_i)\f$ where \f$s_i\f$s is the spin variable
   * This is only resized if a spinor calculation is used
   */
  ValueMatrix dspin_psiM;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector psiV, psiV_temp;
  GradVector dpsiV;
  ValueVector d2psiV;
  ValueVector workV1, workV2;
  //spin  derivative of single-particle orbitals. Only resized if a spinor calculation
  ValueVector dspin_psiV;

  ValueMatrix dotProducts;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueType* FirstAddressOfGrads;
  ValueType* LastAddressOfGrads;
  ValueType* FirstAddressOfdpsiM;
  ValueType* LastAddressOfdpsiM;

  /// determinant ratios with respect to the reference determinant
  ValueVector ratios_to_ref_;
  /// new determinant ratios with respect to the updated reference determinant upon a proposed move
  ValueVector new_ratios_to_ref_;
  /// new value of the reference determinant over the old value upon a proposed move
  ValueType curRatio;
  /// log value of the reference determinant
  LogValueType log_value_ref_det_;
  /// store determinant grads (old and new)
  GradMatrix grads, new_grads;
  /// store determinant lapls (old and new)
  ValueMatrix lapls, new_lapls;
  // additional storage for spin derivatives. Only resized if the calculation uses spinors
  ValueMatrix spingrads, new_spingrads;


  /* mmorales:
   *  i decided to stored the excitation information of all determinants in the following
   *  compact form: (all these integers are found consecutively in the array)
   *    For each determinant:
   *     -n : number of excitations
   *     -i1,i2,...,in : occupied orbital to be replaced (these must be numbers from 0:Nptcl-1)
   *     -a1,a2,...,an : excited states that replace the orbitals (these can be anything)
   */
  std::shared_ptr<std::vector<int>> detData;
  std::shared_ptr<std::vector<std::pair<int, int>>> uniquePairs;
  std::shared_ptr<std::vector<RealType>> DetSigns;
  MultiDiracDeterminantCalculator<ValueType> DetCalculator;
};


} // namespace qmcplusplus
#endif
