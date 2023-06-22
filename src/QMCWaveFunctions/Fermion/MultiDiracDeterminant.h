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
#include "QMCWaveFunctions/Fermion/SmallMatrixDetCalculator.h"
#include "Message/Communicate.h"
#include "Numerics/DeterminantOperators.h"
#include "ResourceCollection.h"
//#include "CPU/BLAS.hpp"

namespace qmcplusplus
{
class MultiDiracDeterminant : public WaveFunctionComponent
{
public:
  NewTimer &inverse_timer, &buildTable_timer, &table2ratios_timer, &evalWalker_timer, &evalOrbValue_timer,
      &evalOrbVGL_timer;
  NewTimer &updateInverse_timer, &calculateRatios_timer, &calculateGradRatios_timer, &updateRatios_timer;
  NewTimer &evaluateDetsForPtclMove_timer, &evaluateDetsAndGradsForPtclMove_timer, &evaluateGrads_timer;
  NewTimer &offload_timer, &transferH2D_timer, &transferD2H_timer;

  // Optimizable parameter
  opt_variables_type myVars;

  template<typename DT>
  using OffloadVector = Vector<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using UnpinnedOffloadVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using UnpinnedOffloadMatrix = Matrix<DT, OffloadAllocator<DT>>;

  using ValueVector = SPOSet::ValueVector;
  using ValueMatrix = SPOSet::ValueMatrix;

  struct MultiDiracDetMultiWalkerResource : public Resource
  {
    MultiDiracDetMultiWalkerResource() : Resource("MultiDiracDeterminant") {}
    MultiDiracDetMultiWalkerResource(const MultiDiracDetMultiWalkerResource&) : MultiDiracDetMultiWalkerResource() {}

    std::unique_ptr<Resource> makeClone() const override
    {
      return std::make_unique<MultiDiracDetMultiWalkerResource>(*this);
    }

    void resizeConstants(size_t nw)
    {
      if (nw > czero_vec.size())
      {
        czero_vec.resize(nw);
        cone_vec.resize(nw);
        cminus_one_vec.resize(nw);
        std::fill_n(czero_vec.data(), nw, 0);
        std::fill_n(cone_vec.data(), nw, 1);
        std::fill_n(cminus_one_vec.data(), nw, -1);
        czero_vec.updateTo();
        cone_vec.updateTo();
        cminus_one_vec.updateTo();
      }
      else
      {
        czero_vec.resize(nw);
        cone_vec.resize(nw);
        cminus_one_vec.resize(nw);
      }
    }

    OffloadVector<ValueType> czero_vec;
    OffloadVector<ValueType> cone_vec;
    OffloadVector<ValueType> cminus_one_vec;

    OffloadVector<ValueType*> workV1_deviceptr_list;
    OffloadVector<ValueType*> workV2_deviceptr_list;
    OffloadVector<ValueType*> psiV_temp_deviceptr_list;
    OffloadVector<ValueType*> psiMinv_temp_deviceptr_list;
    OffloadVector<ValueType*> dpsiMinv_deviceptr_list;

    OffloadVector<ValueType*> psiV_deviceptr_list;
    OffloadVector<GradType*> dpsiV_deviceptr_list;
    OffloadVector<ValueType*> TpsiM_deviceptr_list;
    OffloadVector<ValueType*> psiM_deviceptr_list;
    OffloadVector<ValueType*> psiMinv_deviceptr_list;
    OffloadVector<GradType*> dpsiM_deviceptr_list;

    // pointer lists used by mw_buildTableMatrix_calculateRatios_impl
    OffloadVector<ValueType*> psiinv_deviceptr_list;
    OffloadVector<ValueType*> psi_deviceptr_list;
    OffloadVector<ValueType*> table_matrix_deviceptr_list;
    OffloadVector<ValueType*> ratios_deviceptr_list;

    OffloadVector<ValueType> curRatio_list;
    OffloadVector<ValueType> inv_curRatio_list;

    OffloadVector<ValueType> det0_grad_list;
    OffloadVector<GradType> ratioGradRef_list;
  };

  //lookup table mapping the unique determinants to their element position in C2_node vector
  std::vector<std::vector<int>> lookup_tbl;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  MultiDiracDeterminant(std::unique_ptr<SPOSet>&& spos, bool spinor, int first, int nel);

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

  std::string getClassName() const override { return "MultiDiracDeterminant"; }

  bool isFermionic() const final { return true; }
  inline bool isOptimizable() const final { return Phi->isOptimizable(); }

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) final
  {
    Phi->extractOptimizableObjectRefs(opt_obj_refs);
  }

  inline void checkOutVariables(const opt_variables_type& active) override
  {
    if (Phi->isOptimizable())
      Phi->checkOutVariables(active);
  }

  /// create optimizable orbital rotation parameters
  void buildOptVariables(std::vector<size_t>& C2node);
  ///helper function to buildOptVariables
  int build_occ_vec(const OffloadVector<int>& data, const size_t nel, const size_t nmo, std::vector<int>& occ_vec);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override
  {}

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi,
                           const MultiDiracDeterminant& pseudo_dn,
                           const ValueType& psiCurrent,
                           const std::vector<ValueType>& Coeff,
                           const std::vector<size_t>& C2node_up,
                           const std::vector<size_t>& C2node_dn);

  void evaluateDerivativesWF(ParticleSet& P,
                             const opt_variables_type& optvars,
                             Vector<ValueType>& dlogpsi,
                             const MultiDiracDeterminant& pseudo_dn,
                             const PsiValueType& psiCurrent,
                             const std::vector<ValueType>& Coeff,
                             const std::vector<size_t>& C2node_up,
                             const std::vector<size_t>& C2node_dn);


  void registerData(ParticleSet& P, WFBufferType& buf) override;

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat) override;

  static void mw_accept_rejectMove(const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list,
                                   const RefVectorWithLeader<ParticleSet>& p_list,
                                   int iat,
                                   const std::vector<bool>& isAccepted);

  void createResource(ResourceCollection& collection) const override;
  void acquireResource(ResourceCollection& collection,
                       const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list) const;
  void releaseResource(ResourceCollection& collection,
                       const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list) const;

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

  /** create necessary structures related to unique determinants
   * sort configlist_unsorted by excitation level abd store the results in ciConfigList (class member)
   * ciConfigList shouldn't change during a simulation after it is sorted here
   *
   * @param ref_det_id id of the reference determinant before sorting
   * @param configlist_unsorted config list to be loaded.
   * @param C2nodes_unsorted mapping from overall det index to unique det (configlist_unsorted) index
   * @param C2nodes_sorted mapping from overall det index to unique det (ciConfigList) index
   */
  void createDetData(const int ref_det_id,
                     const std::vector<ci_configuration2>& configlist_unsorted,
                     const std::vector<size_t>& C2nodes_unsorted,
                     std::vector<size_t>& C2nodes_sorted);

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
                                                 int iat,
                                                 UnpinnedOffloadMatrix<ValueType>& mw_grads);
  /// evaluate the value and gradients of all the unique determinants with one electron moved. Used by the table method. Includes Spin Gradient data
  void evaluateDetsAndGradsForPtclMoveWithSpin(const ParticleSet& P, int iat);


  /// evaluate the gradients of all the unique determinants with one electron moved. Used by the table method
  void evaluateGrads(ParticleSet& P, int iat);
  /// multi walker version of mw_evaluateGrads
  void static mw_evaluateGrads(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                               const RefVectorWithLeader<ParticleSet>& P_list,
                               int iat,
                               UnpinnedOffloadMatrix<ValueType>& mw_grads);
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

  const OffloadVector<ValueType>& getRatiosToRefDet() const { return ratios_to_ref_; }
  const OffloadVector<ValueType>& getNewRatiosToRefDet() const { return new_ratios_to_ref_; }
  const Matrix<GradType>& getGrads() const { return grads; }
  const Matrix<GradType>& getNewGrads() const { return new_grads; }
  const Matrix<ValueType>& getLapls() const { return lapls; }
  const Matrix<ValueType>& getNewLapls() const { return new_lapls; }
  const Matrix<ValueType>& getSpinGrads() const { return spingrads; }
  const Matrix<ValueType>& getNewSpinGrads() const { return new_spingrads; }

  PsiValueType getRefDetRatio() const { return static_cast<PsiValueType>(curRatio); }
  LogValueType getLogValueRefDet() const { return log_value_ref_det_; }

private:
  void mw_InverseUpdateByColumn(MultiDiracDetMultiWalkerResource& mw_res,
                                const int working_index,
                                const OffloadVector<ValueType>& curRatio_list,
                                const OffloadVector<ValueType*>& psiV_deviceptr_list,
                                const OffloadVector<ValueType*>& psiMinv_deviceptr_list,
                                const size_t psiMinv_rows) const;

  /** update ratios with respect to the reference deteriminant for a given excitation level
   * @param ext_level excitation level
   * @param det_offset offset of the determinant id
   * @param data_offset offset of the "data" structure
   * @param sign of determinants
   * @param table_matrix_list list of table_matrix
   *
   * this is a general implementation. Support abitrary excitation level
   */
  void mw_updateRatios_generic(const int ext_level,
                               const size_t det_offset,
                               const size_t data_offset,
                               SmallMatrixDetCalculator<ValueType>& det_calculator,
                               const OffloadVector<int>& data,
                               const OffloadVector<RealType>& sign,
                               const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
                               const RefVector<OffloadVector<ValueType>>& ratios_list) const;

  /** update ratios of the reference deteriminant
   * @param det0_list list of reference det value
   */
  void mw_updateRatios_det0(const OffloadVector<ValueType>& det0_list,
                            const OffloadVector<ValueType*>& ratios_deviceptr_list) const;

  /** update ratios with respect to the reference deteriminant for a given excitation level
   * @param det_offset offset of the determinant id
   * @param data_offset offset of the "data" structure
   * @param sign of determinants
   * @param table_matrix_list list of table_matrix
   *
   * this is intended to be customized based on EXT_LEVEL
   */
  template<unsigned EXT_LEVEL>
  void mw_updateRatios(const size_t det_offset,
                       const size_t data_offset,
                       const OffloadVector<int>& data,
                       const OffloadVector<RealType>& sign,
                       const OffloadVector<ValueType*>& table_matrix_deviceptr_list,
                       const size_t num_table_matrix_cols,
                       const OffloadVector<ValueType*>& ratios_deviceptr_list) const;

  /** Function to calculate the ratio of the excited determinant to the reference determinant in CustomizedMatrixDet following the paper by Clark et al. JCP 135(24), 244105
   *@param nw Number of walkers in the batch
   *@param ref ID of the reference determinant
   *@param det0_list takes lists of ValueType(1) for the value or RatioGrad/curRatio for the gradients
   *@param psiinv_list
   *@param psi_list
   *@param data  (Shared by all determinants)
   *@param pairs is the number of unique determinants (std::pair[Nb_unique_alpha][Nb_unique_beta]) (Shared by all determinants)
   *@param sign (Shared by all determinants)
   *@param table_matrix_list stores all the dot products between 2 determinants (I,J)
   *@param ratio_list returned computed ratios
   */
  void mw_buildTableMatrix_calculateRatios_impl(MultiDiracDetMultiWalkerResource& mw_res,
                                                int ref,
                                                const OffloadVector<ValueType>& det0_list,
                                                const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
                                                const RefVector<OffloadMatrix<ValueType>>& psi_list,
                                                const OffloadVector<int>& data,
                                                const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                                const OffloadVector<RealType>& sign,
                                                const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
                                                const RefVector<OffloadVector<ValueType>>& ratios_list);

  /** Function to calculate the ratio of the excited determinant to the reference determinant in CustomizedMatrixDet following the paper by Clark et al. JCP 135(24), 244105
   *@param ref ID of the reference determinant
   *@param det0 take ValueType(1) for the value or RatioGrad/curRatio for the gradients
   *@param ratios returned computed ratios
   *@param psiinv
   *@param psi
   *@param table_matrix stores all the dot products between 2 determinants (I,J)
   *@param data  (Shared by all determinants)
   *@param pairs is the number of unique determinants (std::pair[Nb_unique_alpha][Nb_unique_beta]) (Shared by all determinants)
   *@param sign (Shared by all determinants)
   */
  void buildTableMatrix_calculateRatios_impl(int ref,
                                             ValueType det0,
                                             ValueType* restrict ratios,
                                             const OffloadMatrix<ValueType>& psiinv,
                                             const OffloadMatrix<ValueType>& psi,
                                             OffloadMatrix<ValueType>& table_matrix,
                                             const OffloadVector<int>& data,
                                             const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                             const OffloadVector<RealType>& sign);

  /** compute the ratio of the excited determinants to the reference determinant
   * @param ratios the output.
   */
  void buildTableMatrix_calculateRatios(int ref,
                                        const OffloadMatrix<ValueType>& psiinv,
                                        const OffloadMatrix<ValueType>& psi,
                                        const OffloadVector<int>& data,
                                        const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                        const OffloadVector<RealType>& sign,
                                        OffloadMatrix<ValueType>& table_matrix,
                                        OffloadVector<ValueType>& ratios);

  void mw_buildTableMatrix_calculateRatios(MultiDiracDetMultiWalkerResource& mw_res,
                                           int ref,
                                           const OffloadVector<ValueType>& det0_list,
                                           const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
                                           const RefVector<OffloadMatrix<ValueType>>& psi_list,
                                           const OffloadVector<int>& data,
                                           const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                           const OffloadVector<RealType>& sign,
                                           const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
                                           const RefVector<OffloadVector<ValueType>>& ratios_list);

  /** Function to calculate the ratio of the gradients of the excited determinant to the reference determinant in CustomizedMatrixDet following the paper by Clark et al. JCP 135(24), 244105
   *@param ref ID of the reference determinant
   *@param psiinv
   *@param psi
   *@param data  (Shared by all determinants)
   *@param pairs is the number of unique determinants (std::pair[Nb_unique_alpha][Nb_unique_beta]) (Shared by all determinants)
   *@param sign (Shared by all determinants)
   *@param det0_grad gradient value taking RatioGrad/curRatio 
   *@param table_matrix stores all the dot products between 2 determinants (I,J)
   *@param dx dimension (OHMMS_DIM)
   *@param iat atom ID 
   *@param grads returned computed gradients
   */
  void buildTableMatrix_calculateGradRatios(int ref,
                                            const OffloadMatrix<ValueType>& psiinv,
                                            const OffloadMatrix<ValueType>& psi,
                                            const OffloadVector<int>& data,
                                            const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                            const OffloadVector<RealType>& sign,
                                            const ValueType& det0_grad,
                                            OffloadMatrix<ValueType>& table_matrix,
                                            int dx,
                                            int iat,
                                            Matrix<GradType>& grads);

  /** Function to calculate the ratio of the gradients of the excited determinant to the reference determinant in CustomizedMatrixDet following the paper by Clark et al. JCP 135(24), 244105
   *@param nw Number of walkers in the batch
   *@param ref ID of the reference determinant
   *@param iat atom ID 
   *@param dx dimension (OHMMS_DIM)
   *@param getNumDets Number of determinants
   *@param psiinv_list
   *@param psi_list
   *@param data  (Shared by all determinants)
   *@param pairs is the number of unique determinants (std::pair[Nb_unique_alpha][Nb_unique_beta]) (Shared by all determinants)
   *@param sign (Shared by all determinants)
   *@param WorkSpace_list list refering to det.WorkSpace  
   *@param table_matrix_list stores all the dot products between 2 determinants (I,J)
   *@param ratios_list returned computed list of gradients
   */
  void mw_buildTableMatrix_calculateGradRatios(MultiDiracDetMultiWalkerResource& mw_res,
                                               int ref,
                                               int iat,
                                               int dx,
                                               int getNumDets,
                                               const OffloadVector<ValueType>& det0_grad_list,
                                               const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
                                               const RefVector<OffloadMatrix<ValueType>>& psi_list,
                                               const OffloadVector<int>& data,
                                               const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
                                               const OffloadVector<RealType>& sign,
                                               const RefVector<OffloadVector<ValueType>>& WorkSpace_list,
                                               const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
                                               UnpinnedOffloadMatrix<ValueType>& mw_grads);

  void buildTableMatrix_calculateRatiosValueMatrixOneParticle(
      int ref,
      const OffloadMatrix<ValueType>& psiinv,
      const OffloadMatrix<ValueType>& psi,
      const OffloadVector<int>& data,
      const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
      const OffloadVector<RealType>& sign,
      OffloadMatrix<ValueType>& table_matrix,
      int iat,
      Matrix<ValueType>& ratios);

  ///reset the size: with the number of particles
  void resize();

  ///a set of single-particle orbitals used to fill in the  values of the matrix
  const std::unique_ptr<SPOSet> Phi;
  ///number of single-particle orbitals which belong to this Dirac determinant
  const int NumOrbitals;
  ///index of the first particle with respect to the particle set
  const int FirstIndex;
  ///number of particles which belong to this Dirac determinant
  const int NumPtcls;
  ///index of the last particle with respect to the particle set
  const int LastIndex;
  ///use shared_ptr
  std::shared_ptr<std::vector<ci_configuration2>> ciConfigList;
  /// all the unique determinants are sorted, the id of the reference det id is always 0
  static constexpr int ReferenceDeterminant = 0;
  /// reference determinant occupation
  std::shared_ptr<OffloadVector<size_t>> refdet_occup;
  // flag to determine if spin arrays need to be resized and used. Set by ParticleSet::is_spinor_ in SlaterDetBuilder
  const bool is_spinor_;

  /// psiM(i,j) \f$= \psi_j({\bf r}_i)\f$
  /// TpsiM(i,j) \f$= psiM(j,i) \f$
  OffloadMatrix<ValueType> psiM, TpsiM;
  /// inverse Dirac determinant matrix of the reference det
  OffloadMatrix<ValueType> psiMinv, psiMinv_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  OffloadMatrix<GradType> dpsiM;
  // temporaty storage
  OffloadMatrix<ValueType> dpsiMinv;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  OffloadMatrix<ValueType> d2psiM;
  /* dspin_psiM(i,j) \f$= \partial_{s_i} \psi_j({\bf r}_i,s_i)\f$ where \f$s_i\f$s is the spin variable
   * This is only resized if a spinor calculation is used
   */
  ValueMatrix dspin_psiM;

  /// value of single-particle orbital for particle-by-particle update
  //ValueVector psiV, psiV_temp;
  OffloadVector<ValueType> psiV, psiV_temp;
  OffloadVector<GradType> dpsiV;
  OffloadVector<ValueType> d2psiV;
  OffloadVector<ValueType> workV1, workV2;
  //spin  derivative of single-particle orbitals. Only resized if a spinor calculation
  ValueVector dspin_psiV;

  OffloadMatrix<ValueType> table_matrix;

  OffloadVector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueType* FirstAddressOfGrads;
  ValueType* LastAddressOfGrads;
  ValueType* FirstAddressOfdpsiM;
  ValueType* LastAddressOfdpsiM;

  /// determinant ratios with respect to the reference determinant
  OffloadVector<ValueType> ratios_to_ref_;
  /// new determinant ratios with respect to the updated reference determinant upon a proposed move
  OffloadVector<ValueType> new_ratios_to_ref_;
  /// new value of the reference determinant over the old value upon a proposed move
  ValueType curRatio;
  /// log value of the reference determinant
  LogValueType log_value_ref_det_;
  /// store determinant grads (old and new)
  Matrix<GradType> grads, new_grads;
  /// store determinant lapls (old and new)
  Matrix<ValueType> lapls, new_lapls;
  // additional storage for spin derivatives. Only resized if the calculation uses spinors
  Matrix<ValueType> spingrads, new_spingrads;


  /* mmorales:
   *  i decided to stored the excitation information of all determinants in the following
   *  compact form: (all these integers are found consecutively in the array)
   *    For each determinant:
   *     -n : number of excitations
   *     -i1,i2,...,in : occupied orbital to be replaced (these must be numbers from 0:Nptcl-1)
   *     -a1,a2,...,an : excited states that replace the orbitals (these can be anything)
   */
  std::shared_ptr<OffloadVector<int>> detData;
  std::shared_ptr<VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>> uniquePairs;
  std::shared_ptr<OffloadVector<RealType>> DetSigns;
  /** number of unique determinants at each excitation level (relative to reference)
   *  {1, n_singles, n_doubles, n_triples, ...}
   */
  std::shared_ptr<std::vector<int>> ndets_per_excitation_level_;
  SmallMatrixDetCalculator<ValueType> det_calculator_;

  /// for matrices with leading dimensions <= MaxSmallDet, compute determinant with direct expansion.
  static constexpr size_t MaxSmallDet = 5;

  ResourceHandle<MultiDiracDetMultiWalkerResource> mw_res_handle_;
};


} // namespace qmcplusplus
#endif
