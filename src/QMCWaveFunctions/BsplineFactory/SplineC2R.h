//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file SplineC2R.h
 *
 * class to handle complex splines to real orbitals with splines of arbitrary precision
 * splines storage and computation is offloaded to accelerators using OpenMP target
 */
#ifndef QMCPLUSPLUS_SPLINE_C2R_OMPTARGET_H
#define QMCPLUSPLUS_SPLINE_C2R_OMPTARGET_H

#include <cstddef>
#include <memory>
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBsplineOffloadMapper.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "Utilities/TimerManager.h"
#include <ResourceHandle.h>
#include "SplineOMPTargetMultiWalkerMem.h"
#include <cstdint>

namespace qmcplusplus
{
template<typename T>
class MultiBsplineOffloadMapper;

/** class to match std::complex<ST> spline with BsplineSet::ValueType (real) SPOs with OpenMP offload
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of phase vectors
 * The internal storage of complex spline coefficients uses double sized real arrays of ST type, aligned and padded.
 * Calling assign_v assign_vgl should be restricted to the actual number of complex splines (kPoints.size()).
 * The first nComplexBands complex splines produce 2 real orbitals.
 * The rest complex splines produce 1 real orbital.
 * All the output orbitals are real (C2R). The maximal number of output orbitals is OrbitalSetSize.
 */
template<typename ST>
class SplineC2R : public BsplineSet
{
public:
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using TT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector;
  using BsplineSet::GradVector;
  using BsplineSet::HessVector;
  using BsplineSet::ValueVector;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

  template<typename DT>
  using OffloadVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using OffloadPosVector = VectorSoaContainer<DT, 3, OffloadAllocator<DT>>;

private:
  /// timer for offload portion
  NewTimer& offload_timer_;
  ///number of complex bands
  int nComplexBands;

  std::shared_ptr<OffloadVector<ST>> mKK_offload;
  std::shared_ptr<OffloadPosVector<ST>> myKcart_offload;
  const std::shared_ptr<OffloadVector<ST>> GGt_offload;
  const std::shared_ptr<OffloadVector<ST>> prim_lattice_G_offload;

  ResourceHandle<SplineOMPTargetMultiWalkerMem<ST, TT>> mw_mem_handle_;

  ///team private ratios for reduction, numVP x numTeams
  Matrix<TT, OffloadPinnedAllocator<TT>> ratios_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<TT, OffloadPinnedAllocator<TT>> results_scratch;
  ///psiinv and position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<TT, OffloadPinnedAllocator<TT>> psiinv_pos_copy;
  ///position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ST, OffloadPinnedAllocator<ST>> multi_pos_copy;

  void evaluateVGLMultiPos(const Vector<ST, OffloadPinnedAllocator<ST>>& multi_pos_copy,
                           Vector<ST, OffloadPinnedAllocator<ST>>& offload_scratch,
                           Vector<TT, OffloadPinnedAllocator<TT>>& results_scratch,
                           const RefVector<ValueVector>& psi_v_list,
                           const RefVector<GradVector>& dpsi_v_list,
                           const RefVector<ValueVector>& d2psi_v_list) const;

protected:
  ///multi bspline set
  const std::shared_ptr<MultiBsplineBase<ST>> SplineInst;
  /// multi bspline set offload mapper
  const std::shared_ptr<MultiBsplineOffloadMapper<ST>> offload_mapper_;
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2R(const std::string& my_name,
                     size_t size,
                     const Lattice& prim_lattice,
                     std::unique_ptr<MultiBsplineBase<ST>>&& multi_spline,
                     bool use_offload = false);

  SplineC2R(const SplineC2R& in);

  virtual std::string getClassName() const override { return "SplineC2R"; }
  virtual std::string getKeyword() const override { return "SplineC2R"; }
  virtual bool isOMPoffload() const override { return bool(offload_mapper_); }

  void createResource(ResourceCollection& collection) const override
  { auto resource_index = collection.addResource(std::make_unique<SplineOMPTargetMultiWalkerMem<ST, TT>>()); }

  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override;

  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const override;

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<SplineC2R>(*this); }

  void resizeStorage(size_t n) override;

  /// this routine can not be called from threaded region
  void finalizeConstruction() override;

  /** remap kPoints to pack the double copy */
  void resize_kpoints() override;

  void assign_v(const PointType& r, const vContainer_type& myV, ValueVector& psi, int first, int last) const;

  virtual void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override;

  virtual void evaluateDetRatios(const VirtualParticleSet& VP,
                                 ValueVector& psi,
                                 const ValueVector& psiinv,
                                 std::vector<ValueType>& ratios) override;

  virtual void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                    const RefVector<ValueVector>& psi_list,
                                    const std::vector<const ValueType*>& invRow_ptr_list,
                                    std::vector<std::vector<ValueType>>& ratios_list) const override;

  /** assign_vgl
   */
  void assign_vgl(const PointType& r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi, int first, int last)
      const;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  virtual void evaluateVGL(const ParticleSet& P,
                           const int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           ValueVector& d2psi) override;

  virtual void mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& sa_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector>& psi_v_list,
                              const RefVector<GradVector>& dpsi_v_list,
                              const RefVector<ValueVector>& d2psi_v_list) const override;

  virtual void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                              const RefVectorWithLeader<ParticleSet>& P_list,
                                              int iat,
                                              const std::vector<const ValueType*>& invRow_ptr_list,
                                              OffloadMWVGLArray& phi_vgl_v,
                                              std::vector<ValueType>& ratios,
                                              std::vector<GradType>& grads) const override;

  void assign_vgh(const PointType& r,
                  ValueVector& psi,
                  GradVector& dpsi,
                  HessVector& grad_grad_psi,
                  int first,
                  int last) const;

  virtual void evaluateVGH(const ParticleSet& P,
                           const int iat,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& grad_grad_psi) override;

  void assign_vghgh(const PointType& r,
                    ValueVector& psi,
                    GradVector& dpsi,
                    HessVector& grad_grad_psi,
                    GGGVector& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  virtual void evaluateVGHGH(const ParticleSet& P,
                             const int iat,
                             ValueVector& psi,
                             GradVector& dpsi,
                             HessVector& grad_grad_psi,
                             GGGVector& grad_grad_grad_psi) override;

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix& logdet,
                                    GradMatrix& dlogdet,
                                    ValueMatrix& d2logdet) override;

  friend class SplineSetReader<ST>;
  friend class BsplineReader;
};

extern template class SplineC2R<float>;
extern template class SplineC2R<double>;

} // namespace qmcplusplus
#endif
