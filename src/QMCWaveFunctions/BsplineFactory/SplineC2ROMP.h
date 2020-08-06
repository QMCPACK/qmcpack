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


/** @file SplineC2ROMP.h
 *
 * class to handle complex splines to real orbitals with splines of arbitrary precision
 * splines storage and computation is offloaded to accelerators using OpenMP target
 */
#ifndef QMCPLUSPLUS_SPLINE_C2R_OMP_H
#define QMCPLUSPLUS_SPLINE_C2R_OMP_H

#include <memory>
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include <OhmmsSoA/VectorSoaContainer.h>
#include <spline2/MultiBspline.hpp>
#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Utilities/FairDivide.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{
/** class to match std::complex<ST> spline with BsplineSet::ValueType (real) SPOs
 * @tparam ST precision of spline
 *
 * Requires temporage storage and multiplication of phase vectors
 * Internal storage use double sized arrays of ST type, aligned and padded.
 */
template<typename ST>
class SplineC2ROMP : public BsplineSet
{
public:
  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;

  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = TinyVector<ST, 3>;
  using SingleSplineType = UBspline_3d_d;
  // types for evaluation results
  using TT = typename BsplineSet::ValueType;
  using BsplineSet::GGGVector_t;
  using BsplineSet::GradVector_t;
  using BsplineSet::HessVector_t;
  using BsplineSet::ValueVector_t;

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
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///number of complex bands
  int nComplexBands;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST, OffloadAllocator<ST>, OffloadAllocator<SplineType>>> SplineInst;

  std::shared_ptr<OffloadVector<ST>> mKK;
  std::shared_ptr<OffloadPosVector<ST>> myKcart;
  std::shared_ptr<OffloadVector<ST>> GGt_offload;
  std::shared_ptr<OffloadVector<ST>> PrimLattice_G_offload;

  ///team private ratios for reduction, numVP x numTeams
  Matrix<TT, OffloadPinnedAllocator<TT>> ratios_private;
  ///team private ratios and grads for reduction, numVP x numTeams
  Matrix<TT, OffloadPinnedAllocator<TT>> rg_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<TT, OffloadPinnedAllocator<TT>> results_scratch;
  ///psiinv and position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<TT, OffloadPinnedAllocator<TT>> psiinv_pos_copy;
  ///psiinv and position scratch space of multiple walkers, used to avoid allocation on the fly and faster transfer
  Vector<TT, OffloadPinnedAllocator<TT>> mw_psiinv_pos_copy;
  ///position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ST, OffloadPinnedAllocator<ST>> multi_pos_copy;
  ///multi purpose H2D buffer for mw_evaluateVGLandDetRatioGrads
  Matrix<char, OffloadPinnedAllocator<char>> buffer_H2D;
  ///multi purpose H2D buffer for mw_evaluateDetRatios
  Vector<char, OffloadPinnedAllocator<char>> det_ratios_buffer_H2D;

  void evaluateVGLMultiPos(const Vector<ST, OffloadPinnedAllocator<ST>>& multi_pos_copy,
                           const RefVector<ValueVector_t>& psi_v_list,
                           const RefVector<GradVector_t>& dpsi_v_list,
                           const RefVector<ValueVector_t>& d2psi_v_list);

protected:
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2ROMP()
      : BsplineSet(true),
        offload_timer_(*TimerManager.createTimer("SplineC2ROMP::offload", timer_level_fine)),
        nComplexBands(0),
        GGt_offload(std::make_shared<OffloadVector<ST>>(9)),
        PrimLattice_G_offload(std::make_shared<OffloadVector<ST>>(9))
  {
    is_complex = true;
    className  = "SplineC2ROMP";
    KeyWord    = "SplineC2R";
  }

  virtual SPOSet* makeClone() const override { return new SplineC2ROMP(*this); }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    init_base(n);
    size_t npad = getAlignedSize<ST>(2 * n);
    myV.resize(npad);
    myG.resize(npad);
    myL.resize(npad);
    myH.resize(npad);
    mygH.resize(npad);
  }

  void bcast_tables(Communicate* comm) { chunked_bcast(comm, SplineInst->getSplinePtr()); }

  void gather_tables(Communicate* comm)
  {
    if (comm->size() == 1)
      return;
    const int Nbands      = kPoints.size();
    const int Nbandgroups = comm->size();
    offset.resize(Nbandgroups + 1, 0);
    FairDivideLow(Nbands, Nbandgroups, offset);

    for (size_t ib = 0; ib < offset.size(); ib++)
      offset[ib] = offset[ib] * 2;
    gatherv(comm, SplineInst->getSplinePtr(), SplineInst->getSplinePtr()->z_stride, offset);
  }

  template<typename GT, typename BCT>
  void create_spline(GT& xyz_g, BCT& xyz_bc)
  {
    resize_kpoints();
    SplineInst = std::make_shared<MultiBspline<ST, OffloadAllocator<ST>, OffloadAllocator<SplineType>>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  /// this routine can not be called from threaded region
  void finalizeConstruction() override
  {
    // map the SplineInst->getSplinePtr() structure to GPU
    auto* MultiSpline    = SplineInst->getSplinePtr();
    auto* restrict coefs = MultiSpline->coefs;
    // attach pointers on the device to achieve deep copy
    PRAGMA_OFFLOAD("omp target map(always, to: MultiSpline[0:1], coefs[0:MultiSpline->coefs_size])")
    {
      MultiSpline->coefs = coefs;
    }

    // transfer static data to GPU
    auto* mKK_ptr = mKK->data();
    PRAGMA_OFFLOAD("omp target update to(mKK_ptr[0:mKK->size()])")
    auto* myKcart_ptr = myKcart->data();
    PRAGMA_OFFLOAD("omp target update to(myKcart_ptr[0:myKcart->capacity()*3])")
    for (size_t i = 0; i < 9; i++)
    {
      (*GGt_offload)[i]           = GGt[i];
      (*PrimLattice_G_offload)[i] = PrimLattice.G[i];
    }
    auto* PrimLattice_G_ptr = PrimLattice_G_offload->data();
    PRAGMA_OFFLOAD("omp target update to(PrimLattice_G_ptr[0:9])")
    auto* GGt_ptr = GGt_offload->data();
    PRAGMA_OFFLOAD("omp target update to(GGt_ptr[0:9])")
  }

  inline void flush_zero() { SplineInst->flush_zero(); }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
#ifndef QMC_CUDA
    // GPU CUDA code doesn't allow a change of the ordering
    nComplexBands = this->remap_kpoints();
#endif
    int nk  = kPoints.size();
    mKK     = std::make_shared<OffloadVector<ST>>(nk);
    myKcart = std::make_shared<OffloadPosVector<ST>>(nk);
    for (size_t i = 0; i < nk; ++i)
    {
      (*mKK)[i]     = -dot(kPoints[i], kPoints[i]);
      (*myKcart)(i) = kPoints[i];
    }
  }

  void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level);

  bool read_splines(hdf_archive& h5f);

  bool write_splines(hdf_archive& h5f);

  void assign_v(const PointType& r, const vContainer_type& myV, ValueVector_t& psi, int first, int last) const;

  virtual void evaluateValue(const ParticleSet& P, const int iat, ValueVector_t& psi) override;

  virtual void evaluateDetRatios(const VirtualParticleSet& VP,
                                 ValueVector_t& psi,
                                 const ValueVector_t& psiinv,
                                 std::vector<ValueType>& ratios) override;

  virtual void mw_evaluateDetRatios(const RefVector<SPOSet>& spo_list,
                                    const RefVector<const VirtualParticleSet>& vp_list,
                                    const RefVector<ValueVector_t>& psi_list,
                                    const std::vector<const ValueType*>& invRow_ptr_list,
                                    std::vector<std::vector<ValueType>>& ratios_list) override;

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void evaluateVGL(const ParticleSet& P,
                           const int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           ValueVector_t& d2psi) override;

  virtual void mw_evaluateVGL(const RefVector<SPOSet>& sa_list,
                              const RefVector<ParticleSet>& P_list,
                              int iat,
                              const RefVector<ValueVector_t>& psi_v_list,
                              const RefVector<GradVector_t>& dpsi_v_list,
                              const RefVector<ValueVector_t>& d2psi_v_list) override;

  virtual void mw_evaluateVGLandDetRatioGrads(const RefVector<SPOSet>& spo_list,
                                              const RefVector<ParticleSet>& P_list,
                                              int iat,
                                              const std::vector<const ValueType*>& invRow_ptr_list,
                                              VGLVector_t& phi_vgl_v,
                                              std::vector<ValueType>& ratios,
                                              std::vector<GradType>& grads) override;

  void assign_vgh(const PointType& r,
                  ValueVector_t& psi,
                  GradVector_t& dpsi,
                  HessVector_t& grad_grad_psi,
                  int first,
                  int last) const;

  virtual void evaluateVGH(const ParticleSet& P,
                           const int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           HessVector_t& grad_grad_psi) override;

  void assign_vghgh(const PointType& r,
                    ValueVector_t& psi,
                    GradVector_t& dpsi,
                    HessVector_t& grad_grad_psi,
                    GGGVector_t& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const;

  virtual void evaluateVGHGH(const ParticleSet& P,
                             const int iat,
                             ValueVector_t& psi,
                             GradVector_t& dpsi,
                             HessVector_t& grad_grad_psi,
                             GGGVector_t& grad_grad_grad_psi) override;

  virtual void evaluate_notranspose(const ParticleSet& P,
                                    int first,
                                    int last,
                                    ValueMatrix_t& logdet,
                                    GradMatrix_t& dlogdet,
                                    ValueMatrix_t& d2logdet) override;

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend class BsplineReaderBase;
};

extern template class SplineC2ROMP<float>;
extern template class SplineC2ROMP<double>;

} // namespace qmcplusplus
#endif
