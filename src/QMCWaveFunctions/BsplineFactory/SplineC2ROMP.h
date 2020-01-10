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
#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>
#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Utilities/FairDivide.h"

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

private:
  ///primitive cell
  CrystalLattice<ST, 3> PrimLattice;
  ///\f$GGt=G^t G \f$, transformation for tensor in LatticeUnit to CartesianUnit, e.g. Hessian
  Tensor<ST, 3> GGt;
  ///number of complex bands
  int nComplexBands;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST, OffloadAllocator<ST>>> SplineInst;

  vContainer_type mKK;
  VectorSoaContainer<ST, 3> myKcart;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<TT, OffloadPinnedAllocator<TT>> ratios_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<TT, OffloadPinnedAllocator<TT>> results_scratch;
  ///psiinv and position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<TT, OffloadPinnedAllocator<TT>> psiinv_pos_copy;
  ///position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<ST, OffloadPinnedAllocator<ST>> mw_pos_copy;
  ///the following pointers are used for keep and access the data on device
  ///cloned objects copy the pointer by value without the need of mapping to the device
  ///Thus master_PrimLattice_G_ptr is different from PrimLattice.G.data() in cloned objects
  ///mKK data pointer
  const ST* master_mKK_ptr;
  ///myKcart data pointer
  const ST* master_myKcart_ptr;
  ///PrimLattice.G data pointer
  const ST* master_PrimLattice_G_ptr;
  ///GGt data pointer
  const ST* master_GGt_ptr;

protected:
  /// intermediate result vectors
  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

public:
  SplineC2ROMP() : nComplexBands(0)
  {
    is_complex = true;
    className  = "SplineC2ROMP";
    KeyWord    = "SplineC2R";
  }

  ~SplineC2ROMP()
  {
    if (SplineInst.use_count() == 1)
    {
      // clean up mapping by the last owner
      const auto* MultiSpline = SplineInst->getSplinePtr();
      PRAGMA_OFFLOAD("omp target exit data map(delete:MultiSpline[0:1])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_mKK_ptr[0:mKK.size()])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_myKcart_ptr[0:myKcart.capacity()*3])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_PrimLattice_G_ptr[0:9])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_GGt_ptr[0:9])")
    }
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
    SplineInst = std::make_shared<MultiBspline<ST, OffloadAllocator<ST>>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  /// this routine can not be called from threaded region
  void finalizeConstruction() override
  {
    // map the SplineInst->getSplinePtr() structure to GPU
    auto* MultiSpline = SplineInst->getSplinePtr();
    PRAGMA_OFFLOAD("omp target enter data map(alloc:MultiSpline[0:1])")
    auto* restrict coefs = MultiSpline->coefs;
    // attach pointers on the device to achieve deep copy
    PRAGMA_OFFLOAD("omp target map(always, to: MultiSpline[0:1], coefs[0:MultiSpline->coefs_size])")
    {
      MultiSpline->coefs = coefs;
    }

    // transfer static data to GPU
    master_mKK_ptr = mKK.data();
    PRAGMA_OFFLOAD("omp target enter data map(alloc:master_mKK_ptr[0:mKK.size()])")
    PRAGMA_OFFLOAD("omp target update to(master_mKK_ptr[0:mKK.size()])")
    master_myKcart_ptr = myKcart.data();
    PRAGMA_OFFLOAD("omp target enter data map(alloc:master_myKcart_ptr[0:myKcart.capacity()*3])")
    PRAGMA_OFFLOAD("omp target update to(master_myKcart_ptr[0:myKcart.capacity()*3])")
    master_PrimLattice_G_ptr = PrimLattice.G.data();
    PRAGMA_OFFLOAD("omp target enter data map(alloc:master_PrimLattice_G_ptr[0:9])")
    PRAGMA_OFFLOAD("omp target update to(master_PrimLattice_G_ptr[0:9])")
    master_GGt_ptr = GGt.data();
    PRAGMA_OFFLOAD("omp target enter data map(alloc:master_GGt_ptr[0:9])")
    PRAGMA_OFFLOAD("omp target update to(master_GGt_ptr[0:9])")
    /* debug pointers
    std::cout << "Ye debug mapping" << std::endl;
    std::cout << "SplineInst = " << SplineInst << std::endl;
    std::cout << "MultiSpline = " << MultiSpline << std::endl;
    std::cout << "master_mKK_ptr = " << master_mKK_ptr << std::endl;
    std::cout << "master_myKcart_ptr = " << master_myKcart_ptr << std::endl;
    std::cout << "master_PrimLattice_G_ptr = " << master_PrimLattice_G_ptr << std::endl;
    std::cout << "master_GGt_ptr = " << master_GGt_ptr << std::endl;
    std::cout << "this = " << this << std::endl;
    */
  }

  inline void flush_zero() { SplineInst->flush_zero(); }

  /** remap kPoints to pack the double copy */
  inline void resize_kpoints()
  {
#ifndef QMC_CUDA
    // GPU CUDA code doesn't allow a change of the ordering
    nComplexBands = this->remap_kpoints();
#endif
    int nk = kPoints.size();
    mKK.resize(nk);
    myKcart.resize(nk);
    for (size_t i = 0; i < nk; ++i)
    {
      mKK[i]     = -dot(kPoints[i], kPoints[i]);
      myKcart(i) = kPoints[i];
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

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  void assign_vgl_from_l(const PointType& r, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void evaluateVGL(const ParticleSet& P,
                           const int iat,
                           ValueVector_t& psi,
                           GradVector_t& dpsi,
                           ValueVector_t& d2psi) override;

  virtual void mw_evaluateVGL(const std::vector<SPOSet*>& sa_list,
                              const std::vector<ParticleSet*>& P_list,
                              int iat,
                              const std::vector<ValueVector_t*>& psi_v_list,
                              const std::vector<GradVector_t*>& dpsi_v_list,
                              const std::vector<ValueVector_t*>& d2psi_v_list) override;

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

  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend class BsplineReaderBase;
};

extern template class SplineC2ROMP<float>;
extern template class SplineC2ROMP<double>;

} // namespace qmcplusplus
#endif
