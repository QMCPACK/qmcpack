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
 * Adoptor classes to handle complex-to-(real,complex) with arbitrary precision
 */
#ifndef QMCPLUSPLUS_EINSPLINE_C2R_OMP_H
#define QMCPLUSPLUS_EINSPLINE_C2R_OMP_H

#include <memory>
#include <OhmmsSoA/Container.h>
#include <spline2/MultiBspline.hpp>
#include <spline2/MultiBsplineEval.hpp>
#include <spline2/MultiBsplineEval_OMPoffload.hpp>
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorBase.h"
#include "OpenMP/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
namespace C2R
{
template<typename ST, typename TT>
inline void assign_v(ST x,
                     ST y,
                     ST z,
                     TT* restrict results_scratch_ptr,
                     size_t orb_size,
                     const ST* restrict offload_scratch_ptr,
                     const ST* restrict myKcart_ptr,
                     size_t myKcart_padded_size,
                     size_t first_spo,
                     int nComplexBands,
                     int first,
                     int last)
{
  // protect last
  if (last > orb_size)
    last = orb_size;

  const ST* restrict kx = myKcart_ptr;
  const ST* restrict ky = myKcart_ptr + myKcart_padded_size;
  const ST* restrict kz = myKcart_ptr + myKcart_padded_size * 2;

  const ST* restrict val = offload_scratch_ptr;
  TT* restrict psi_s     = results_scratch_ptr;

#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd
#endif
  for (size_t j = first; j < last; j++)
  {
    const size_t jr = j << 1;
    const size_t ji = jr + 1;
    //phase
    ST s, c, p = -(x * kx[j] + y * ky[j] + z * kz[j]);
    sincos(p, &s, &c);

    const ST val_r        = val[jr];
    const ST val_i        = val[ji];
    const size_t psiIndex = first_spo + j + (j < nComplexBands ? j : nComplexBands);
    psi_s[psiIndex]       = val_r * c - val_i * s;
    if (j < nComplexBands)
      psi_s[psiIndex + 1] = val_i * c + val_r * s;
  }
}

/** assign_vgl
   */
template<typename ST, typename TT>
inline void assign_vgl(ST x,
                       ST y,
                       ST z,
                       TT* restrict results_scratch_ptr,
                       const ST* mKK_ptr,
                       size_t orb_size,
                       const ST* restrict offload_scratch_ptr,
                       size_t spline_padded_size,
                       const ST symGGt[6],
                       const ST G[9],
                       const ST* myKcart_ptr,
                       size_t myKcart_padded_size,
                       size_t first_spo,
                       int nComplexBands,
                       int first,
                       int last)
{
  // protect last
  if (last > orb_size)
    last = orb_size;

  constexpr ST two(2);
  const ST &g00 = G[0], &g01 = G[1], &g02 = G[2],
           &g10 = G[3], &g11 = G[4], &g12 = G[5],
           &g20 = G[6], &g21 = G[7], &g22 = G[8];

  const ST* restrict k0 = myKcart_ptr;
  const ST* restrict k1 = myKcart_ptr + myKcart_padded_size;
  const ST* restrict k2 = myKcart_ptr + myKcart_padded_size * 2;

  const ST* restrict val = offload_scratch_ptr;
  const ST* restrict g0  = offload_scratch_ptr + spline_padded_size;
  const ST* restrict g1  = offload_scratch_ptr + spline_padded_size * 2;
  const ST* restrict g2  = offload_scratch_ptr + spline_padded_size * 3;
  const ST* restrict h00 = offload_scratch_ptr + spline_padded_size * 4;
  const ST* restrict h01 = offload_scratch_ptr + spline_padded_size * 5;
  const ST* restrict h02 = offload_scratch_ptr + spline_padded_size * 6;
  const ST* restrict h11 = offload_scratch_ptr + spline_padded_size * 7;
  const ST* restrict h12 = offload_scratch_ptr + spline_padded_size * 8;
  const ST* restrict h22 = offload_scratch_ptr + spline_padded_size * 9;

  TT* restrict psi   = results_scratch_ptr;
  TT* restrict dpsi  = results_scratch_ptr + orb_size;
  TT* restrict d2psi = results_scratch_ptr + orb_size * 4;
#ifdef ENABLE_OFFLOAD
#pragma omp for
#else
#pragma omp simd
#endif
  for (size_t j = first; j < last; j++)
  {
    const size_t jr = j << 1;
    const size_t ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = val[jr];
    const ST val_i = val[ji];

    //phase
    ST s, c, p = -(x * kX + y * kY + z * kZ);
    sincos(p, &s, &c);

    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
    const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
    const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

    const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
    const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
    const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

    // \f$\nabla \psi_r + {\bf k}\psi_i\f$
    const ST gX_r = dX_r + val_i * kX;
    const ST gY_r = dY_r + val_i * kY;
    const ST gZ_r = dZ_r + val_i * kZ;
    const ST gX_i = dX_i - val_r * kX;
    const ST gY_i = dY_i - val_r * kY;
    const ST gZ_i = dZ_i - val_r * kZ;

    const ST lcart_r = SymTrace(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], symGGt);
    const ST lcart_i = SymTrace(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], symGGt);
    const ST lap_r   = lcart_r + mKK_ptr[j] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
    const ST lap_i   = lcart_i + mKK_ptr[j] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);

    const size_t psiIndex = first_spo + j + (j < nComplexBands ? j : nComplexBands);
    //this will be fixed later
    psi[psiIndex]   = c * val_r - s * val_i;
    d2psi[psiIndex] = c * lap_r - s * lap_i;
    //this will go way with Determinant
    dpsi[psiIndex * 3]     = c * gX_r - s * gX_i;
    dpsi[psiIndex * 3 + 1] = c * gY_r - s * gY_i;
    dpsi[psiIndex * 3 + 2] = c * gZ_r - s * gZ_i;

    if (j < nComplexBands)
    {
      psi[psiIndex + 1]      = c * val_i + s * val_r;
      d2psi[psiIndex + 1]    = c * lap_i + s * lap_r;
      dpsi[psiIndex * 3 + 3] = c * gX_i + s * gX_r;
      dpsi[psiIndex * 3 + 4] = c * gY_i + s * gY_r;
      dpsi[psiIndex * 3 + 5] = c * gZ_i + s * gZ_r;
    }
  }
}
} // namespace C2R

/** adoptor class to match std::complex<ST> spline with TT real SPOs
 * @tparam ST precision of spline
 * @tparam TT precision of SPOs
 * @tparam D dimension
 *
 * Requires temporage storage and multiplication of phase vectors
 * Internal storage use double sized arrays of ST type, aligned and padded.
 */
template<typename ST, typename TT>
struct SplineC2ROMP : public SplineAdoptorBase<ST, 3>
{
  static const int ALIGN   = QMC_CLINE;
  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;

  static const int D     = 3;
  using Base             = SplineAdoptorBase<ST, 3>;
  using SplineType       = typename bspline_traits<ST, 3>::SplineType;
  using BCType           = typename bspline_traits<ST, 3>::BCType;
  using DataType         = ST;
  using PointType        = typename Base::PointType;
  using SingleSplineType = typename Base::SingleSplineType;

  using vContainer_type  = Vector<ST, aligned_allocator<ST>>;
  using gContainer_type  = VectorSoaContainer<ST, 3>;
  using hContainer_type  = VectorSoaContainer<ST, 6>;
  using ghContainer_type = VectorSoaContainer<ST, 10>;

  using Base::first_spo;
  using Base::GGt;
  using Base::kPoints;
  using Base::last_spo;
  using Base::MakeTwoCopies;
  using Base::offset;
  using Base::PrimLattice;

  ///number of complex bands
  int nComplexBands;
  ///multi bspline set
  std::shared_ptr<MultiBspline<ST, ALIGN, OffloadAllocator<ST>>> SplineInst;

  vContainer_type mKK;
  VectorSoaContainer<ST, 3> myKcart;

  vContainer_type myV;
  vContainer_type myL;
  gContainer_type myG;
  hContainer_type myH;
  ghContainer_type mygH;

  ///thread private ratios for reduction when using nested threading, numVP x numThread
  Matrix<TT, OffloadPinnedAllocator<TT>> ratios_private;
  ///offload scratch space, dynamically resized to the maximal need
  Vector<ST, OffloadPinnedAllocator<ST>> offload_scratch;
  ///result scratch space, dynamically resized to the maximal need
  Vector<TT, OffloadPinnedAllocator<TT>> results_scratch;
  ///psiinv and position scratch space, used to avoid allocation on the fly and faster transfer
  Vector<TT, OffloadPinnedAllocator<TT>> psiinv_pos_copy;
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

  SplineC2ROMP() : Base(), nComplexBands(0)
  {
    this->is_complex   = true;
    this->is_soa_ready = true;
    this->AdoptorName  = "SplineC2ROMPAdoptor";
    this->KeyWord      = "SplineC2ROMP";
  }

  ~SplineC2ROMP()
  {
    if (SplineInst.use_count() == 1)
    {
      // clean up mapping by the last owner
      const auto* MultiSpline =  SplineInst->getSplinePtr();
      PRAGMA_OFFLOAD("omp target exit data map(delete:MultiSpline[0:1])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_mKK_ptr[0:mKK.size()])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_myKcart_ptr[0:myKcart.capacity()*3])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_PrimLattice_G_ptr[0:9])")
      PRAGMA_OFFLOAD("omp target exit data map(delete:master_GGt_ptr[0:9])")
    }
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    Base::init_base(n);
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
    SplineInst = std::make_shared<MultiBspline<ST, ALIGN, OffloadAllocator<ST>>>();
    SplineInst->create(xyz_g, xyz_bc, myV.size());

    app_log() << "MEMORY " << SplineInst->sizeInByte() / (1 << 20) << " MB allocated "
              << "for the coefficients in 3D spline orbital representation" << std::endl;
  }

  /// this routine can not be called from threaded region
  void finalizeConstruction()
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

  inline void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int level)
  {
    SplineInst->copy_spline(spline_r, 2 * ispline);
    SplineInst->copy_spline(spline_i, 2 * ispline + 1);
  }

  bool read_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o << "spline_" << SplineAdoptorBase<ST, D>::MyIndex;
    einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
    return h5f.readEntry(bigtable, o.str().c_str()); //"spline_0");
  }

  bool write_splines(hdf_archive& h5f)
  {
    std::ostringstream o;
    o << "spline_" << SplineAdoptorBase<ST, D>::MyIndex;
    einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
    return h5f.writeEntry(bigtable, o.str().c_str()); //"spline_0");
  }

  template<typename VV>
  inline void assign_v(const PointType& r, const vContainer_type& myV, VV& psi, int first, int last) const
  {
    // protect last
    last = last > kPoints.size() ? kPoints.size() : last;

    const ST x = r[0], y = r[1], z = r[2];
    const ST* restrict kx = myKcart.data(0);
    const ST* restrict ky = myKcart.data(1);
    const ST* restrict kz = myKcart.data(2);

    TT* restrict psi_s = psi.data() + first_spo;
#pragma omp simd
    for (size_t j = first; j < std::min(nComplexBands, last); j++)
    {
      ST s, c;
      const size_t jr = j << 1;
      const size_t ji = jr + 1;
      const ST val_r  = myV[jr];
      const ST val_i  = myV[ji];
      sincos(-(x * kx[j] + y * ky[j] + z * kz[j]), &s, &c);
      psi_s[jr] = val_r * c - val_i * s;
      psi_s[ji] = val_i * c + val_r * s;
    }

    psi_s += nComplexBands;
#pragma omp simd
    for (size_t j = std::max(nComplexBands, first); j < last; j++)
    {
      ST s, c;
      const ST val_r = myV[2 * j];
      const ST val_i = myV[2 * j + 1];
      sincos(-(x * kx[j] + y * ky[j] + z * kz[j]), &s, &c);
      psi_s[j] = val_r * c - val_i * s;
    }
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r = P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));

    if (true)
    {
#pragma omp parallel
      {
        int first, last;
        FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

        spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
        assign_v(r, myV, psi, first / 2, last / 2);
      }
    }
    else
    {
      const int ChunkSizePerTeam = 128;
      const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

      const auto padded_size = myV.size();
      if (offload_scratch.size() < padded_size)
        offload_scratch.resize(padded_size);

      // Ye: need to extract sizes and pointers before entering target region
      const auto orb_size       = psi.size();
      const auto* spline_ptr    = SplineInst->getSplinePtr();
      auto* offload_scratch_ptr = offload_scratch.data();
      auto* psi_ptr             = psi.data();
      const auto x = r[0], y = r[1], z = r[2];
      const auto rux = ru[0], ruy = ru[1], ruz = ru[2];
      const auto myKcart_padded_size = myKcart.capacity();
      auto* myKcart_ptr              = master_myKcart_ptr;
      const size_t first_spo_local   = first_spo;
      const int nComplexBands_local  = nComplexBands;

      PRAGMA_OFFLOAD("omp target teams distribute num_teams(NumTeams) thread_limit(ChunkSizePerTeam) \
                  map(always, from: psi_ptr[0:orb_size])")
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, rux, ruy, ruz, ix, iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel")
        {
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, a, b, c, offload_scratch_ptr + first, first, last);
          C2R::assign_v(x, y, z, psi_ptr, orb_size, offload_scratch_ptr, myKcart_ptr, myKcart_padded_size,
                        first_spo_local, nComplexBands_local, first / 2, last / 2);
        }
      }
    }
  }

  template<typename VV, typename RT>
  inline void evaluateDetRatios(const VirtualParticleSet& VP, VV& psi, const VV& psiinv, std::vector<RT>& ratios)
  {
    const int nVP = VP.getTotalNum();
    if(psiinv_pos_copy.size() < psiinv.size() + nVP * 6)
      psiinv_pos_copy.resize(psiinv.size() + nVP * 6);

    // stage psiinv to psiinv_pos_copy
    std::copy_n(psiinv.data(), psiinv.size(), psiinv_pos_copy.data());

    // pack particle positions
    auto* restrict pos_scratch = psiinv_pos_copy.data() + psiinv.size();
    for (int iat = 0; iat < nVP; ++iat)
    {
      const PointType& r = VP.activeR(iat);
      PointType ru(PrimLattice.toUnit_floor(r));
      pos_scratch[iat * 6]     = r[0];
      pos_scratch[iat * 6 + 1] = r[1];
      pos_scratch[iat * 6 + 2] = r[2];
      pos_scratch[iat * 6 + 3] = ru[0];
      pos_scratch[iat * 6 + 4] = ru[1];
      pos_scratch[iat * 6 + 5] = ru[2];
    }

    const int ChunkSizePerTeam = 128;
    const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
    if (ratios_private.size() < NumTeams * nVP)
      ratios_private.resize(nVP, NumTeams);
    const auto padded_size = myV.size();
    if (offload_scratch.size() < padded_size * nVP)
      offload_scratch.resize(padded_size * nVP);
    const auto orb_size = psiinv.size();
    if (results_scratch.size() < orb_size * nVP)
      results_scratch.resize(orb_size * nVP);

    // Ye: need to extract sizes and pointers before entering target region
    const auto* spline_ptr         = SplineInst->getSplinePtr();
    auto* offload_scratch_ptr      = offload_scratch.data();
    auto* results_scratch_ptr      = results_scratch.data();
    const auto myKcart_padded_size = myKcart.capacity();
    auto* myKcart_ptr              = master_myKcart_ptr;
    auto* psiinv_ptr               = psiinv_pos_copy.data();
    auto* ratios_private_ptr       = ratios_private.data();
    const size_t first_spo_local   = first_spo;
    const int nComplexBands_local  = nComplexBands;

    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*nVP) thread_limit(ChunkSizePerTeam) \
                map(always, to: psiinv_ptr[0:psiinv_pos_copy.size()]) \
                map(always, from: ratios_private_ptr[0:NumTeams*nVP])")
    for (int iat = 0; iat < nVP; iat++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first      = ChunkSizePerTeam * team_id;
        const int last       = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;
        const int first_cplx = first / 2;
        const int last_cplx  = orb_size < last / 2 ? orb_size : last / 2;
        const int first_real = first_cplx + std::min(nComplexBands_local, first_cplx);
        const int last_real  = last_cplx + std::min(nComplexBands_local, last_cplx);
        auto* restrict offload_scratch_iat_ptr = offload_scratch_ptr + padded_size * iat;
        auto* restrict psi_iat_ptr             = results_scratch_ptr + orb_size * iat;
        auto* restrict pos_scratch             = psiinv_ptr + orb_size;

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr,
                                              ST(pos_scratch[iat * 6 + 3]),
                                              ST(pos_scratch[iat * 6 + 4]),
                                              ST(pos_scratch[iat * 6 + 5]),
                                              ix, iy, iz, a, b, c);

        TT sum(0);
        PRAGMA_OFFLOAD("omp parallel")
        {
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, a, b, c,
                                             offload_scratch_iat_ptr + first, first,
                                             last);
          C2R::assign_v(ST(pos_scratch[iat * 6]), ST(pos_scratch[iat * 6 + 1]), ST(pos_scratch[iat * 6 + 2]),
                        psi_iat_ptr, orb_size, offload_scratch_iat_ptr, myKcart_ptr, myKcart_padded_size,
                        first_spo_local, nComplexBands_local, first / 2, last / 2);

          PRAGMA_OFFLOAD("omp for reduction(+:sum)")
          for (int i = first_real; i < last_real; i++)
            sum += psi_iat_ptr[i] * psiinv_ptr[i];
        }
        ratios_private_ptr[iat * NumTeams + team_id] = sum;
      }

    // do the reduction manually
    for (int iat = 0; iat < nVP; ++iat)
    {
      ratios[iat] = TT(0);
      for (int tid = 0; tid < NumTeams; tid++)
        ratios[iat] += ratios_private[iat][tid];
    }
  }

  /** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
  template<typename VV, typename GV>
  inline void assign_vgl_from_l(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
  {
    constexpr ST two(2);
    const ST x = r[0], y = r[1], z = r[2];

    const ST* restrict k0 = myKcart.data(0);
    ASSUME_ALIGNED(k0);
    const ST* restrict k1 = myKcart.data(1);
    ASSUME_ALIGNED(k1);
    const ST* restrict k2 = myKcart.data(2);
    ASSUME_ALIGNED(k2);

    const ST* restrict g0 = myG.data(0);
    ASSUME_ALIGNED(g0);
    const ST* restrict g1 = myG.data(1);
    ASSUME_ALIGNED(g1);
    const ST* restrict g2 = myG.data(2);
    ASSUME_ALIGNED(g2);

    const size_t N = kPoints.size();

#pragma omp simd
    for (size_t j = 0; j < nComplexBands; j++)
    {
      const size_t jr = j << 1;
      const size_t ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g0[jr];
      const ST dY_r = g1[jr];
      const ST dZ_r = g2[jr];

      const ST dX_i = g0[ji];
      const ST dY_i = g1[ji];
      const ST dZ_i = g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r = dX_r + val_i * kX;
      const ST gY_r = dY_r + val_i * kY;
      const ST gZ_r = dZ_r + val_i * kZ;
      const ST gX_i = dX_i - val_r * kX;
      const ST gY_i = dY_i - val_r * kY;
      const ST gZ_i = dZ_i - val_r * kZ;

      const ST lap_r = myL[jr] + mKK[j] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
      const ST lap_i = myL[ji] + mKK[j] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);

      //this will be fixed later
      const size_t psiIndex = first_spo + jr;
      psi[psiIndex]         = c * val_r - s * val_i;
      psi[psiIndex + 1]     = c * val_i + s * val_r;
      d2psi[psiIndex]       = c * lap_r - s * lap_i;
      d2psi[psiIndex + 1]   = c * lap_i + s * lap_r;
      //this will go way with Determinant
      dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
      dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
      dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;
      dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
      dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
      dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;
    }

#pragma omp simd
    for (size_t j = nComplexBands; j < N; j++)
    {
      const size_t jr = j << 1;
      const size_t ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g0[jr];
      const ST dY_r = g1[jr];
      const ST dZ_r = g2[jr];

      const ST dX_i = g0[ji];
      const ST dY_i = g1[ji];
      const ST dZ_i = g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r         = dX_r + val_i * kX;
      const ST gY_r         = dY_r + val_i * kY;
      const ST gZ_r         = dZ_r + val_i * kZ;
      const ST gX_i         = dX_i - val_r * kX;
      const ST gY_i         = dY_i - val_r * kY;
      const ST gZ_i         = dZ_i - val_r * kZ;
      const size_t psiIndex = first_spo + nComplexBands + j;
      psi[psiIndex]         = c * val_r - s * val_i;
      //this will be fixed later
      dpsi[psiIndex][0] = c * gX_r - s * gX_i;
      dpsi[psiIndex][1] = c * gY_r - s * gY_i;
      dpsi[psiIndex][2] = c * gZ_r - s * gZ_i;

      const ST lap_r  = myL[jr] + mKK[j] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
      const ST lap_i  = myL[ji] + mKK[j] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);
      d2psi[psiIndex] = c * lap_r - s * lap_i;
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r = P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));

    const int ChunkSizePerTeam = 128;
    const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    const auto padded_size = myV.size();
    if (offload_scratch.size() < padded_size * 10)
      offload_scratch.resize(padded_size * 10);
    const auto orb_size = psi.size();
    if (results_scratch.size() < orb_size * 5)
      results_scratch.resize(orb_size * 5);

    // Ye: need to extract sizes and pointers before entering target region
    const auto* spline_ptr    = SplineInst->getSplinePtr();
    auto* offload_scratch_ptr = offload_scratch.data();
    auto* results_scratch_ptr = results_scratch.data();
    const auto x = r[0], y = r[1], z = r[2];
    const auto rux = ru[0], ruy = ru[1], ruz = ru[2];
    const auto myKcart_padded_size = myKcart.capacity();
    auto* mKK_ptr                  = master_mKK_ptr;
    auto* GGt_ptr                  = master_GGt_ptr;
    auto* PrimLattice_G_ptr        = master_PrimLattice_G_ptr;
    auto* myKcart_ptr              = master_myKcart_ptr;
    const size_t first_spo_local   = first_spo;
    const int nComplexBands_local  = nComplexBands;

    PRAGMA_OFFLOAD("omp target teams distribute num_teams(NumTeams) thread_limit(ChunkSizePerTeam) \
                map(always, from: results_scratch_ptr[0:orb_size*5])")
    for (int team_id = 0; team_id < NumTeams; team_id++)
    {
      const int first = ChunkSizePerTeam * team_id;
      const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;

      int ix, iy, iz;
      ST a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
      spline2::computeLocationAndFractional(spline_ptr, rux, ruy, ruz, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

      const ST G[9] = {PrimLattice_G_ptr[0], PrimLattice_G_ptr[1], PrimLattice_G_ptr[2],
                       PrimLattice_G_ptr[3], PrimLattice_G_ptr[4], PrimLattice_G_ptr[5],
                       PrimLattice_G_ptr[6], PrimLattice_G_ptr[7], PrimLattice_G_ptr[8]};
      const ST symGGt[6] = {GGt_ptr[0], GGt_ptr[1] + GGt_ptr[3], GGt_ptr[2] + GGt_ptr[6],
                            GGt_ptr[4], GGt_ptr[5] + GGt_ptr[7], GGt_ptr[8]};

      PRAGMA_OFFLOAD("omp parallel")
      {
        spline2offload::evaluate_vgh_impl_v2(spline_ptr,
                                             ix, iy, iz,
                                             a, b, c,
                                             da, db, dc,
                                             d2a, d2b, d2c,
                                             offload_scratch_ptr + first,
                                             offload_scratch_ptr + padded_size + first,
                                             offload_scratch_ptr + padded_size * 4 + first, padded_size, first, last);
        C2R::assign_vgl(x, y, z, results_scratch_ptr, mKK_ptr, orb_size, offload_scratch_ptr, padded_size, symGGt,
                        G, myKcart_ptr, myKcart_padded_size, first_spo_local, nComplexBands_local,
                        first / 2, last / 2);
      }
    }

    for (size_t i = 0; i < orb_size; i++)
    {
      psi[i]     = results_scratch[i];
      dpsi[i][0] = results_scratch[orb_size + i * 3];
      dpsi[i][1] = results_scratch[orb_size + i * 3 + 1];
      dpsi[i][2] = results_scratch[orb_size + i * 3 + 2];
      d2psi[i]   = results_scratch[orb_size * 4 + i];
    }
  }

  template<typename VV, typename GV, typename GGV>
  void assign_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi, int first, int last) const
  {
    // protect last
    last = last > kPoints.size() ? kPoints.size() : last;

    const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
             g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
             g22 = PrimLattice.G(8);
    const ST x = r[0], y = r[1], z = r[2];

    const ST* restrict k0 = myKcart.data(0);
    const ST* restrict k1 = myKcart.data(1);
    const ST* restrict k2 = myKcart.data(2);

    const ST* restrict g0  = myG.data(0);
    const ST* restrict g1  = myG.data(1);
    const ST* restrict g2  = myG.data(2);
    const ST* restrict h00 = myH.data(0);
    const ST* restrict h01 = myH.data(1);
    const ST* restrict h02 = myH.data(2);
    const ST* restrict h11 = myH.data(3);
    const ST* restrict h12 = myH.data(4);
    const ST* restrict h22 = myH.data(5);

#pragma omp simd
    for (size_t j = first; j < std::min(nComplexBands, last); j++)
    {
      int jr = j << 1;
      int ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
      const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
      const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

      const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
      const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
      const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r = dX_r + val_i * kX;
      const ST gY_r = dY_r + val_i * kY;
      const ST gZ_r = dZ_r + val_i * kZ;
      const ST gX_i = dX_i - val_r * kX;
      const ST gY_i = dY_i - val_r * kY;
      const ST gZ_i = dZ_i - val_r * kZ;

      const size_t psiIndex = first_spo + jr;

      psi[psiIndex]     = c * val_r - s * val_i;
      dpsi[psiIndex][0] = c * gX_r - s * gX_i;
      dpsi[psiIndex][1] = c * gY_r - s * gY_i;
      dpsi[psiIndex][2] = c * gZ_r - s * gZ_i;

      psi[psiIndex + 1]     = c * val_i + s * val_r;
      dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
      dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
      dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;

      const ST h_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02) +
          kX * (gX_i + dX_i);
      const ST h_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12) +
          kX * (gY_i + dY_i);
      const ST h_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22) +
          kX * (gZ_i + dZ_i);
      const ST h_yx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g00, g01, g02) +
          kY * (gX_i + dX_i);
      const ST h_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12) +
          kY * (gY_i + dY_i);
      const ST h_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22) +
          kY * (gZ_i + dZ_i);
      const ST h_zx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g00, g01, g02) +
          kZ * (gX_i + dX_i);
      const ST h_zy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g10, g11, g12) +
          kZ * (gY_i + dY_i);
      const ST h_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22) +
          kZ * (gZ_i + dZ_i);

      const ST h_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02) -
          kX * (gX_r + dX_r);
      const ST h_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12) -
          kX * (gY_r + dY_r);
      const ST h_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22) -
          kX * (gZ_r + dZ_r);
      const ST h_yx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g00, g01, g02) -
          kY * (gX_r + dX_r);
      const ST h_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12) -
          kY * (gY_r + dY_r);
      const ST h_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22) -
          kY * (gZ_r + dZ_r);
      const ST h_zx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g00, g01, g02) -
          kZ * (gX_r + dX_r);
      const ST h_zy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g10, g11, g12) -
          kZ * (gY_r + dY_r);
      const ST h_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22) -
          kZ * (gZ_r + dZ_r);

      grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
      grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][3] = c * h_yx_r - s * h_yx_i;
      grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
      grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][6] = c * h_zx_r - s * h_zx_i;
      grad_grad_psi[psiIndex][7] = c * h_zy_r - s * h_zy_i;
      grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

      grad_grad_psi[psiIndex + 1][0] = c * h_xx_i + s * h_xx_r;
      grad_grad_psi[psiIndex + 1][1] = c * h_xy_i + s * h_xy_r;
      grad_grad_psi[psiIndex + 1][2] = c * h_xz_i + s * h_xz_r;
      grad_grad_psi[psiIndex + 1][3] = c * h_yx_i + s * h_yx_r;
      grad_grad_psi[psiIndex + 1][4] = c * h_yy_i + s * h_yy_r;
      grad_grad_psi[psiIndex + 1][5] = c * h_yz_i + s * h_yz_r;
      grad_grad_psi[psiIndex + 1][6] = c * h_zx_i + s * h_zx_r;
      grad_grad_psi[psiIndex + 1][7] = c * h_zy_i + s * h_zy_r;
      grad_grad_psi[psiIndex + 1][8] = c * h_zz_i + s * h_zz_r;
    }

#pragma omp simd
    for (size_t j = std::max(nComplexBands, first); j < last; j++)
    {
      int jr = j << 1;
      int ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
      const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
      const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

      const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
      const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
      const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r = dX_r + val_i * kX;
      const ST gY_r = dY_r + val_i * kY;
      const ST gZ_r = dZ_r + val_i * kZ;
      const ST gX_i = dX_i - val_r * kX;
      const ST gY_i = dY_i - val_r * kY;
      const ST gZ_i = dZ_i - val_r * kZ;

      const size_t psiIndex = first_spo + nComplexBands + j;

      psi[psiIndex]     = c * val_r - s * val_i;
      dpsi[psiIndex][0] = c * gX_r - s * gX_i;
      dpsi[psiIndex][1] = c * gY_r - s * gY_i;
      dpsi[psiIndex][2] = c * gZ_r - s * gZ_i;

      const ST h_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02) +
          kX * (gX_i + dX_i);
      const ST h_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12) +
          kX * (gY_i + dY_i);
      const ST h_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22) +
          kX * (gZ_i + dZ_i);
      const ST h_yx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g00, g01, g02) +
          kY * (gX_i + dX_i);
      const ST h_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12) +
          kY * (gY_i + dY_i);
      const ST h_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22) +
          kY * (gZ_i + dZ_i);
      const ST h_zx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g00, g01, g02) +
          kZ * (gX_i + dX_i);
      const ST h_zy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g10, g11, g12) +
          kZ * (gY_i + dY_i);
      const ST h_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22) +
          kZ * (gZ_i + dZ_i);

      const ST h_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02) -
          kX * (gX_r + dX_r);
      const ST h_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12) -
          kX * (gY_r + dY_r);
      const ST h_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22) -
          kX * (gZ_r + dZ_r);
      const ST h_yx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g00, g01, g02) -
          kY * (gX_r + dX_r);
      const ST h_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12) -
          kY * (gY_r + dY_r);
      const ST h_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22) -
          kY * (gZ_r + dZ_r);
      const ST h_zx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g00, g01, g02) -
          kZ * (gX_r + dX_r);
      const ST h_zy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g10, g11, g12) -
          kZ * (gY_r + dY_r);
      const ST h_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22) -
          kZ * (gZ_r + dZ_r);

      grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
      grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][3] = c * h_yx_r - s * h_yx_i;
      grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
      grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][6] = c * h_zx_r - s * h_zx_i;
      grad_grad_psi[psiIndex][7] = c * h_zy_r - s * h_zy_i;
      grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;
    }
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r = P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
#pragma omp parallel
    {
      int first, last;
      FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
      assign_vgh(r, psi, dpsi, grad_grad_psi, first / 2, last / 2);
    }
  }

  template<typename VV, typename GV, typename GGV, typename GGGV>
  void assign_vghgh(const PointType& r,
                    VV& psi,
                    GV& dpsi,
                    GGV& grad_grad_psi,
                    GGGV& grad_grad_grad_psi,
                    int first = 0,
                    int last  = -1) const
  {
    // protect last
    last = last < 0 ? kPoints.size() : (last > kPoints.size() ? kPoints.size() : last);

    const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
             g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
             g22 = PrimLattice.G(8);
    const ST x = r[0], y = r[1], z = r[2];

    const ST* restrict k0 = myKcart.data(0);
    const ST* restrict k1 = myKcart.data(1);
    const ST* restrict k2 = myKcart.data(2);

    const ST* restrict g0  = myG.data(0);
    const ST* restrict g1  = myG.data(1);
    const ST* restrict g2  = myG.data(2);
    const ST* restrict h00 = myH.data(0);
    const ST* restrict h01 = myH.data(1);
    const ST* restrict h02 = myH.data(2);
    const ST* restrict h11 = myH.data(3);
    const ST* restrict h12 = myH.data(4);
    const ST* restrict h22 = myH.data(5);

    const ST* restrict gh000 = mygH.data(0);
    const ST* restrict gh001 = mygH.data(1);
    const ST* restrict gh002 = mygH.data(2);
    const ST* restrict gh011 = mygH.data(3);
    const ST* restrict gh012 = mygH.data(4);
    const ST* restrict gh022 = mygH.data(5);
    const ST* restrict gh111 = mygH.data(6);
    const ST* restrict gh112 = mygH.data(7);
    const ST* restrict gh122 = mygH.data(8);
    const ST* restrict gh222 = mygH.data(9);

//SIMD doesn't work quite right yet.  Comment out until further debugging.
#pragma omp simd
    for (size_t j = first; j < std::min(nComplexBands, last); j++)
    {
      int jr = j << 1;
      int ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
      const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
      const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

      const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
      const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
      const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r = dX_r + val_i * kX;
      const ST gY_r = dY_r + val_i * kY;
      const ST gZ_r = dZ_r + val_i * kZ;
      const ST gX_i = dX_i - val_r * kX;
      const ST gY_i = dY_i - val_r * kY;
      const ST gZ_i = dZ_i - val_r * kZ;

      const size_t psiIndex = first_spo + jr;
      psi[psiIndex]         = c * val_r - s * val_i;
      dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
      dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
      dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;

      psi[psiIndex + 1]     = c * val_i + s * val_r;
      dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
      dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
      dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;

      //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
      const ST f_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02);
      const ST f_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12);
      const ST f_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22);
      const ST f_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12);
      const ST f_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22);
      const ST f_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22);

      const ST f_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02);
      const ST f_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12);
      const ST f_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22);
      const ST f_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12);
      const ST f_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22);
      const ST f_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22);

      const ST h_xx_r = f_xx_r + 2 * kX * dX_i - kX * kX * val_r;
      const ST h_xy_r = f_xy_r + (kX * dY_i + kY * dX_i) - kX * kY * val_r;
      const ST h_xz_r = f_xz_r + (kX * dZ_i + kZ * dX_i) - kX * kZ * val_r;
      const ST h_yy_r = f_yy_r + 2 * kY * dY_i - kY * kY * val_r;
      const ST h_yz_r = f_yz_r + (kY * dZ_i + kZ * dY_i) - kY * kZ * val_r;
      const ST h_zz_r = f_zz_r + 2 * kZ * dZ_i - kZ * kZ * val_r;

      const ST h_xx_i = f_xx_i - 2 * kX * dX_r - kX * kX * val_i;
      const ST h_xy_i = f_xy_i - (kX * dY_r + kY * dX_r) - kX * kY * val_i;
      const ST h_xz_i = f_xz_i - (kX * dZ_r + kZ * dX_r) - kX * kZ * val_i;
      const ST h_yy_i = f_yy_i - 2 * kY * dY_r - kY * kY * val_i;
      const ST h_yz_i = f_yz_i - (kZ * dY_r + kY * dZ_r) - kZ * kY * val_i;
      const ST h_zz_i = f_zz_i - 2 * kZ * dZ_r - kZ * kZ * val_i;

      grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
      grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][3] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
      grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][6] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][7] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

      grad_grad_psi[psiIndex + 1][0] = c * h_xx_i + s * h_xx_r;
      grad_grad_psi[psiIndex + 1][1] = c * h_xy_i + s * h_xy_r;
      grad_grad_psi[psiIndex + 1][2] = c * h_xz_i + s * h_xz_r;
      grad_grad_psi[psiIndex + 1][3] = c * h_xy_i + s * h_xy_r;
      grad_grad_psi[psiIndex + 1][4] = c * h_yy_i + s * h_yy_r;
      grad_grad_psi[psiIndex + 1][5] = c * h_yz_i + s * h_yz_r;
      grad_grad_psi[psiIndex + 1][6] = c * h_xz_i + s * h_xz_r;
      grad_grad_psi[psiIndex + 1][7] = c * h_yz_i + s * h_yz_r;
      grad_grad_psi[psiIndex + 1][8] = c * h_zz_i + s * h_zz_r;

      //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
      // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

      const ST f3_xxx_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g00, g01, g02);
      const ST f3_xxy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g10, g11, g12);
      const ST f3_xxz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g20, g21, g22);
      const ST f3_xyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g10, g11, g12);
      const ST f3_xyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g20, g21, g22);
      const ST f3_xzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g20, g21, g22, g20, g21, g22);
      const ST f3_yyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g10, g11, g12);
      const ST f3_yyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g20, g21, g22);
      const ST f3_yzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g20, g21, g22, g20, g21, g22);
      const ST f3_zzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g20, g21, g22, g20, g21, g22, g20, g21, g22);

      const ST f3_xxx_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g00, g01, g02);
      const ST f3_xxy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g10, g11, g12);
      const ST f3_xxz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g20, g21, g22);
      const ST f3_xyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g10, g11, g12);
      const ST f3_xyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g20, g21, g22);
      const ST f3_xzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g20, g21, g22, g20, g21, g22);
      const ST f3_yyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g10, g11, g12);
      const ST f3_yyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g20, g21, g22);
      const ST f3_yzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g20, g21, g22, g20, g21, g22);
      const ST f3_zzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g20, g21, g22, g20, g21, g22, g20, g21, g22);

      //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
      const ST gh_xxx_r = f3_xxx_r + 3 * kX * f_xx_i - 3 * kX * kX * dX_r - kX * kX * kX * val_i;
      const ST gh_xxx_i = f3_xxx_i - 3 * kX * f_xx_r - 3 * kX * kX * dX_i + kX * kX * kX * val_r;
      const ST gh_xxy_r =
          f3_xxy_r + (kY * f_xx_i + 2 * kX * f_xy_i) - (kX * kX * dY_r + 2 * kX * kY * dX_r) - kX * kX * kY * val_i;
      const ST gh_xxy_i =
          f3_xxy_i - (kY * f_xx_r + 2 * kX * f_xy_r) - (kX * kX * dY_i + 2 * kX * kY * dX_i) + kX * kX * kY * val_r;
      const ST gh_xxz_r =
          f3_xxz_r + (kZ * f_xx_i + 2 * kX * f_xz_i) - (kX * kX * dZ_r + 2 * kX * kZ * dX_r) - kX * kX * kZ * val_i;
      const ST gh_xxz_i =
          f3_xxz_i - (kZ * f_xx_r + 2 * kX * f_xz_r) - (kX * kX * dZ_i + 2 * kX * kZ * dX_i) + kX * kX * kZ * val_r;
      const ST gh_xyy_r =
          f3_xyy_r + (2 * kY * f_xy_i + kX * f_yy_i) - (2 * kX * kY * dY_r + kY * kY * dX_r) - kX * kY * kY * val_i;
      const ST gh_xyy_i =
          f3_xyy_i - (2 * kY * f_xy_r + kX * f_yy_r) - (2 * kX * kY * dY_i + kY * kY * dX_i) + kX * kY * kY * val_r;
      const ST gh_xyz_r = f3_xyz_r + (kX * f_yz_i + kY * f_xz_i + kZ * f_xy_i) -
          (kX * kY * dZ_r + kY * kZ * dX_r + kZ * kX * dY_r) - kX * kY * kZ * val_i;
      const ST gh_xyz_i = f3_xyz_i - (kX * f_yz_r + kY * f_xz_r + kZ * f_xy_r) -
          (kX * kY * dZ_i + kY * kZ * dX_i + kZ * kX * dY_i) + kX * kY * kZ * val_r;
      const ST gh_xzz_r =
          f3_xzz_r + (2 * kZ * f_xz_i + kX * f_zz_i) - (2 * kX * kZ * dZ_r + kZ * kZ * dX_r) - kX * kZ * kZ * val_i;
      const ST gh_xzz_i =
          f3_xzz_i - (2 * kZ * f_xz_r + kX * f_zz_r) - (2 * kX * kZ * dZ_i + kZ * kZ * dX_i) + kX * kZ * kZ * val_r;
      const ST gh_yyy_r = f3_yyy_r + 3 * kY * f_yy_i - 3 * kY * kY * dY_r - kY * kY * kY * val_i;
      const ST gh_yyy_i = f3_yyy_i - 3 * kY * f_yy_r - 3 * kY * kY * dY_i + kY * kY * kY * val_r;
      const ST gh_yyz_r =
          f3_yyz_r + (kZ * f_yy_i + 2 * kY * f_yz_i) - (kY * kY * dZ_r + 2 * kY * kZ * dY_r) - kY * kY * kZ * val_i;
      const ST gh_yyz_i =
          f3_yyz_i - (kZ * f_yy_r + 2 * kY * f_yz_r) - (kY * kY * dZ_i + 2 * kY * kZ * dY_i) + kY * kY * kZ * val_r;
      const ST gh_yzz_r =
          f3_yzz_r + (2 * kZ * f_yz_i + kY * f_zz_i) - (2 * kY * kZ * dZ_r + kZ * kZ * dY_r) - kY * kZ * kZ * val_i;
      const ST gh_yzz_i =
          f3_yzz_i - (2 * kZ * f_yz_r + kY * f_zz_r) - (2 * kY * kZ * dZ_i + kZ * kZ * dY_i) + kY * kZ * kZ * val_r;
      const ST gh_zzz_r = f3_zzz_r + 3 * kZ * f_zz_i - 3 * kZ * kZ * dZ_r - kZ * kZ * kZ * val_i;
      const ST gh_zzz_i = f3_zzz_i - 3 * kZ * f_zz_r - 3 * kZ * kZ * dZ_i + kZ * kZ * kZ * val_r;

      grad_grad_grad_psi[psiIndex][0][0] = c * gh_xxx_r - s * gh_xxx_i;
      grad_grad_grad_psi[psiIndex][0][1] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][0][2] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][0][3] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][0][4] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][0][5] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][0][6] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][0][7] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][0][8] = c * gh_xzz_r - s * gh_xzz_i;

      grad_grad_grad_psi[psiIndex][1][0] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][1][1] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][1][2] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][1][3] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][1][4] = c * gh_yyy_r - s * gh_yyy_i;
      grad_grad_grad_psi[psiIndex][1][5] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][1][6] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][1][7] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][1][8] = c * gh_yzz_r - s * gh_yzz_i;

      grad_grad_grad_psi[psiIndex][2][0] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][2][1] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][2][2] = c * gh_xzz_r - s * gh_xzz_i;
      grad_grad_grad_psi[psiIndex][2][3] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][2][4] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][2][5] = c * gh_yzz_r - s * gh_yzz_i;
      grad_grad_grad_psi[psiIndex][2][6] = c * gh_xzz_r - s * gh_xzz_i;
      grad_grad_grad_psi[psiIndex][2][7] = c * gh_yzz_r - s * gh_yzz_i;
      grad_grad_grad_psi[psiIndex][2][8] = c * gh_zzz_r - s * gh_zzz_i;

      grad_grad_grad_psi[psiIndex + 1][0][0] = c * gh_xxx_i + s * gh_xxx_r;
      grad_grad_grad_psi[psiIndex + 1][0][1] = c * gh_xxy_i + s * gh_xxy_r;
      grad_grad_grad_psi[psiIndex + 1][0][2] = c * gh_xxz_i + s * gh_xxz_r;
      grad_grad_grad_psi[psiIndex + 1][0][3] = c * gh_xxy_i + s * gh_xxy_r;
      grad_grad_grad_psi[psiIndex + 1][0][4] = c * gh_xyy_i + s * gh_xyy_r;
      grad_grad_grad_psi[psiIndex + 1][0][5] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][0][6] = c * gh_xxz_i + s * gh_xxz_r;
      grad_grad_grad_psi[psiIndex + 1][0][7] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][0][8] = c * gh_xzz_i + s * gh_xzz_r;

      grad_grad_grad_psi[psiIndex + 1][1][0] = c * gh_xxy_i + s * gh_xxy_r;
      grad_grad_grad_psi[psiIndex + 1][1][1] = c * gh_xyy_i + s * gh_xyy_r;
      grad_grad_grad_psi[psiIndex + 1][1][2] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][1][3] = c * gh_xyy_i + s * gh_xyy_r;
      grad_grad_grad_psi[psiIndex + 1][1][4] = c * gh_yyy_i + s * gh_yyy_r;
      grad_grad_grad_psi[psiIndex + 1][1][5] = c * gh_yyz_i + s * gh_yyz_r;
      grad_grad_grad_psi[psiIndex + 1][1][6] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][1][7] = c * gh_yyz_i + s * gh_yyz_r;
      grad_grad_grad_psi[psiIndex + 1][1][8] = c * gh_yzz_i + s * gh_yzz_r;

      grad_grad_grad_psi[psiIndex + 1][2][0] = c * gh_xxz_i + s * gh_xxz_r;
      grad_grad_grad_psi[psiIndex + 1][2][1] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][2][2] = c * gh_xzz_i + s * gh_xzz_r;
      grad_grad_grad_psi[psiIndex + 1][2][3] = c * gh_xyz_i + s * gh_xyz_r;
      grad_grad_grad_psi[psiIndex + 1][2][4] = c * gh_yyz_i + s * gh_yyz_r;
      grad_grad_grad_psi[psiIndex + 1][2][5] = c * gh_yzz_i + s * gh_yzz_r;
      grad_grad_grad_psi[psiIndex + 1][2][6] = c * gh_xzz_i + s * gh_xzz_r;
      grad_grad_grad_psi[psiIndex + 1][2][7] = c * gh_yzz_i + s * gh_yzz_r;
      grad_grad_grad_psi[psiIndex + 1][2][8] = c * gh_zzz_i + s * gh_zzz_r;
    }
#pragma omp simd
    for (size_t j = std::max(nComplexBands, first); j < last; j++)
    {
      int jr = j << 1;
      int ji = jr + 1;

      const ST kX    = k0[j];
      const ST kY    = k1[j];
      const ST kZ    = k2[j];
      const ST val_r = myV[jr];
      const ST val_i = myV[ji];

      //phase
      ST s, c;
      sincos(-(x * kX + y * kY + z * kZ), &s, &c);

      //dot(PrimLattice.G,myG[j])
      const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
      const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
      const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

      const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
      const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
      const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

      // \f$\nabla \psi_r + {\bf k}\psi_i\f$
      const ST gX_r = dX_r + val_i * kX;
      const ST gY_r = dY_r + val_i * kY;
      const ST gZ_r = dZ_r + val_i * kZ;
      const ST gX_i = dX_i - val_r * kX;
      const ST gY_i = dY_i - val_r * kY;
      const ST gZ_i = dZ_i - val_r * kZ;

      const size_t psiIndex = first_spo + nComplexBands + j;
      psi[psiIndex]         = c * val_r - s * val_i;
      dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
      dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
      dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;

      //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
      const ST f_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02);
      const ST f_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12);
      const ST f_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22);
      const ST f_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12);
      const ST f_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22);
      const ST f_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22);

      const ST f_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02);
      const ST f_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12);
      const ST f_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22);
      const ST f_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12);
      const ST f_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22);
      const ST f_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22);

      const ST h_xx_r = f_xx_r + 2 * kX * dX_i - kX * kX * val_r;
      const ST h_xy_r = f_xy_r + (kX * dY_i + kY * dX_i) - kX * kY * val_r;
      const ST h_xz_r = f_xz_r + (kX * dZ_i + kZ * dX_i) - kX * kZ * val_r;
      const ST h_yy_r = f_yy_r + 2 * kY * dY_i - kY * kY * val_r;
      const ST h_yz_r = f_yz_r + (kY * dZ_i + kZ * dY_i) - kY * kZ * val_r;
      const ST h_zz_r = f_zz_r + 2 * kZ * dZ_i - kZ * kZ * val_r;

      const ST h_xx_i = f_xx_i - 2 * kX * dX_r - kX * kX * val_i;
      const ST h_xy_i = f_xy_i - (kX * dY_r + kY * dX_r) - kX * kY * val_i;
      const ST h_xz_i = f_xz_i - (kX * dZ_r + kZ * dX_r) - kX * kZ * val_i;
      const ST h_yy_i = f_yy_i - 2 * kY * dY_r - kY * kY * val_i;
      const ST h_yz_i = f_yz_i - (kZ * dY_r + kY * dZ_r) - kZ * kY * val_i;
      const ST h_zz_i = f_zz_i - 2 * kZ * dZ_r - kZ * kZ * val_i;

      grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
      grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][3] = c * h_xy_r - s * h_xy_i;
      grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
      grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][6] = c * h_xz_r - s * h_xz_i;
      grad_grad_psi[psiIndex][7] = c * h_yz_r - s * h_yz_i;
      grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

      //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
      // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

      const ST f3_xxx_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g00, g01, g02);
      const ST f3_xxy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g10, g11, g12);
      const ST f3_xxz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g20, g21, g22);
      const ST f3_xyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g10, g11, g12);
      const ST f3_xyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g20, g21, g22);
      const ST f3_xzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g20, g21, g22, g20, g21, g22);
      const ST f3_yyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g10, g11, g12);
      const ST f3_yyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g20, g21, g22);
      const ST f3_yzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g20, g21, g22, g20, g21, g22);
      const ST f3_zzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                      gh112[jr], gh122[jr], gh222[jr], g20, g21, g22, g20, g21, g22, g20, g21, g22);

      const ST f3_xxx_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g00, g01, g02);
      const ST f3_xxy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g10, g11, g12);
      const ST f3_xxz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g20, g21, g22);
      const ST f3_xyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g10, g11, g12);
      const ST f3_xyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g20, g21, g22);
      const ST f3_xzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g20, g21, g22, g20, g21, g22);
      const ST f3_yyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g10, g11, g12);
      const ST f3_yyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g20, g21, g22);
      const ST f3_yzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g20, g21, g22, g20, g21, g22);
      const ST f3_zzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                      gh112[ji], gh122[ji], gh222[ji], g20, g21, g22, g20, g21, g22, g20, g21, g22);

      //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
      const ST gh_xxx_r = f3_xxx_r + 3 * kX * f_xx_i - 3 * kX * kX * dX_r - kX * kX * kX * val_i;
      const ST gh_xxx_i = f3_xxx_i - 3 * kX * f_xx_r - 3 * kX * kX * dX_i + kX * kX * kX * val_r;
      const ST gh_xxy_r =
          f3_xxy_r + (kY * f_xx_i + 2 * kX * f_xy_i) - (kX * kX * dY_r + 2 * kX * kY * dX_r) - kX * kX * kY * val_i;
      const ST gh_xxy_i =
          f3_xxy_i - (kY * f_xx_r + 2 * kX * f_xy_r) - (kX * kX * dY_i + 2 * kX * kY * dX_i) + kX * kX * kY * val_r;
      const ST gh_xxz_r =
          f3_xxz_r + (kZ * f_xx_i + 2 * kX * f_xz_i) - (kX * kX * dZ_r + 2 * kX * kZ * dX_r) - kX * kX * kZ * val_i;
      const ST gh_xxz_i =
          f3_xxz_i - (kZ * f_xx_r + 2 * kX * f_xz_r) - (kX * kX * dZ_i + 2 * kX * kZ * dX_i) + kX * kX * kZ * val_r;
      const ST gh_xyy_r =
          f3_xyy_r + (2 * kY * f_xy_i + kX * f_yy_i) - (2 * kX * kY * dY_r + kY * kY * dX_r) - kX * kY * kY * val_i;
      const ST gh_xyy_i =
          f3_xyy_i - (2 * kY * f_xy_r + kX * f_yy_r) - (2 * kX * kY * dY_i + kY * kY * dX_i) + kX * kY * kY * val_r;
      const ST gh_xyz_r = f3_xyz_r + (kX * f_yz_i + kY * f_xz_i + kZ * f_xy_i) -
          (kX * kY * dZ_r + kY * kZ * dX_r + kZ * kX * dY_r) - kX * kY * kZ * val_i;
      const ST gh_xyz_i = f3_xyz_i - (kX * f_yz_r + kY * f_xz_r + kZ * f_xy_r) -
          (kX * kY * dZ_i + kY * kZ * dX_i + kZ * kX * dY_i) + kX * kY * kZ * val_r;
      const ST gh_xzz_r =
          f3_xzz_r + (2 * kZ * f_xz_i + kX * f_zz_i) - (2 * kX * kZ * dZ_r + kZ * kZ * dX_r) - kX * kZ * kZ * val_i;
      const ST gh_xzz_i =
          f3_xzz_i - (2 * kZ * f_xz_r + kX * f_zz_r) - (2 * kX * kZ * dZ_i + kZ * kZ * dX_i) + kX * kZ * kZ * val_r;
      const ST gh_yyy_r = f3_yyy_r + 3 * kY * f_yy_i - 3 * kY * kY * dY_r - kY * kY * kY * val_i;
      const ST gh_yyy_i = f3_yyy_i - 3 * kY * f_yy_r - 3 * kY * kY * dY_i + kY * kY * kY * val_r;
      const ST gh_yyz_r =
          f3_yyz_r + (kZ * f_yy_i + 2 * kY * f_yz_i) - (kY * kY * dZ_r + 2 * kY * kZ * dY_r) - kY * kY * kZ * val_i;
      const ST gh_yyz_i =
          f3_yyz_i - (kZ * f_yy_r + 2 * kY * f_yz_r) - (kY * kY * dZ_i + 2 * kY * kZ * dY_i) + kY * kY * kZ * val_r;
      const ST gh_yzz_r =
          f3_yzz_r + (2 * kZ * f_yz_i + kY * f_zz_i) - (2 * kY * kZ * dZ_r + kZ * kZ * dY_r) - kY * kZ * kZ * val_i;
      const ST gh_yzz_i =
          f3_yzz_i - (2 * kZ * f_yz_r + kY * f_zz_r) - (2 * kY * kZ * dZ_i + kZ * kZ * dY_i) + kY * kZ * kZ * val_r;
      const ST gh_zzz_r = f3_zzz_r + 3 * kZ * f_zz_i - 3 * kZ * kZ * dZ_r - kZ * kZ * kZ * val_i;
      const ST gh_zzz_i = f3_zzz_i - 3 * kZ * f_zz_r - 3 * kZ * kZ * dZ_i + kZ * kZ * kZ * val_r;
      //[x][xx] //These are the unique entries
      grad_grad_grad_psi[psiIndex][0][0] = c * gh_xxx_r - s * gh_xxx_i;
      grad_grad_grad_psi[psiIndex][0][1] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][0][2] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][0][3] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][0][4] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][0][5] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][0][6] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][0][7] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][0][8] = c * gh_xzz_r - s * gh_xzz_i;

      grad_grad_grad_psi[psiIndex][1][0] = c * gh_xxy_r - s * gh_xxy_i;
      grad_grad_grad_psi[psiIndex][1][1] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][1][2] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][1][3] = c * gh_xyy_r - s * gh_xyy_i;
      grad_grad_grad_psi[psiIndex][1][4] = c * gh_yyy_r - s * gh_yyy_i;
      grad_grad_grad_psi[psiIndex][1][5] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][1][6] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][1][7] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][1][8] = c * gh_yzz_r - s * gh_yzz_i;

      grad_grad_grad_psi[psiIndex][2][0] = c * gh_xxz_r - s * gh_xxz_i;
      grad_grad_grad_psi[psiIndex][2][1] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][2][2] = c * gh_xzz_r - s * gh_xzz_i;
      grad_grad_grad_psi[psiIndex][2][3] = c * gh_xyz_r - s * gh_xyz_i;
      grad_grad_grad_psi[psiIndex][2][4] = c * gh_yyz_r - s * gh_yyz_i;
      grad_grad_grad_psi[psiIndex][2][5] = c * gh_yzz_r - s * gh_yzz_i;
      grad_grad_grad_psi[psiIndex][2][6] = c * gh_xzz_r - s * gh_xzz_i;
      grad_grad_grad_psi[psiIndex][2][7] = c * gh_yzz_r - s * gh_yzz_i;
      grad_grad_grad_psi[psiIndex][2][8] = c * gh_zzz_r - s * gh_zzz_i;
    }
  }

  template<typename VV, typename GV, typename GGV, typename GGGV>
  void evaluate_vghgh(const ParticleSet& P,
                      const int iat,
                      VV& psi,
                      GV& dpsi,
                      GGV& grad_grad_psi,
                      GGGV& grad_grad_grad_psi)
  {
    const PointType& r = P.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
#pragma omp parallel
    {
      int first, last;
      FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      spline2::evaluate3d_vghgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, mygH, first, last);
      assign_vghgh(r, psi, dpsi, grad_grad_psi, grad_grad_grad_psi, first / 2, last / 2);
    }
  }
};

} // namespace qmcplusplus
#endif
