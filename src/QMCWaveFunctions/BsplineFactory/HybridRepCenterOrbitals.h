//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

/** @file HybridRepCenterOrbitals.h
 *
 * Hybrid representation atomic centered orbitals
 */
#ifndef QMCPLUSPLUS_HYBRIDREP_CENTER_ORBITALS_H
#define QMCPLUSPLUS_HYBRIDREP_CENTER_ORBITALS_H

#include "Particle/DistanceTable.h"
#include "Particle/VirtualParticleSet.h"
#include "Numerics/SoaSphericalTensor.h"
#include "spline2/MultiBspline1D.hpp"
#include "Numerics/SmoothFunctions.hpp"
#include "hdf/hdf_archive.h"

namespace qmcplusplus
{
template<class BSPLINESPO>
class HybridRepSetReader;

template<typename ST>
class AtomicOrbitals
{
public:
  static const int D           = 3;
  using AtomicSplineType       = typename bspline_traits<ST, 1>::SplineType;
  using AtomicBCType           = typename bspline_traits<ST, 1>::BCType;
  using AtomicSingleSplineType = UBspline_1d_d;
  using PointType              = TinyVector<ST, D>;
  using value_type             = ST;

  using vContainer_type = aligned_vector<ST>;

private:
  // near core cutoff
  ST rmin;
  // far from core cutoff, rmin_sqrt>=rmin
  ST rmin_sqrt;
  ST cutoff, cutoff_buffer, spline_radius, non_overlapping_radius;
  int spline_npoints, BaseN;
  int NumBands, Npad;
  PointType center_pos;
  const int lmax, lm_tot;
  SoaSphericalTensor<ST> Ylm;
  vContainer_type l_vals;
  vContainer_type r_power_minus_l;
  ///1D spline of radial functions of all the orbitals
  std::shared_ptr<MultiBspline1D<ST>> SplineInst;
  ///coef copy for orbital rotation
  std::shared_ptr<std::vector<ST>> coef_copy_;

  vContainer_type localV, localG, localL;

public:
  AtomicOrbitals(int Lmax) : lmax(Lmax), lm_tot((Lmax + 1) * (Lmax + 1)), Ylm(Lmax)
  {
    r_power_minus_l.resize(lm_tot);
    l_vals.resize(lm_tot);
    for (int l = 0; l <= lmax; l++)
      for (int m = -l; m <= l; m++)
        l_vals[l * (l + 1) + m] = l;
    rmin      = std::exp(std::log(std::numeric_limits<ST>::min()) / std::max(Lmax, 1));
    rmin      = std::max(rmin, std::numeric_limits<ST>::epsilon());
    rmin_sqrt = std::max(rmin, std::sqrt(std::numeric_limits<ST>::epsilon()));
  }

  // accessing functions, const only
  ST getCutoff() const { return cutoff; }
  ST getCutoffBuffer() const { return cutoff_buffer; }
  ST getSplineRadius() const { return spline_radius; }
  ST getNonOverlappingRadius() const { return non_overlapping_radius; }
  int getSplineNpoints() const { return spline_npoints; }
  int getLmax() const { return lmax; }
  const PointType& getCenterPos() const { return center_pos; }

  inline void resizeStorage(size_t Nb)
  {
    NumBands = Nb;
    Npad     = getAlignedSize<ST>(Nb);
    localV.resize(Npad * lm_tot);
    localG.resize(Npad * lm_tot);
    localL.resize(Npad * lm_tot);
    create_spline();
  }

  void bcast_tables(Communicate* comm) { chunked_bcast(comm, SplineInst->getSplinePtr()); }

  void gather_tables(Communicate* comm, std::vector<int>& offset)
  {
    gatherv(comm, SplineInst->getSplinePtr(), Npad, offset);
  }

  template<typename PT, typename VT>
  inline void set_info(const PT& R,
                       const VT& cutoff_in,
                       const VT& cutoff_buffer_in,
                       const VT& spline_radius_in,
                       const VT& non_overlapping_radius_in,
                       const int spline_npoints_in)
  {
    center_pos[0]          = R[0];
    center_pos[1]          = R[1];
    center_pos[2]          = R[2];
    cutoff                 = cutoff_in;
    cutoff_buffer          = cutoff_buffer_in;
    spline_radius          = spline_radius_in;
    spline_npoints         = spline_npoints_in;
    non_overlapping_radius = non_overlapping_radius_in;
    BaseN                  = spline_npoints + 2;
  }

  inline void create_spline()
  {
    AtomicBCType bc;
    bc.lCode = FLAT;
    bc.rCode = NATURAL;
    Ugrid grid;
    grid.start = 0.0;
    grid.end   = spline_radius;
    grid.num   = spline_npoints;
    SplineInst = std::make_shared<MultiBspline1D<ST>>();
    SplineInst->create(grid, bc, lm_tot * Npad);
  }

  inline size_t getSplineSizeInBytes() const { return SplineInst->sizeInByte(); }

  inline void flush_zero() { SplineInst->flush_zero(); }

  inline void set_spline(AtomicSingleSplineType* spline, int lm, int ispline)
  {
    SplineInst->copy_spline(spline, lm * Npad + ispline, 0, BaseN);
  }

  bool read_splines(hdf_archive& h5f)
  {
    einspline_engine<AtomicSplineType> bigtable(SplineInst->getSplinePtr());
    int lmax_in = 0, spline_npoints_in = 0;
    ST spline_radius_in;
    if (!h5f.readEntry(lmax_in, "l_max") || lmax_in != lmax)
      return false;
    if (!h5f.readEntry(spline_radius_in, "spline_radius") || spline_radius_in != spline_radius)
      return false;
    if (!h5f.readEntry(spline_npoints_in, "spline_npoints") || spline_npoints_in != spline_npoints)
      return false;
    return h5f.readEntry(bigtable, "radial_spline");
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success = true;
    success      = success && h5f.writeEntry(spline_radius, "spline_radius");
    success      = success && h5f.writeEntry(spline_npoints, "spline_npoints");
    success      = success && h5f.writeEntry(lmax, "l_max");
    success      = success && h5f.writeEntry(center_pos, "position");
    einspline_engine<AtomicSplineType> bigtable(SplineInst->getSplinePtr());
    success = success && h5f.writeEntry(bigtable, "radial_spline");
    return success;
  }

  //evaluate only V
  template<typename VV>
  inline void evaluate_v(const ST& r, const PointType& dr, VV& myV)
  {
    if (r > std::numeric_limits<ST>::epsilon())
      Ylm.evaluateV(dr[0] / r, dr[1] / r, dr[2] / r);
    else
      Ylm.evaluateV(0, 0, 1);
    const ST* restrict Ylm_v = Ylm[0];

    constexpr ST czero(0);
    ST* restrict val       = myV.data();
    ST* restrict local_val = localV.data();
    std::fill(myV.begin(), myV.end(), czero);

    SplineInst->evaluate(r, localV);

    for (size_t lm = 0; lm < lm_tot; lm++)
    {
#pragma omp simd aligned(val, local_val : QMC_SIMD_ALIGNMENT)
      for (size_t ib = 0; ib < myV.size(); ib++)
        val[ib] += Ylm_v[lm] * local_val[ib];
      local_val += Npad;
    }
  }

  template<typename DISPL, typename VM>
  inline void evaluateValues(const DISPL& Displacements, const int center_idx, const ST& r, VM& multi_myV)
  {
    if (r <= std::numeric_limits<ST>::epsilon())
      Ylm.evaluateV(0, 0, 1);
    const ST* restrict Ylm_v = Ylm[0];

    const size_t m = multi_myV.cols();
    constexpr ST czero(0);
    std::fill(multi_myV.begin(), multi_myV.end(), czero);
    SplineInst->evaluate(r, localV);

    for (int ivp = 0; ivp < Displacements.size(); ivp++)
    {
      PointType dr = Displacements[ivp][center_idx];
      if (r > std::numeric_limits<ST>::epsilon())
        Ylm.evaluateV(-dr[0] / r, -dr[1] / r, -dr[2] / r);

      ST* restrict val       = multi_myV[ivp];
      ST* restrict local_val = localV.data();
      for (size_t lm = 0; lm < lm_tot; lm++)
      {
#pragma omp simd aligned(val, local_val : QMC_SIMD_ALIGNMENT)
        for (size_t ib = 0; ib < m; ib++)
          val[ib] += Ylm_v[lm] * local_val[ib];
        local_val += Npad;
      }
    }
  }

  //evaluate VGL
  template<typename VV, typename GV>
  inline void evaluate_vgl(const ST& r, const PointType& dr, VV& myV, GV& myG, VV& myL)
  {
    ST drx, dry, drz, rhatx, rhaty, rhatz, rinv;
    if (r > rmin)
    {
      rinv = 1.0 / r;
    }
    else
    {
      rinv = 0;
    }
    drx   = dr[0];
    dry   = dr[1];
    drz   = dr[2];
    rhatx = drx * rinv;
    rhaty = dry * rinv;
    rhatz = drz * rinv;

    Ylm.evaluateVGL(drx, dry, drz);
    const ST* restrict Ylm_v  = Ylm[0];
    const ST* restrict Ylm_gx = Ylm[1];
    const ST* restrict Ylm_gy = Ylm[2];
    const ST* restrict Ylm_gz = Ylm[3];

    ST* restrict g0 = myG.data(0);
    ST* restrict g1 = myG.data(1);
    ST* restrict g2 = myG.data(2);
    constexpr ST czero(0), cone(1), chalf(0.5);
    std::fill(myV.begin(), myV.end(), czero);
    std::fill(g0, g0 + Npad, czero);
    std::fill(g1, g1 + Npad, czero);
    std::fill(g2, g2 + Npad, czero);
    std::fill(myL.begin(), myL.end(), czero);
    ST* restrict val        = myV.data();
    ST* restrict lapl       = myL.data();
    ST* restrict local_val  = localV.data();
    ST* restrict local_grad = localG.data();
    ST* restrict local_lapl = localL.data();

    SplineInst->evaluate_vgl(r, localV, localG, localL);

    if (r > rmin_sqrt)
    {
      // far from core
      r_power_minus_l[0] = cone;
      ST r_power_temp    = cone;
      for (int l = 1; l <= lmax; l++)
      {
        r_power_temp *= rinv;
        for (int m = -l, lm = l * l; m <= l; m++, lm++)
          r_power_minus_l[lm] = r_power_temp;
      }

      for (size_t lm = 0; lm < lm_tot; lm++)
      {
        const ST& l_val      = l_vals[lm];
        const ST& r_power    = r_power_minus_l[lm];
        const ST Ylm_rescale = Ylm_v[lm] * r_power;
        const ST rhat_dot_G  = (rhatx * Ylm_gx[lm] + rhaty * Ylm_gy[lm] + rhatz * Ylm_gz[lm]) * r_power;
#pragma omp simd aligned(val, g0, g1, g2, lapl, local_val, local_grad, local_lapl : QMC_SIMD_ALIGNMENT)
        for (size_t ib = 0; ib < myV.size(); ib++)
        {
          const ST local_v = local_val[ib];
          const ST local_g = local_grad[ib];
          const ST local_l = local_lapl[ib];
          // value
          const ST Vpart = l_val * rinv * local_v;
          val[ib] += Ylm_rescale * local_v;

          // grad
          const ST factor1 = local_g * Ylm_rescale;
          const ST factor2 = local_v * r_power;
          const ST factor3 = -Vpart * Ylm_rescale;
          g0[ib] += factor1 * rhatx + factor2 * Ylm_gx[lm] + factor3 * rhatx;
          g1[ib] += factor1 * rhaty + factor2 * Ylm_gy[lm] + factor3 * rhaty;
          g2[ib] += factor1 * rhatz + factor2 * Ylm_gz[lm] + factor3 * rhatz;

          // laplacian
          lapl[ib] += (local_l + (local_g * (2 - l_val) - Vpart) * rinv) * Ylm_rescale + (local_g - Vpart) * rhat_dot_G;
        }
        local_val += Npad;
        local_grad += Npad;
        local_lapl += Npad;
      }
    }
    else if (r > rmin)
    {
      // the possibility of reaching here is very very low
      std::cout << "Warning: an electron is very close to an ion, distance=" << r << " be careful!" << std::endl;
      // near core, kill divergence in the laplacian
      r_power_minus_l[0] = cone;
      ST r_power_temp    = cone;
      for (int l = 1; l <= lmax; l++)
      {
        r_power_temp *= rinv;
        for (int m = -l, lm = l * l; m <= l; m++, lm++)
          r_power_minus_l[lm] = r_power_temp;
      }

      for (size_t lm = 0; lm < lm_tot; lm++)
      {
        const ST& l_val      = l_vals[lm];
        const ST& r_power    = r_power_minus_l[lm];
        const ST Ylm_rescale = Ylm_v[lm] * r_power;
        const ST rhat_dot_G  = (Ylm_gx[lm] * rhatx + Ylm_gy[lm] * rhaty + Ylm_gz[lm] * rhatz) * r_power * r;
#pragma omp simd aligned(val, g0, g1, g2, lapl, local_val, local_grad, local_lapl : QMC_SIMD_ALIGNMENT)
        for (size_t ib = 0; ib < myV.size(); ib++)
        {
          const ST local_v = local_val[ib];
          const ST local_g = local_grad[ib];
          const ST local_l = local_lapl[ib];
          // value
          const ST Vpart = Ylm_rescale * local_v;
          val[ib] += Vpart;

          // grad
          const ST factor1 = local_g * Ylm_rescale;
          const ST factor2 = local_v * r_power;
          const ST factor3 = -l_val * Vpart * rinv;
          g0[ib] += factor1 * rhatx + factor2 * Ylm_gx[lm] + factor3 * rhatx;
          g1[ib] += factor1 * rhaty + factor2 * Ylm_gy[lm] + factor3 * rhaty;
          g2[ib] += factor1 * rhatz + factor2 * Ylm_gz[lm] + factor3 * rhatz;

          // laplacian
          lapl[ib] += local_l * (cone - chalf * l_val) * (3 * Ylm_rescale + rhat_dot_G);
        }
        local_val += Npad;
        local_grad += Npad;
        local_lapl += Npad;
      }
    }
    else
    {
      std::cout << "Warning: an electron is on top of an ion!" << std::endl;
      // strictly zero

#pragma omp simd aligned(val, lapl, local_val, local_lapl : QMC_SIMD_ALIGNMENT)
      for (size_t ib = 0; ib < myV.size(); ib++)
      {
        // value
        val[ib] = Ylm_v[0] * local_val[ib];

        // laplacian
        lapl[ib] = local_lapl[ib] * static_cast<ST>(3) * Ylm_v[0];
      }
      local_val += Npad;
      local_grad += Npad;
      local_lapl += Npad;
      if (lm_tot > 0)
      {
        //std::cout << std::endl;
        for (size_t lm = 1; lm < 4; lm++)
        {
#pragma omp simd aligned(g0, g1, g2, local_grad : QMC_SIMD_ALIGNMENT)
          for (size_t ib = 0; ib < myV.size(); ib++)
          {
            const ST local_g = local_grad[ib];
            // grad
            g0[ib] += local_g * Ylm_gx[lm];
            g1[ib] += local_g * Ylm_gy[lm];
            g2[ib] += local_g * Ylm_gz[lm];
          }
          local_grad += Npad;
        }
      }
    }
  }

  template<typename VV, typename GV, typename HT>
  void evaluate_vgh(const ST& r, const PointType& dr, VV& myV, GV& myG, HT& myH)
  {
    //Needed to do tensor product here
    APP_ABORT("AtomicOrbitals::evaluate_vgh");
  }

  void storeParamsBeforeRotation()
  {
    const auto spline_ptr     = SplineInst->getSplinePtr();
    const auto coefs_tot_size = spline_ptr->coefs_size;
    coef_copy_                = std::make_shared<std::vector<ST>>(coefs_tot_size);
    std::copy_n(spline_ptr->coefs, coefs_tot_size, coef_copy_->begin());
  }

  template<typename VM>
  void applyRotation(const VM& rot_mat, bool use_stored_copy)
  {
    // SplineInst is a MultiBspline. See src/spline2/MultiBspline.hpp
    const auto spline_ptr = SplineInst->getSplinePtr();
    assert(spline_ptr != nullptr);
    const auto spl_coefs      = spline_ptr->coefs;
    const auto Nsplines       = spline_ptr->num_splines; // May include padding
    const auto coefs_tot_size = spline_ptr->coefs_size;
    const auto BasisSetSize   = coefs_tot_size / Nsplines;
    const auto TrueNOrbs      = rot_mat.size1(); // == Nsplines - padding

    if (!use_stored_copy)
    {
      assert(coef_copy_ != nullptr);
      std::copy_n(spl_coefs, coefs_tot_size, coef_copy_->begin());
    }

    // There are TrueNOrb splines for every l,m pair, and also padding
#ifdef QMC_COMPLEX
    for (size_t lidx = 0; lidx < lm_tot; lidx++)
      for (size_t i = 0; i < BasisSetSize; i++)
        for (size_t j = 0; j < TrueNOrbs; j++)
        {
          const auto cur_elem = i * Nsplines + lidx * Npad + 2 * j;
          ST newval_r{0.};
          ST newval_i{0.};
          for (size_t k = 0; k < TrueNOrbs; k++)
          {
            const auto index = i * Nsplines + lidx * Npad + 2 * k;
            ST zr            = (*coef_copy_)[index];
            ST zi            = (*coef_copy_)[index + 1];
            ST wr            = rot_mat[k][j].real();
            ST wi            = rot_mat[k][j].imag();
            newval_r += zr * wr - zi * wi;
            newval_i += zr * wi + zi * wr;
          }
          spl_coefs[cur_elem]     = newval_r;
          spl_coefs[cur_elem + 1] = newval_i;
        }
#else
    for (size_t lidx = 0; lidx < lm_tot; lidx++)
      for (size_t i = 0; i < BasisSetSize; i++)
        for (size_t j = 0; j < TrueNOrbs; j++)
        {
          const auto cur_elem = i * Nsplines + lidx * Npad + j;
          ST newval{0.};
          for (size_t k = 0; k < TrueNOrbs; k++)
          {
            const auto index = i * Nsplines + lidx * Npad + k;
            newval += (*coef_copy_)[index] * rot_mat[k][j];
          }
          spl_coefs[cur_elem] = newval;
        }
#endif
  }
};

template<typename ST>
class HybridRepCenterOrbitals
{
public:
  static const int D = 3;
  using PointType    = typename AtomicOrbitals<ST>::PointType;
  using RealType     = typename DistanceTable::RealType;
  using PosType      = typename DistanceTable::PosType;

  enum class Region
  {
    INSIDE, // within the buffer shell
    BUFFER, // in the buffer region
    INTER   // interstitial area
  };

  struct LocationSmoothingInfo
  {
    ///r from distance table
    RealType dist_r;
    ///dr from distance table
    PosType dist_dr;
    ///for APBC
    PointType r_image;
    /// region of the location
    Region region;
    ///smooth function value
    RealType f;
    ///smooth function first derivative
    RealType df_dr;
    ///smooth function second derivative
    RealType d2f_dr2;
  };

private:
  ///atomic centers
  std::vector<AtomicOrbitals<ST>> AtomicCenters;
  ///table index
  int myTableID;
  ///mapping supercell to primitive cell
  std::vector<int> Super2Prim;
  ///smoothing schemes
  enum class smoothing_schemes
  {
    CONSISTENT = 0,
    SMOOTHALL,
    SMOOTHPARTIAL
  } smooth_scheme;
  /// smoothing function
  smoothing_functions smooth_func_id;

  /// select a region (within the buffer shell, in the buffer, interstitial region) and compute the smoothing function if in the buffer.
  inline void selectRegionAndComputeSmoothing(const ST& cutoff_buffer,
                                              const ST& cutoff,
                                              LocationSmoothingInfo& info) const
  {
    const RealType r = info.dist_r;
    if (r < cutoff_buffer)
      info.region = Region::INSIDE;
    else if (r < cutoff)
    {
      constexpr RealType cone(1);
      const RealType scale = cone / (cutoff - cutoff_buffer);
      const RealType x     = (r - cutoff_buffer) * scale;
      info.f               = smoothing(smooth_func_id, x, info.df_dr, info.d2f_dr2);
      info.df_dr *= scale;
      info.d2f_dr2 *= scale * scale;
      info.region = Region::BUFFER;
    }
    else
      info.region = Region::INTER;
  }

public:
  HybridRepCenterOrbitals() {}

  void storeParamsBeforeRotation()
  {
    for (auto& atomic_center : AtomicCenters)
      atomic_center.storeParamsBeforeRotation();
  }

  template<typename VM>
  void applyRotation(const VM& rot_mat, bool use_stored_copy)
  {
    for (auto& atomic_center : AtomicCenters)
      atomic_center.applyRotation(rot_mat, use_stored_copy);
  }

  void set_info(const ParticleSet& ions, ParticleSet& els, const std::vector<int>& mapping)
  {
    myTableID  = els.addTable(ions, DTModes::NEED_VP_FULL_TABLE_ON_HOST);
    Super2Prim = mapping;
  }

  inline void resizeStorage(size_t Nb)
  {
    size_t SplineCoefsBytes = 0;

    for (int ic = 0; ic < AtomicCenters.size(); ic++)
    {
      AtomicCenters[ic].resizeStorage(Nb);
      SplineCoefsBytes += AtomicCenters[ic].getSplineSizeInBytes();
    }

    app_log() << "MEMORY " << SplineCoefsBytes / (1 << 20) << " MB allocated "
              << "for the atomic radial splines in hybrid orbital representation" << std::endl;
  }

  void bcast_tables(Communicate* comm)
  {
    for (int ic = 0; ic < AtomicCenters.size(); ic++)
      AtomicCenters[ic].bcast_tables(comm);
  }

  void gather_atomic_tables(Communicate* comm, std::vector<int>& offset)
  {
    if (comm->size() == 1)
      return;
    for (int ic = 0; ic < AtomicCenters.size(); ic++)
      AtomicCenters[ic].gather_tables(comm, offset);
  }

  inline void flush_zero()
  {
    for (int ic = 0; ic < AtomicCenters.size(); ic++)
      AtomicCenters[ic].flush_zero();
  }

  bool read_splines(hdf_archive& h5f)
  {
    bool success = true;
    size_t ncenter;

    try
    {
      h5f.push("atomic_centers", false);
    }
    catch (...)
    {
      success = false;
    }
    success = success && h5f.readEntry(ncenter, "number_of_centers");
    if (!success)
      return success;
    if (ncenter != AtomicCenters.size())
      success = false;
    // read splines of each center
    for (int ic = 0; ic < AtomicCenters.size(); ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      try
      {
        h5f.push(gname.str().c_str(), false);
      }
      catch (...)
      {
        success = false;
      }
      success = success && AtomicCenters[ic].read_splines(h5f);
      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  bool write_splines(hdf_archive& h5f)
  {
    bool success = true;
    int ncenter  = AtomicCenters.size();
    try
    {
      h5f.push("atomic_centers", true);
    }
    catch (...)
    {
      success = false;
    }
    success = success && h5f.writeEntry(ncenter, "number_of_centers");
    // write splines of each center
    for (int ic = 0; ic < AtomicCenters.size(); ic++)
    {
      std::ostringstream gname;
      gname << "center_" << ic;
      try
      {
        h5f.push(gname.str().c_str(), true);
      }
      catch (...)
      {
        success = false;
      }
      success = success && AtomicCenters[ic].write_splines(h5f);
      h5f.pop();
    }
    h5f.pop();
    return success;
  }

  template<typename Cell>
  inline int get_bc_sign(const PointType& r,
                         const PointType& r_image,
                         const Cell& PrimLattice,
                         TinyVector<int, D>& HalfG) const
  {
    int bc_sign          = 0;
    PointType shift_unit = PrimLattice.toUnit(r - r_image);
    for (int i = 0; i < D; i++)
    {
      ST img = round(shift_unit[i]);
      bc_sign += HalfG[i] * (int)img;
    }
    return bc_sign;
  }

  //evaluate only V
  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& myV, LocationSmoothingInfo& info)
  {
    const auto& ei_dist  = P.getDistTableAB(myTableID);
    const int center_idx = ei_dist.get_first_neighbor(iat, info.dist_r, info.dist_dr, P.getActivePtcl() == iat);
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    selectRegionAndComputeSmoothing(myCenter.getCutoffBuffer(), myCenter.getCutoff(), info);
    if (info.region != Region::INTER)
    {
      PointType dr(-info.dist_dr[0], -info.dist_dr[1], -info.dist_dr[2]);
      info.r_image = myCenter.getCenterPos() + dr;
      myCenter.evaluate_v(info.dist_r, dr, myV);
    }
  }

  /* check if the batched algorithm is safe to operate
   * @param VP virtual particle set
   * @return true if it is safe
   *
   * When the reference electron in the NLPP evaluation has a distance larger than the non overlapping radius of the reference center.
   * Some qudrature points may get its SPOs evaluated from the nearest center which is not the reference center.
   * The batched algorthm forces the evaluation on the reference center and introduce some error.
   * In this case, the non-batched algorithm should be used.
   */
  bool is_VP_batching_safe(const VirtualParticleSet& VP) const
  {
    const int center_idx = VP.refSourcePtcl;
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    return VP.getRefPS().getDistTableAB(myTableID).getDistRow(VP.refPtcl)[center_idx] <
        myCenter.getNonOverlappingRadius();
  }

  // C2C, C2R cases
  template<typename VM>
  inline void evaluateValuesC2X(const VirtualParticleSet& VP, VM& multi_myV, LocationSmoothingInfo& info)
  {
    const int center_idx = VP.refSourcePtcl;
    info.dist_r          = VP.getRefPS().getDistTableAB(myTableID).getDistRow(VP.refPtcl)[center_idx];
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    selectRegionAndComputeSmoothing(myCenter.getCutoffBuffer(), myCenter.getCutoff(), info);
    if (info.region != Region::INTER)
      myCenter.evaluateValues(VP.getDistTableAB(myTableID).getDisplacements(), center_idx, info.dist_r, multi_myV);
  }

  // R2R case
  template<typename VM, typename Cell, typename SV>
  inline void evaluateValuesR2R(const VirtualParticleSet& VP,
                                const Cell& PrimLattice,
                                TinyVector<int, D>& HalfG,
                                VM& multi_myV,
                                SV& bc_signs,
                                LocationSmoothingInfo& info)
  {
    const int center_idx = VP.refSourcePtcl;
    info.dist_r          = VP.getRefPS().getDistTableAB(myTableID).getDistRow(VP.refPtcl)[center_idx];
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    selectRegionAndComputeSmoothing(myCenter.getCutoffBuffer(), myCenter.getCutoff(), info);
    if (info.region != Region::INTER)
    {
      const auto& displ = VP.getDistTableAB(myTableID).getDisplacements();
      for (int ivp = 0; ivp < VP.getTotalNum(); ivp++)
        bc_signs[ivp] = get_bc_sign(VP.R[ivp], myCenter.getCenterPos() - displ[ivp][center_idx], PrimLattice, HalfG);
      myCenter.evaluateValues(displ, center_idx, info.dist_r, multi_myV);
    }
  }

  //evaluate only VGL
  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& myV, GV& myG, VV& myL, LocationSmoothingInfo& info)
  {
    const auto& ei_dist  = P.getDistTableAB(myTableID);
    const int center_idx = ei_dist.get_first_neighbor(iat, info.dist_r, info.dist_dr, P.getActivePtcl() == iat);
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    selectRegionAndComputeSmoothing(myCenter.getCutoffBuffer(), myCenter.getCutoff(), info);
    if (info.region != Region::INTER)
    {
      const PointType dr(-info.dist_dr[0], -info.dist_dr[1], -info.dist_dr[2]);
      info.r_image = myCenter.getCenterPos() + dr;
      myCenter.evaluate_vgl(info.dist_r, dr, myV, myG, myL);
    }
  }

  //evaluate only VGH
  template<typename VV, typename GV, typename HT>
  inline void evaluate_vgh(const ParticleSet& P, const int iat, VV& myV, GV& myG, HT& myH, LocationSmoothingInfo& info)
  {
    const auto& ei_dist  = P.getDistTableAB(myTableID);
    const int center_idx = ei_dist.get_first_neighbor(iat, info.dist_r, info.dist_dr, P.getActivePtcl() == iat);
    auto& myCenter       = AtomicCenters[Super2Prim[center_idx]];
    selectRegionAndComputeSmoothing(myCenter.getCutoffBuffer(), myCenter.getCutoff(), info);
    if (info.region != Region::INTER)
    {
      const PointType dr(-info.dist_dr[0], -info.dist_dr[1], -info.dist_dr[2]);
      info.r_image = myCenter.getCenterPos() + dr;
      myCenter.evaluate_vgh(info.dist_r, dr, myV, myG, myH);
    }
  }

  // interpolate buffer region, value only
  template<typename VV>
  inline void interpolate_buffer_v(VV& psi, const VV& psi_AO, const RealType f) const
  {
    constexpr RealType cone(1);
    for (size_t i = 0; i < psi.size(); i++)
      psi[i] = psi_AO[i] * f + psi[i] * (cone - f);
  }

  // interpolate buffer region, value, gradients and laplacian
  template<typename VV, typename GV>
  inline void interpolate_buffer_vgl(VV& psi,
                                     GV& dpsi,
                                     VV& d2psi,
                                     const VV& psi_AO,
                                     const GV& dpsi_AO,
                                     const VV& d2psi_AO,
                                     const LocationSmoothingInfo& info) const
  {
    constexpr RealType cone(1), ctwo(2);
    const RealType rinv(1.0 / info.dist_r);
    auto& dist_dr = info.dist_dr;
    auto& f       = info.f;
    auto& df_dr   = info.df_dr;
    auto& d2f_dr2 = info.d2f_dr2;
    if (smooth_scheme == smoothing_schemes::CONSISTENT)
      for (size_t i = 0; i < psi.size(); i++)
      { // psi, dpsi, d2psi are all consistent
        d2psi[i] = d2psi_AO[i] * f + d2psi[i] * (cone - f) + df_dr * rinv * ctwo * dot(dpsi[i] - dpsi_AO[i], dist_dr) +
            (psi_AO[i] - psi[i]) * (d2f_dr2 + ctwo * rinv * df_dr);
        dpsi[i] = dpsi_AO[i] * f + dpsi[i] * (cone - f) + df_dr * rinv * dist_dr * (psi[i] - psi_AO[i]);
        psi[i]  = psi_AO[i] * f + psi[i] * (cone - f);
      }
    else if (smooth_scheme == smoothing_schemes::SMOOTHALL)
      for (size_t i = 0; i < psi.size(); i++)
      {
        d2psi[i] = d2psi_AO[i] * f + d2psi[i] * (cone - f);
        dpsi[i]  = dpsi_AO[i] * f + dpsi[i] * (cone - f);
        psi[i]   = psi_AO[i] * f + psi[i] * (cone - f);
      }
    else if (smooth_scheme == smoothing_schemes::SMOOTHPARTIAL)
      for (size_t i = 0; i < psi.size(); i++)
      { // dpsi, d2psi are consistent but psi is not.
        d2psi[i] = d2psi_AO[i] * f + d2psi[i] * (cone - f) + df_dr * rinv * ctwo * dot(dpsi[i] - dpsi_AO[i], dist_dr);
        dpsi[i]  = dpsi_AO[i] * f + dpsi[i] * (cone - f);
        psi[i]   = psi_AO[i] * f + psi[i] * (cone - f);
      }
    else
      throw std::runtime_error("Unknown smooth scheme!");
  }

  template<class BSPLINESPO>
  friend class qmcplusplus::HybridRepSetReader;
};

extern template class AtomicOrbitals<float>;
extern template class AtomicOrbitals<double>;
extern template class HybridRepCenterOrbitals<float>;
extern template class HybridRepCenterOrbitals<double>;
} // namespace qmcplusplus
#endif
