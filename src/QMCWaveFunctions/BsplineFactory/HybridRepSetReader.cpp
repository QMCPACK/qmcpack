//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file
 *
 * derived from SplineSetReader
 */

#include "HybridRepSetReader.h"
#include "Numerics/Quadrature.h"
#include "Numerics/Bessel.h"
#include "OhmmsData/AttributeSet.h"
#include "CPU/math.hpp"
#include "Concurrency/OpenMP.h"
#include <Timer.h>
#include "OneSplineOrbData.hpp"

namespace qmcplusplus
{
template<typename ST, typename LT>
struct Gvectors
{
  using PosType   = TinyVector<ST, 3>;
  using ValueType = std::complex<ST>;

  const LT& Lattice;
  std::vector<PosType> gvecs_cart; //Cartesian.
  std::vector<ST> gmag;
  const size_t NumGvecs;

  Gvectors(const std::vector<TinyVector<int, 3>>& gvecs_in,
           const LT& Lattice_in,
           const TinyVector<int, 3>& HalfG,
           size_t first,
           size_t last)
      : Lattice(Lattice_in), NumGvecs(last - first)
  {
    gvecs_cart.resize(NumGvecs);
    gmag.resize(NumGvecs);
#pragma omp parallel for
    for (size_t ig = 0; ig < NumGvecs; ig++)
    {
      TinyVector<ST, 3> gvec_shift;
      gvec_shift     = gvecs_in[ig + first] + HalfG * 0.5;
      gvecs_cart[ig] = Lattice.k_cart(gvec_shift);
      gmag[ig]       = std::sqrt(dot(gvecs_cart[ig], gvecs_cart[ig]));
    }
  }

  template<typename YLM_ENGINE, typename VVT>
  void calc_Ylm_G(const size_t ig, YLM_ENGINE& Ylm, VVT& YlmG) const
  {
    PosType Ghat(0.0, 0.0, 1.0);
    if (gmag[ig] > 0)
      Ghat = gvecs_cart[ig] / gmag[ig];
    Ylm.evaluateV(Ghat[0], Ghat[1], Ghat[2], YlmG.data());
  }

  template<typename VVT>
  inline void calc_jlm_G(const int lmax, ST& r, const size_t ig, VVT& j_lm_G) const
  {
    bessel_steed_array_cpu(lmax, gmag[ig] * r, j_lm_G.data());
    for (size_t l = lmax; l > 0; l--)
      for (size_t lm = l * l; lm < (l + 1) * (l + 1); lm++)
        j_lm_G[lm] = j_lm_G[l];
  }

  template<typename PT, typename VT>
  inline void calc_phase_shift(const PT& RSoA, const size_t ig, VT& phase_shift_real, VT& phase_shift_imag) const
  {
    const ST* restrict px = RSoA.data(0);
    const ST* restrict py = RSoA.data(1);
    const ST* restrict pz = RSoA.data(2);
    ST* restrict v_r      = phase_shift_real.data();
    ST* restrict v_i      = phase_shift_imag.data();
    const ST& gv_x        = gvecs_cart[ig][0];
    const ST& gv_y        = gvecs_cart[ig][1];
    const ST& gv_z        = gvecs_cart[ig][2];

#pragma omp simd aligned(px, py, pz, v_r, v_i : QMC_SIMD_ALIGNMENT)
    for (size_t iat = 0; iat < RSoA.size(); iat++)
      qmcplusplus::sincos(px[iat] * gv_x + py[iat] * gv_y + pz[iat] * gv_z, v_i + iat, v_r + iat);
  }

  template<typename PT>
  ValueType evaluate_psi_r(const Vector<std::complex<double>>& cG, const PT& pos)
  {
    assert(cG.size() == NumGvecs);
    std::complex<ST> val(0.0, 0.0);
    for (size_t ig = 0; ig < NumGvecs; ig++)
    {
      ST s, c;
      qmcplusplus::sincos(dot(gvecs_cart[ig], pos), &s, &c);
      ValueType pw0(c, s);
      val += cG[ig] * pw0;
    }
    return val;
  }

  template<typename PT>
  void evaluate_psi_r(const Vector<std::complex<double>>& cG, const PT& pos, ValueType& phi, ValueType& d2phi)
  {
    assert(cG.size() == NumGvecs);
    d2phi = phi = 0.0;
    for (size_t ig = 0; ig < NumGvecs; ig++)
    {
      ST s, c;
      qmcplusplus::sincos(dot(gvecs_cart[ig], pos), &s, &c);
      ValueType pw0(c, s);
      phi += cG[ig] * pw0;
      d2phi += cG[ig] * pw0 * (-dot(gvecs_cart[ig], gvecs_cart[ig]));
    }
  }

  double evaluate_KE(const Vector<std::complex<double>>& cG)
  {
    assert(cG.size() == NumGvecs);
    double KE = 0;
    for (size_t ig = 0; ig < NumGvecs; ig++)
      KE += dot(gvecs_cart[ig], gvecs_cart[ig]) * (cG[ig].real() * cG[ig].real() + cG[ig].imag() * cG[ig].imag());
    return KE / 2.0;
  }
};

template<typename SA>
HybridRepSetReader<SA>::HybridRepSetReader(EinsplineSetBuilder* e) : BsplineReader(e), spline_reader_(e)
{}

template<typename SA>
std::unique_ptr<SPOSet> HybridRepSetReader<SA>::create_spline_set(const std::string& my_name,
                                                                  int spin,
                                                                  const BandInfoGroup& bandgroup)
{
  auto bspline = std::make_unique<SA>(my_name);
  app_log() << "  ClassName = " << bspline->getClassName() << std::endl;
  // set info for Hybrid
  initialize_hybridrep_atomic_centers(*bspline);
  bool foundspline = spline_reader_.createSplineDataSpaceLookforDumpFile(bandgroup, *bspline);
  typename SA::HYBRIDBASE& hybrid_center_orbs = *bspline;
  hybrid_center_orbs.resizeStorage(bspline->myV.size());
  if (foundspline)
  {
    Timer now;
    hdf_archive h5f(myComm);
    const auto splinefile = getSplineDumpFileName(bandgroup);
    h5f.open(splinefile, H5F_ACC_RDONLY);
    foundspline = bspline->read_splines(h5f);
    if (foundspline)
      app_log() << "  Successfully restored 3D B-spline coefficients from " << splinefile << ". The reading time is "
                << now.elapsed() << " sec." << std::endl;
  }

  if (!foundspline)
  {
    hybrid_center_orbs.flush_zero();
    initialize_hybrid_pio_gather(spin, bandgroup, *bspline);

    if (saveSplineCoefs && myComm->rank() == 0)
    {
      Timer now;
      const std::string splinefile(getSplineDumpFileName(bandgroup));
      hdf_archive h5f;
      h5f.create(splinefile);
      std::string classname = bspline->getClassName();
      h5f.write(classname, "class_name");
      int sizeD = sizeof(DataType);
      h5f.write(sizeD, "sizeof");
      bspline->write_splines(h5f);
      h5f.close();
      app_log() << "  Stored spline coefficients in " << splinefile << " for potential reuse. The writing time is "
                << now.elapsed() << " sec." << std::endl;
    }
  }

  {
    Timer now;
    bspline->bcast_tables(myComm);
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
  }
  return bspline;
}

template<typename SA>
void HybridRepSetReader<SA>::initialize_hybridrep_atomic_centers(SA& bspline) const
{
  auto& mybuilder = spline_reader_.mybuilder;

  OhmmsAttributeSet a;
  std::string scheme_name("Consistent");
  std::string s_function_name("LEKS2018");
  a.add(scheme_name, "smoothing_scheme");
  a.add(s_function_name, "smoothing_function");
  a.put(mybuilder->XMLRoot);
  // assign smooth_scheme
  if (scheme_name == "Consistent")
    bspline.smooth_scheme = SA::smoothing_schemes::CONSISTENT;
  else if (scheme_name == "SmoothAll")
    bspline.smooth_scheme = SA::smoothing_schemes::SMOOTHALL;
  else if (scheme_name == "SmoothPartial")
    bspline.smooth_scheme = SA::smoothing_schemes::SMOOTHPARTIAL;
  else
    APP_ABORT("initialize_hybridrep_atomic_centers wrong smoothing_scheme name! Only allows Consistent, SmoothAll or "
              "SmoothPartial.");

  // assign smooth_function
  if (s_function_name == "LEKS2018")
    bspline.smooth_func_id = smoothing_functions::LEKS2018;
  else if (s_function_name == "coscos")
    bspline.smooth_func_id = smoothing_functions::COSCOS;
  else if (s_function_name == "linear")
    bspline.smooth_func_id = smoothing_functions::LINEAR;
  else
    APP_ABORT(
        "initialize_hybridrep_atomic_centers wrong smoothing_function name! Only allows LEKS2018, coscos or linear.");
  app_log() << "Hybrid orbital representation uses " << scheme_name << " smoothing scheme and " << s_function_name
            << " smoothing function." << std::endl;

  bspline.set_info(*(mybuilder->SourcePtcl), mybuilder->TargetPtcl, mybuilder->Super2Prim);
  auto& centers = bspline.AtomicCenters;
  auto& ACInfo  = mybuilder->AtomicCentersInfo;
  // load atomic center info only when it is not initialized
  if (centers.size() == 0)
  {
    bool success = true;
    app_log() << "Reading atomic center info for hybrid representation" << std::endl;
    for (int center_idx = 0; center_idx < ACInfo.Ncenters; center_idx++)
    {
      const int my_GroupID = ACInfo.GroupID[center_idx];
      if (ACInfo.cutoff[center_idx] < 0)
      {
        app_error() << "Hybrid orbital representation needs parameter 'cutoff_radius' for atom " << center_idx
                    << std::endl;
        success = false;
      }

      if (ACInfo.inner_cutoff[center_idx] < 0)
      {
        const double inner_cutoff = std::max(ACInfo.cutoff[center_idx] - 0.3, 0.0);
        app_log() << "Hybrid orbital representation setting 'inner_cutoff' to " << inner_cutoff << " for group "
                  << my_GroupID << " as atom " << center_idx << std::endl;
        // overwrite the inner_cutoff of all the atoms of the same species
        for (int id = 0; id < ACInfo.Ncenters; id++)
          if (my_GroupID == ACInfo.GroupID[id])
            ACInfo.inner_cutoff[id] = inner_cutoff;
      }
      else if (ACInfo.inner_cutoff[center_idx] > ACInfo.cutoff[center_idx])
      {
        app_error() << "Hybrid orbital representation 'inner_cutoff' must be smaller than 'spline_radius' for atom "
                    << center_idx << std::endl;
        success = false;
      }

      if (ACInfo.cutoff[center_idx] > 0)
      {
        if (ACInfo.lmax[center_idx] < 0)
        {
          app_error() << "Hybrid orbital representation needs parameter 'lmax' for atom " << center_idx << std::endl;
          success = false;
        }

        if (ACInfo.spline_radius[center_idx] < 0 && ACInfo.spline_npoints[center_idx] < 0)
        {
          app_log() << "Parameters 'spline_radius' and 'spline_npoints' for group " << my_GroupID << " as atom "
                    << center_idx << " are not specified." << std::endl;
          const double delta     = std::min(0.02, ACInfo.cutoff[center_idx] / 4.0);
          const int n_grid_point = std::ceil((ACInfo.cutoff[center_idx] + 1e-4) / delta) + 3;
          for (int id = 0; id < ACInfo.Ncenters; id++)
            if (my_GroupID == ACInfo.GroupID[id])
            {
              ACInfo.spline_npoints[id] = n_grid_point;
              ACInfo.spline_radius[id]  = (n_grid_point - 1) * delta;
            }
          app_log() << "  Based on default grid point distance " << delta << std::endl;
          app_log() << "  Setting 'spline_npoints' to " << ACInfo.spline_npoints[center_idx] << std::endl;
          app_log() << "  Setting 'spline_radius' to " << ACInfo.spline_radius[center_idx] << std::endl;
        }
        else
        {
          if (ACInfo.spline_radius[center_idx] < 0)
          {
            app_error() << "Hybrid orbital representation needs parameter 'spline_radius' for atom " << center_idx
                        << std::endl;
            success = false;
          }

          if (ACInfo.spline_npoints[center_idx] < 0)
          {
            app_error() << "Hybrid orbital representation needs parameter 'spline_npoints' for atom " << center_idx
                        << std::endl;
            success = false;
          }
        }

        // check maximally allowed cutoff_radius
        double max_allowed_cutoff = ACInfo.spline_radius[center_idx] -
            2.0 * ACInfo.spline_radius[center_idx] / (ACInfo.spline_npoints[center_idx] - 1);
        if (success && ACInfo.cutoff[center_idx] > max_allowed_cutoff)
        {
          app_error() << "Hybrid orbital representation requires cutoff_radius<=" << max_allowed_cutoff
                      << " calculated by spline_radius-2*spline_radius/(spline_npoints-1) for atom " << center_idx
                      << std::endl;
          success = false;
        }
      }
      else
      {
        // no atomic regions for this atom type
        ACInfo.spline_radius[center_idx]  = 0.0;
        ACInfo.spline_npoints[center_idx] = 0;
        ACInfo.lmax[center_idx]           = 0;
      }
    }

    if (!success)
      spline_reader_.myComm->barrier_and_abort(
          "initialize_hybridrep_atomic_centers Failed to initialize atomic centers "
          "in hybrid orbital representation!");

    for (int center_idx = 0; center_idx < ACInfo.Ncenters; center_idx++)
    {
      AtomicOrbitals<DataType> oneCenter(ACInfo.lmax[center_idx]);
      oneCenter.set_info(ACInfo.ion_pos[center_idx], ACInfo.cutoff[center_idx], ACInfo.inner_cutoff[center_idx],
                         ACInfo.spline_radius[center_idx], ACInfo.non_overlapping_radius[center_idx],
                         ACInfo.spline_npoints[center_idx]);
      centers.push_back(oneCenter);
    }
  }
}

template<typename SA>
void HybridRepSetReader<SA>::create_atomic_centers_Gspace(const Vector<std::complex<double>>& cG,
                                                          Communicate& band_group_comm,
                                                          const int iorb,
                                                          const std::complex<double>& rotate_phase,
                                                          SA& bspline) const
{
  auto& mybuilder       = spline_reader_.mybuilder;
  double rotate_phase_r = rotate_phase.real();
  double rotate_phase_i = rotate_phase.imag();

  band_group_comm.bcast(rotate_phase_r);
  band_group_comm.bcast(rotate_phase_i);
  //distribute G-vectors over processor groups
  const int Ngvecs      = mybuilder->Gvecs[0].size();
  const int Nprocs      = band_group_comm.size();
  const int Ngvecgroups = std::min(Ngvecs, Nprocs);
  Communicate gvec_group_comm(band_group_comm, Ngvecgroups);
  std::vector<int> gvec_groups(Ngvecgroups + 1, 0);
  FairDivideLow(Ngvecs, Ngvecgroups, gvec_groups);
  const int gvec_first = gvec_groups[gvec_group_comm.getGroupID()];
  const int gvec_last  = gvec_groups[gvec_group_comm.getGroupID() + 1];

  // prepare Gvecs Ylm(G)
  using UnitCellType = typename EinsplineSetBuilder::UnitCellType;
  Gvectors<double, UnitCellType> Gvecs(mybuilder->Gvecs[0], mybuilder->PrimCell, bspline.HalfG, gvec_first, gvec_last);
  // if(band_group_comm.isGroupLeader()) std::cout << "print band=" << iorb << " KE=" << Gvecs.evaluate_KE(cG) << std::endl;

  std::vector<AtomicOrbitals<DataType>>& centers = bspline.AtomicCenters;
  app_log() << "Transforming band " << iorb << " on Rank 0" << std::endl;
  // collect atomic centers by group
  std::vector<int> uniq_species;
  for (int center_idx = 0; center_idx < centers.size(); center_idx++)
  {
    auto& ACInfo         = mybuilder->AtomicCentersInfo;
    const int my_GroupID = ACInfo.GroupID[center_idx];
    int found_idx        = -1;
    for (size_t idx = 0; idx < uniq_species.size(); idx++)
      if (my_GroupID == uniq_species[idx])
      {
        found_idx = idx;
        break;
      }
    if (found_idx < 0)
      uniq_species.push_back(my_GroupID);
  }
  // construct group list
  std::vector<std::vector<int>> group_list(uniq_species.size());
  for (int center_idx = 0; center_idx < centers.size(); center_idx++)
  {
    auto& ACInfo         = mybuilder->AtomicCentersInfo;
    const int my_GroupID = ACInfo.GroupID[center_idx];
    for (size_t idx = 0; idx < uniq_species.size(); idx++)
      if (my_GroupID == uniq_species[idx])
      {
        group_list[idx].push_back(center_idx);
        break;
      }
  }

  for (int group_idx = 0; group_idx < group_list.size(); group_idx++)
  {
    const auto& mygroup        = group_list[group_idx];
    const double spline_radius = centers[mygroup[0]].getSplineRadius();
    const int spline_npoints   = centers[mygroup[0]].getSplineNpoints();
    const int lmax             = centers[mygroup[0]].getLmax();
    const double delta         = spline_radius / static_cast<double>(spline_npoints - 1);
    const int lm_tot           = (lmax + 1) * (lmax + 1);
    const size_t natoms        = mygroup.size();
    const int policy           = lm_tot > natoms ? 0 : 1;

    std::vector<std::complex<double>> i_power(lm_tot);
    // rotate phase is introduced here.
    std::complex<double> i_temp(rotate_phase_r, rotate_phase_i);
    for (size_t l = 0; l <= lmax; l++)
    {
      for (size_t lm = l * l; lm < (l + 1) * (l + 1); lm++)
        i_power[lm] = i_temp;
      i_temp *= std::complex<double>(0.0, 1.0);
    }

    std::vector<Matrix<double>> all_vals(natoms);
    std::vector<std::vector<aligned_vector<double>>> vals_local(spline_npoints * omp_get_max_threads());
    VectorSoaContainer<double, 3> myRSoA(natoms);
    for (size_t idx = 0; idx < natoms; idx++)
    {
      all_vals[idx].resize(spline_npoints, lm_tot * 2);
      all_vals[idx] = 0.0;
      myRSoA(idx)   = centers[mygroup[idx]].getCenterPos();
    }

#pragma omp parallel
    {
      const size_t tid = omp_get_thread_num();
      const size_t nt  = omp_get_num_threads();

      for (int ip = 0; ip < spline_npoints; ip++)
      {
        const size_t ip_idx = tid * spline_npoints + ip;
        if (policy == 1)
        {
          vals_local[ip_idx].resize(lm_tot * 2);
          for (size_t lm = 0; lm < lm_tot * 2; lm++)
          {
            auto& vals = vals_local[ip_idx][lm];
            vals.resize(natoms);
            std::fill(vals.begin(), vals.end(), 0.0);
          }
        }
        else
        {
          vals_local[ip_idx].resize(natoms * 2);
          for (size_t iat = 0; iat < natoms * 2; iat++)
          {
            auto& vals = vals_local[ip_idx][iat];
            vals.resize(lm_tot);
            std::fill(vals.begin(), vals.end(), 0.0);
          }
        }
      }

      const size_t size_pw_tile = 32;
      const size_t num_pw_tiles = (Gvecs.NumGvecs + size_pw_tile - 1) / size_pw_tile;
      aligned_vector<double> j_lm_G(lm_tot, 0.0);
      std::vector<aligned_vector<double>> phase_shift_r(size_pw_tile);
      std::vector<aligned_vector<double>> phase_shift_i(size_pw_tile);
      std::vector<aligned_vector<double>> YlmG(size_pw_tile);
      for (size_t ig = 0; ig < size_pw_tile; ig++)
      {
        phase_shift_r[ig].resize(natoms);
        phase_shift_i[ig].resize(natoms);
        YlmG[ig].resize(lm_tot);
      }
      SoaSphericalTensor<double> Ylm(lmax);

#pragma omp for
      for (size_t tile_id = 0; tile_id < num_pw_tiles; tile_id++)
      {
        const size_t ig_first = tile_id * size_pw_tile;
        const size_t ig_last  = std::min((tile_id + 1) * size_pw_tile, Gvecs.NumGvecs);
        for (size_t ig = ig_first; ig < ig_last; ig++)
        {
          const size_t ig_local = ig - ig_first;
          // calculate phase shift for all the centers of this group
          Gvecs.calc_phase_shift(myRSoA, ig, phase_shift_r[ig_local], phase_shift_i[ig_local]);
          Gvecs.calc_Ylm_G(ig, Ylm, YlmG[ig_local]);
        }

        for (int ip = 0; ip < spline_npoints; ip++)
        {
          double r            = delta * static_cast<double>(ip);
          const size_t ip_idx = tid * spline_npoints + ip;

          for (size_t ig = ig_first; ig < ig_last; ig++)
          {
            const size_t ig_local = ig - ig_first;
            // calculate spherical bessel function
            Gvecs.calc_jlm_G(lmax, r, ig, j_lm_G);
            for (size_t lm = 0; lm < lm_tot; lm++)
              j_lm_G[lm] *= YlmG[ig_local][lm];

            const double cG_r = cG[ig + gvec_first].real();
            const double cG_i = cG[ig + gvec_first].imag();
            if (policy == 1)
            {
              for (size_t lm = 0; lm < lm_tot; lm++)
              {
                double* restrict vals_r         = vals_local[ip_idx][lm * 2].data();
                double* restrict vals_i         = vals_local[ip_idx][lm * 2 + 1].data();
                const double* restrict ps_r_ptr = phase_shift_r[ig_local].data();
                const double* restrict ps_i_ptr = phase_shift_i[ig_local].data();
                double cG_j_r                   = cG_r * j_lm_G[lm];
                double cG_j_i                   = cG_i * j_lm_G[lm];
#pragma omp simd aligned(vals_r, vals_i, ps_r_ptr, ps_i_ptr : QMC_SIMD_ALIGNMENT)
                for (size_t idx = 0; idx < natoms; idx++)
                {
                  const double ps_r = ps_r_ptr[idx];
                  const double ps_i = ps_i_ptr[idx];
                  vals_r[idx] += cG_j_r * ps_r - cG_j_i * ps_i;
                  vals_i[idx] += cG_j_i * ps_r + cG_j_r * ps_i;
                }
              }
            }
            else
            {
              for (size_t idx = 0; idx < natoms; idx++)
              {
                double* restrict vals_r           = vals_local[ip_idx][idx * 2].data();
                double* restrict vals_i           = vals_local[ip_idx][idx * 2 + 1].data();
                const double* restrict j_lm_G_ptr = j_lm_G.data();
                double cG_ps_r = cG_r * phase_shift_r[ig_local][idx] - cG_i * phase_shift_i[ig_local][idx];
                double cG_ps_i = cG_i * phase_shift_r[ig_local][idx] + cG_r * phase_shift_i[ig_local][idx];
#pragma omp simd aligned(vals_r, vals_i, j_lm_G_ptr : QMC_SIMD_ALIGNMENT)
                for (size_t lm = 0; lm < lm_tot; lm++)
                {
                  const double jlm = j_lm_G_ptr[lm];
                  vals_r[lm] += cG_ps_r * jlm;
                  vals_i[lm] += cG_ps_i * jlm;
                }
              }
            }
          }
        }
      }

#pragma omp for collapse(2)
      for (int ip = 0; ip < spline_npoints; ip++)
        for (size_t idx = 0; idx < natoms; idx++)
        {
          double* vals = all_vals[idx][ip];
          for (size_t tid = 0; tid < nt; tid++)
            for (size_t lm = 0; lm < lm_tot; lm++)
            {
              double vals_th_r, vals_th_i;
              const size_t ip_idx = tid * spline_npoints + ip;
              if (policy == 1)
              {
                vals_th_r = vals_local[ip_idx][lm * 2][idx];
                vals_th_i = vals_local[ip_idx][lm * 2 + 1][idx];
              }
              else
              {
                vals_th_r = vals_local[ip_idx][idx * 2][lm];
                vals_th_i = vals_local[ip_idx][idx * 2 + 1][lm];
              }
              const double real_tmp = 4.0 * M_PI * i_power[lm].real();
              const double imag_tmp = 4.0 * M_PI * i_power[lm].imag();
              vals[lm] += vals_th_r * real_tmp - vals_th_i * imag_tmp;
              vals[lm + lm_tot] += vals_th_i * real_tmp + vals_th_r * imag_tmp;
            }
        }
    }
    //app_log() << "Building band " << iorb << " at center " << center_idx << std::endl;

    for (size_t idx = 0; idx < natoms; idx++)
    {
      // reduce all_vals
      band_group_comm.reduce_in_place(all_vals[idx].data(), all_vals[idx].size());
      if (!band_group_comm.isGroupLeader())
        continue;
#pragma omp parallel for
      for (int lm = 0; lm < lm_tot; lm++)
      {
        auto& mycenter = centers[mygroup[idx]];
        aligned_vector<double> splineData_r(spline_npoints);
        UBspline_1d_d* atomic_spline_r = nullptr;
        for (size_t ip = 0; ip < spline_npoints; ip++)
          splineData_r[ip] = all_vals[idx][ip][lm];
        atomic_spline_r = einspline::create(atomic_spline_r, 0.0, spline_radius, spline_npoints, splineData_r.data(),
                                            ((lm == 0) || (lm > 3)));
        if (!bspline.isComplex())
        {
          mycenter.set_spline(atomic_spline_r, lm, iorb);
          einspline::destroy(atomic_spline_r);
        }
        else
        {
          aligned_vector<double> splineData_i(spline_npoints);
          UBspline_1d_d* atomic_spline_i = nullptr;
          for (size_t ip = 0; ip < spline_npoints; ip++)
            splineData_i[ip] = all_vals[idx][ip][lm + lm_tot];
          atomic_spline_i = einspline::create(atomic_spline_i, 0.0, spline_radius, spline_npoints, splineData_i.data(),
                                              ((lm == 0) || (lm > 3)));
          mycenter.set_spline(atomic_spline_r, lm, iorb * 2);
          mycenter.set_spline(atomic_spline_i, lm, iorb * 2 + 1);
          einspline::destroy(atomic_spline_r);
          einspline::destroy(atomic_spline_i);
        }
      }
    }
  }
}

template<typename SA>
void HybridRepSetReader<SA>::initialize_hybrid_pio_gather(const int spin,
                                                          const BandInfoGroup& bandgroup,
                                                          SA& bspline) const
{
  auto& mybuilder = spline_reader_.mybuilder;
  //distribute bands over processor groups
  int Nbands            = bandgroup.getNumDistinctOrbitals();
  const int Nprocs      = myComm->size();
  const int Nbandgroups = std::min(Nbands, Nprocs);
  Communicate band_group_comm(*myComm, Nbandgroups);
  std::vector<int> band_groups(Nbandgroups + 1, 0);
  FairDivideLow(Nbands, Nbandgroups, band_groups);
  int iorb_first = band_groups[band_group_comm.getGroupID()];
  int iorb_last  = band_groups[band_group_comm.getGroupID() + 1];

  app_log() << "Start transforming plane waves to 3D B-splines and atomic radial orbital 1D B-splines." << std::endl;
  OneSplineOrbData oneband(mybuilder->MeshSize, bspline.HalfG, bspline.isComplex());
  hdf_archive h5f(&band_group_comm, false);
  Vector<std::complex<double>> cG(mybuilder->Gvecs[0].size());
  const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
  if (band_group_comm.isGroupLeader())
    h5f.open(mybuilder->H5FileName, H5F_ACC_RDONLY);
  for (int iorb = iorb_first; iorb < iorb_last; iorb++)
  {
    if (band_group_comm.isGroupLeader())
    {
      const auto& cur_band = cur_bands[bspline.BandIndexMap[iorb]];
      const int ti         = cur_band.TwistIndex;
      spline_reader_.readOneOrbitalCoefs(psi_g_path(ti, spin, cur_band.BandIndex), h5f, cG);
      oneband.fft_spline(cG, mybuilder->Gvecs[0], mybuilder->primcell_kpoints[ti], rotate);
      bspline.set_spline(&oneband.get_spline_r(), &oneband.get_spline_i(), cur_band.TwistIndex, iorb, 0);
    }
    band_group_comm.bcast(cG);
    create_atomic_centers_Gspace(cG, band_group_comm, iorb, oneband.getRotatePhase(), bspline);
  }

  myComm->barrier();
  if (band_group_comm.isGroupLeader())
  {
    Timer now;
    bspline.gather_tables(band_group_comm.getGroupLeaderComm());
    app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
  }
}

#if defined(QMC_COMPLEX)
template class HybridRepSetReader<HybridRepCplx<SplineC2C<float>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2COMPTarget<float>>>;
#if !defined(QMC_MIXED_PRECISION)
template class HybridRepSetReader<HybridRepCplx<SplineC2C<double>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2COMPTarget<double>>>;
#endif
#else
template class HybridRepSetReader<HybridRepReal<SplineR2R<float>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2R<float>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2ROMPTarget<float>>>;
#if !defined(QMC_MIXED_PRECISION)
template class HybridRepSetReader<HybridRepReal<SplineR2R<double>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2R<double>>>;
template class HybridRepSetReader<HybridRepCplx<SplineC2ROMPTarget<double>>>;
#endif
#endif
} // namespace qmcplusplus
