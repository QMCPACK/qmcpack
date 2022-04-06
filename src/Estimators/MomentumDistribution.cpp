//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MomentumEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "MomentumDistribution.h"
#include "CPU/e2iphi.h"
#include "TrialWaveFunction.h"

#include <iostream>
#include <numeric>


namespace qmcplusplus
{
MomentumDistribution::MomentumDistribution(MomentumDistributionInput&& mdi,
                                           size_t np,
                                           const PosType& twist_in,
                                           const LatticeType& lattice,
                                           DataLocality dl)
    : OperatorEstBase(dl),
      input_(std::move(mdi)),
      twist(twist_in),
      Lattice(lattice),
      norm_nofK(1.0 / RealType(mdi.get_samples()))
{
  psi_ratios.resize(np);

  my_name_ = input_.get_name();

  //dims of a grid for generating k points (obtained below)
  int kgrid = 0;
  // minimal length as 2 x WS radius.
  RealType min_Length = Lattice.WignerSeitzRadius_G * 4.0 * M_PI;
  PosType vec_length;
  //length of reciprocal lattice vector
  for (int i = 0; i < OHMMS_DIM; i++)
    vec_length[i] = 2.0 * M_PI * std::sqrt(dot(Lattice.Gv[i], Lattice.Gv[i]));
  RealType kmax = input_.get_kmax();
  PosType kmaxs      = {input_.get_kmax0(), input_.get_kmax1(), input_.get_kmax2()};
  RealType sum_kmaxs = kmaxs[0] + kmaxs[1] + kmaxs[2];
  RealType sphere_kmax;
  bool sphere      = input_.get_kmax() > 0.0 ? true : false;
  bool directional = sum_kmaxs > 0.0 ? true : false;
  if (!sphere && !directional)
  {
    // default: kmax = 2 x k_F of polarized non-interacting electron system
    kmax   = 2.0 * std::pow(6.0 * M_PI * M_PI * np / Lattice.Volume, 1.0 / 3);
    sphere = true;
  }
  sphere_kmax = kmax;
  for (int i = 0; i < OHMMS_DIM; i++)
  {
    if (kmaxs[i] > kmax)
      kmax = kmaxs[i];
  }
  kgrid = int(kmax / min_Length) + 1;
  RealType kgrid_squared[OHMMS_DIM];
  for (int i = 0; i < OHMMS_DIM; i++)
    kgrid_squared[i] = kmaxs[i] * kmaxs[i] / vec_length[i] / vec_length[i];
  RealType kmax_squared = sphere_kmax * sphere_kmax;
  std::vector<int> kcount0;
  std::vector<int> kcount1;
  std::vector<int> kcount2;
  kcount0.resize(2 * kgrid + 1, 0);
  kcount1.resize(2 * kgrid + 1, 0);
  kcount2.resize(2 * kgrid + 1, 0);
  for (int i = -kgrid; i < (kgrid + 1); i++)
  {
    for (int j = -kgrid; j < (kgrid + 1); j++)
    {
      for (int k = -kgrid; k < (kgrid + 1); k++)
      {
        PosType ikpt, kpt;
        ikpt[0] = i + twist[0];
        ikpt[1] = j + twist[1];
        ikpt[2] = k + twist[2];
        //convert to Cartesian: note that 2Pi is multiplied
        kpt               = Lattice.k_cart(ikpt);
        bool not_recorded = true;
        // This collects the k-points within the parallelepiped (if enabled)
        if (directional && ikpt[0] * ikpt[0] <= kgrid_squared[0] && ikpt[1] * ikpt[1] <= kgrid_squared[1] &&
            ikpt[2] * ikpt[2] <= kgrid_squared[2])
        {
          kPoints.push_back(kpt);
          kcount0[kgrid + i] = 1;
          kcount1[kgrid + j] = 1;
          kcount2[kgrid + k] = 1;
          not_recorded       = false;
        }
        // This collects the k-points within a sphere (if enabled and the k-point has not been recorded yet)
        if (sphere && not_recorded &&
            kpt[0] * kpt[0] + kpt[1] * kpt[1] + kpt[2] * kpt[2] <=
                kmax_squared) //if (std::sqrt(kx*kx+ky*ky+kz*kz)<=sphere_kmax)
        {
          kPoints.push_back(kpt);
        }
      }
    }
  }
  app_log() << "\n  MomentumDistribution named " << my_name_ << "\n";
  if (sphere && !directional)
  {
    app_log() << "    Using all k-space points with (kx^2+ky^2+kz^2)^0.5 < " << sphere_kmax
              << " for Momentum Distribution." << std::endl;
    app_log() << "    Total number of k-points for Momentum Distribution is " << kPoints.size() << std::endl;
  }
  else if (directional && !sphere)
  {
    int sums[3];
    sums[0] = 0;
    sums[1] = 0;
    sums[2] = 0;
    for (int i = 0; i < 2 * kgrid + 1; i++)
    {
      sums[0] += kcount0[i];
      sums[1] += kcount1[i];
      sums[2] += kcount2[i];
    }
    app_log() << "    Using all k-space points within cut-offs " << input_.get_kmax0() << ", " << input_.get_kmax1() << ", " << input_.get_kmax2()
              << " for Momentum Distribution." << std::endl;
    app_log() << "    Total number of k-points for Momentum Distribution: " << kPoints.size() << std::endl;
    app_log() << "      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log() << "      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
    app_log() << "      Number of grid points in kmax2 direction: " << sums[2] << std::endl;
  }
  else
  {
    int sums[3];
    sums[0] = 0;
    sums[1] = 0;
    sums[2] = 0;
    for (int i = 0; i < 2 * kgrid + 1; i++)
    {
      sums[0] += kcount0[i];
      sums[1] += kcount1[i];
      sums[2] += kcount2[i];
    }
    app_log() << "    Using all k-space points with (kx^2+ky^2+kz^2)^0.5 < " << sphere_kmax << ", and" << std::endl;
    app_log() << "    within the cut-offs " << input_.get_kmax0() << ", " << input_.get_kmax1() << ", " << input_.get_kmax2() << " for Momentum Distribution."
              << std::endl;
    app_log() << "    Total number of k-points for Momentum Distribution is " << kPoints.size() << std::endl;
    app_log() << "    The number of k-points within the cut-off region: " << sums[0] * sums[1] * sums[2] << std::endl;
    app_log() << "      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log() << "      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
    app_log() << "      Number of grid points in kmax2 direction: " << sums[2] << std::endl;
  }
  app_log() << "    Number of samples: " << input_.get_samples() << std::endl;
  app_log() << "    My twist is: " << twist[0] << "  " << twist[1] << "  " << twist[2] << "\n\n";

  // resize arrays
  nofK.resize(kPoints.size());
  kdotp.resize(kPoints.size());
  auto samples = input_.get_samples();
  vPos.resize(samples);
  phases.resize(kPoints.size());
  phases_vPos.resize(samples);
  for (int im = 0; im < samples; im++)
    phases_vPos[im].resize(kPoints.size());
  psi_ratios_all.resize(samples, psi_ratios.size());

  // allocate data storage
  size_t data_size = nofK.size();
  data_.resize(data_size, 0.0);
}

MomentumDistribution::MomentumDistribution(const MomentumDistribution& md, DataLocality dl): MomentumDistribution(md) {
  data_locality_ = dl;
}
 
std::unique_ptr<OperatorEstBase> MomentumDistribution::spawnCrowdClone() const
{
  std::size_t data_size = data_.size();
  auto spawn_data_locality = data_locality_;

  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size = 0;
    throw std::runtime_error("There is no memory savings implementation for MomentumDistribution");
  }

  auto spawn = std::make_unique<MomentumDistribution>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size);
  return spawn;
}

//MomentumDistribution::MomentumDistribution(const MomentumDistribution& md)
//    : OperatorEstBase(md),
//      input_(std::move(md.input_))
//{
//  if (data_locality_ == DataLocality::crowd)
//  {
//    app_log()<<"MD::cons dl crowd\n";
//    size_t data_size = md.data_->size();
//    data_            = createLocalData(data_size, data_locality_);
//  }
//  else if (data_locality_ == DataLocality::rank)
//  {
//    app_log()<<"MD::cons dl rank\n";
//    assert(md.data_locality_ == DataLocality::rank);
//    size_t data_size  = 10; // jtk fix
//    data_locality_    = DataLocality::queue;
//    data_             = createLocalData(data_size, data_locality_);
//  }
//  else
//    app_log()<<"MD::cons dl other\n";
//}

void MomentumDistribution::startBlock(int steps)
{
  //if (data_locality_ == DataLocality::rank)
  //{
  //  app_log()<<"MD::startBlock dl rank\n";
  //  size_t data_size  = 10; // jtk fix
  //  data_->reserve(data_size);
  //  data_->resize(0);
  //}
  //else
  //  app_log()<<"MD::startBlock dl other\n";
}

/** Gets called every step and writes to thread local data.
 *
 */
void MomentumDistribution::accumulate(const RefVector<MCPWalker>& walkers,
                                      const RefVector<ParticleSet>& psets,
                                      const RefVector<TrialWaveFunction>& wfns,
                                      RandomGenerator& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker      = walkers[iw];
    ParticleSet& pset      = psets[iw];
    TrialWaveFunction& psi = wfns[iw];
    RealType weight        = walker.Weight;

    const int np = pset.getTotalNum();
    const int nk = kPoints.size();

    // accumulate weight
    //  (required by all estimators, otherwise inf results)
    walkers_weight_ += weight;

    auto samples = input_.get_samples();
    // compute phase factors
    for (int s = 0; s < samples; ++s)
    {
      PosType newpos;
      for (int i = 0; i < OHMMS_DIM; ++i)
        newpos[i] = rng();
      //make it cartesian
      vPos[s] = Lattice.toCart(newpos);
      pset.makeVirtualMoves(vPos[s]);
      psi.evaluateRatiosAlltoOne(pset, psi_ratios);
      for (int i = 0; i < np; ++i)
        psi_ratios_all[s][i] = psi_ratios[i];

      for (int ik = 0; ik < nk; ++ik)
        kdotp[ik] = -dot(kPoints[ik], vPos[s]);
      eval_e2iphi(nk, kdotp.data(), phases_vPos[s].data(0), phases_vPos[s].data(1));
    }

    // update n(k)
    std::fill_n(nofK.begin(), nk, RealType(0));
    for (int i = 0; i < np; ++i)
    {
      for (int ik = 0; ik < nk; ++ik)
        kdotp[ik] = dot(kPoints[ik], pset.R[i]);
      eval_e2iphi(nk, kdotp.data(), phases.data(0), phases.data(1));
      for (int s = 0; s < samples; ++s)
      {
        const ComplexType one_ratio(psi_ratios_all[s][i]);
        const RealType ratio_c                 = one_ratio.real();
        const RealType ratio_s                 = one_ratio.imag();
        const RealType* restrict phases_c      = phases.data(0);
        const RealType* restrict phases_s      = phases.data(1);
        const RealType* restrict phases_vPos_c = phases_vPos[s].data(0);
        const RealType* restrict phases_vPos_s = phases_vPos[s].data(1);
        RealType* restrict nofK_here           = nofK.data();
#pragma omp simd aligned(nofK_here, phases_c, phases_s, phases_vPos_c, phases_vPos_s : QMC_SIMD_ALIGNMENT)
        for (int ik = 0; ik < nk; ++ik)
          nofK_here[ik] += (phases_c[ik] * phases_vPos_c[ik] - phases_s[ik] * phases_vPos_s[ik]) * ratio_c -
              (phases_s[ik] * phases_vPos_c[ik] + phases_c[ik] * phases_vPos_s[ik]) * ratio_s;
      }
    }

    // accumulate data
    for (int ik = 0; ik < nofK.size(); ++ik)
      data_[ik] += weight * nofK[ik] * norm_nofK;

  }
}


void MomentumDistribution::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a MomentumDistribution with this DataLocality");
  }
}


void MomentumDistribution::registerOperatorEstimator(hid_t gid)
{
  //descriptor for the data, 1-D data
  std::vector<int> ng(1);
  //add nofk
  ng[0] = nofK.size();
  h5desc_.emplace_back(std::make_unique<ObservableHelper>("nofk"));
  auto& h5o = h5desc_.back();
  //h5o.set_dimensions(ng, my_index_);
  h5o->set_dimensions(ng, 0); // JTK: doesn't seem right
  h5o->open(gid);
  h5o->addProperty(const_cast<std::vector<PosType>&>(kPoints), "kpoints");
  h5o->addProperty(const_cast<std::vector<int>&>(kWeights), "kweights");
}


} // namespace qmcplusplus
