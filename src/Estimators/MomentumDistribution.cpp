//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MomentumEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "MomentumDistribution.h"

#include <iostream>
#include <numeric>

namespace qmcplusplus
{
MomentumDistribution::MomentumDistribution(MomentumDistributionInput&& mdi, size_t np, const PosType& twist_in, const LatticeType& lattice, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(mdi)), twist(twist_in), Lattice(lattice)
{
  data_locality_ = dl;
  psi_ratios.resize(np);

  myName = input_.get<std::string>("name");
  M      = input_.get<int>("samples");

  norm_nofK = 1.0 / RealType(M);
  
  //maximum k-value in the k-grid in cartesian coordinates
  auto kmax = input_.get<RealType>("kmax");
  //maximum k-values in the k-grid along the reciprocal cell axis
  auto kmax0 = input_.get<RealType>("kmax0");
  auto kmax1 = input_.get<RealType>("kmax1");
  auto kmax2 = input_.get<RealType>("kmax2");

  //dims of a grid for generating k points (obtained below)
  size_t kgrid = 0;
  // minimal length as 2 x WS radius.
  RealType min_Length = Lattice.WignerSeitzRadius_G * 4.0 * M_PI;
  PosType vec_length;
  //length of reciprocal lattice vector
  for (int i = 0; i < OHMMS_DIM; i++)
    vec_length[i] = 2.0 * M_PI * std::sqrt(dot(Lattice.Gv[i], Lattice.Gv[i]));
  PosType kmaxs = {kmax0,kmax1,kmax2};
  RealType sum_kmaxs = kmaxs[0] + kmaxs[1] + kmaxs[2];
  RealType sphere_kmax;
  bool sphere      = kmax > 0.0 ? true : false;
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
        ikpt[0] = i - twist[0];
        ikpt[1] = j - twist[1];
        ikpt[2] = k - twist[2];
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
  app_log() <<"\n  MomentumDistribution named "<<myName<<"\n";
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
    app_log() << "    Using all k-space points within cut-offs " << kmax0 << ", " << kmax1 << ", " << kmax2
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
    app_log() << "    within the cut-offs " << kmax0 << ", " << kmax1 << ", " << kmax2 << " for Momentum Distribution."
              << std::endl;
    app_log() << "    Total number of k-points for Momentum Distribution is " << kPoints.size() << std::endl;
    app_log() << "    The number of k-points within the cut-off region: " << sums[0] * sums[1] * sums[2] << std::endl;
    app_log() << "      Number of grid points in kmax0 direction: " << sums[0] << std::endl;
    app_log() << "      Number of grid points in kmax1 direction: " << sums[1] << std::endl;
    app_log() << "      Number of grid points in kmax2 direction: " << sums[2] << std::endl;
  }
  app_log() << "    Number of samples: " << M << std::endl;
  app_log() << "    My twist is: " << twist[0] << "  " << twist[1] << "  " << twist[2] << "\n\n";

  // resize arrays
  nofK.resize(kPoints.size());
  kdotp.resize(kPoints.size());
  vPos.resize(M);
  phases.resize(kPoints.size());
  phases_vPos.resize(M);
  for (int im = 0; im < M; im++)
    phases_vPos[im].resize(kPoints.size());
  psi_ratios_all.resize(M, psi_ratios.size());

  // allocate data storage
  size_t data_size  = nofK.size();
  data_             = createLocalData(data_size, data_locality_);

}


MomentumDistribution* MomentumDistribution::clone()
{
  auto* md = new MomentumDistribution(*this);
  if (md->data_locality_ == DataLocality::crowd)
  {
    app_log()<<"MD::clone dl crowd\n";
    size_t data_size = data_->size();
    md->data_        = createLocalData(data_size, data_locality_);
  }
  else if (md->data_locality_ == DataLocality::rank)
  {
    app_log()<<"MD::clone dl rank\n";
    assert(data_locality_ == DataLocality::rank);
    size_t data_size  = 10; // jtk fix
    md->data_locality_    = DataLocality::queue;
    md->data_             = createLocalData(data_size, data_locality_);
  }
  else
    app_log()<<"MD::clone dl other\n";
  
  return md;
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
  if (data_locality_ == DataLocality::rank)
  {
    app_log()<<"MD::startBlock dl rank\n";
    size_t data_size  = 10; // jtk fix
    data_->reserve(data_size);
    data_->resize(0);
  }
  else
    app_log()<<"MD::startBlock dl other\n";

}

/** Gets called every step and writes to thread local data.
 *
 */
void MomentumDistribution::accumulate(const RefVector<MCPWalker>& walkers, const RefVector<ParticleSet>& psets)
{
};


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
  //h5o.set_dimensions(ng, myIndex);
  h5o->set_dimensions(ng, 0); // JTK: doesn't seem right
  h5o->open(gid);
  h5o->addProperty(const_cast<std::vector<PosType>&>(kPoints), "kpoints");
  h5o->addProperty(const_cast<std::vector<int>&>(kWeights), "kweights");    
}


} // namespace qmcplusplus
