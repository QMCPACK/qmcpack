//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////


#include "DensityEstimator.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/LRCoulombSingleton.h"
#include "Particle/DistanceTable.h"
#include "Particle/MCWalkerConfiguration.h"

#include <stdexcept> // std::invalid_argument

namespace qmcplusplus
{
using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
using GridType       = LRCoulombSingleton::GridType;
using RadFunctorType = LRCoulombSingleton::RadFunctorType;


DensityEstimator::DensityEstimator(ParticleSet& elns)
{
  update_mode_.set(COLLECTABLE, 1);
  periodic_ = (elns.getLattice().SuperCellEnum != SUPERCELL_OPEN);
  for (int dim = 0; dim < OHMMS_DIM; ++dim)
  {
    density_max_[dim]  = elns.getLattice().Length[dim];
    scale_factor_[dim] = 1.0 / elns.getLattice().Length[dim];
  }
}

std::string DensityEstimator::getClassName() const { return "DensityEstimator"; }

void DensityEstimator::resetTargetParticleSet(ParticleSet& P) {}

DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
{
  if (t_walker_ == nullptr)
  {
    throw std::invalid_argument("calling DensityEstimator::evaluate needs the weight of the walkers");
  }

  const RealType wgt = t_walker_->Weight;

  if (periodic_)
  {
    for (int iat = 0; iat < P.getTotalNum(); ++iat)
    {
      PosType ru;
      ru    = P.getLattice().toUnit(P.R[iat]);
      int i = static_cast<int>(delta_inv_[0] * (ru[0] - std::floor(ru[0])));
      int j = static_cast<int>(delta_inv_[1] * (ru[1] - std::floor(ru[1])));
      int k = static_cast<int>(delta_inv_[2] * (ru[2] - std::floor(ru[2])));
      P.Collectables[getGridIndex(i, j, k)] += wgt; //1.0;
    }
  }
  else
  {
    for (int iat = 0; iat < P.getTotalNum(); ++iat)
    {
      PosType ru;
      for (int dim = 0; dim < OHMMS_DIM; dim++)
      {
        ru[dim] = (P.R[iat][dim] - density_min_[dim]) * scale_factor_[dim];
      }
      if (ru[0] > 0.0 && ru[1] > 0.0 && ru[2] > 0.0 && ru[0] < 1.0 && ru[1] < 1.0 && ru[2] < 1.0)
      {
        int i = static_cast<int>(delta_inv_[0] * (ru[0] - std::floor(ru[0])));
        int j = static_cast<int>(delta_inv_[1] * (ru[1] - std::floor(ru[1])));
        int k = static_cast<int>(delta_inv_[2] * (ru[2] - std::floor(ru[2])));
        P.Collectables[getGridIndex(i, j, k)] += wgt; //1.0;
      }
    }
  }
  return 0.0;
}

void DensityEstimator::addObservables(PropertySetType& plist) {}

void DensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  //current index
  my_index_ = collectables.current();
  std::vector<RealType> tmp(num_grids_[OHMMS_DIM]);
  collectables.add(tmp.begin(), tmp.end());
}

void DensityEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  int loc = h5desc.size();
  std::vector<int> ng(OHMMS_DIM);
  for (int i = 0; i < OHMMS_DIM; ++i)
    ng[i] = num_grids_[i];
  h5desc.emplace_back(hdf_path{name_});
  auto& h5o = h5desc.back();
  h5o.set_dimensions(ng, my_index_);
}

void DensityEstimator::setObservables(PropertySetType& plist) {}

void DensityEstimator::setParticlePropertyList(PropertySetType& plist, int offset) {}

/** check xml elements
 *
 * <estimator name="density" debug="no" delta="0.1 0.1 0.1"/>
 */
bool DensityEstimator::put(xmlNodePtr cur)
{
  delta_ = 0.1;
  std::vector<double> delta;
  std::string debug("no");
  std::string potential("no");
  OhmmsAttributeSet attrib;
  attrib.add(debug, "debug");
  attrib.add(potential, "potential");
  attrib.add(density_min_[0], "x_min");
  attrib.add(density_min_[1], "y_min");
  attrib.add(density_min_[2], "z_min");
  attrib.add(density_max_[0], "x_max");
  attrib.add(density_max_[1], "y_max");
  attrib.add(density_max_[2], "z_max");
  attrib.add(delta_, "delta");
  attrib.put(cur);
  if (!periodic_)
  {
    for (int dim = 0; dim < OHMMS_DIM; ++dim)
      scale_factor_[dim] = 1.0 / (density_max_[dim] - density_min_[dim]);
  }
  resize();
  return true;
}

bool DensityEstimator::get(std::ostream& os) const
{
  os << name_ << " bin =" << delta_ << " bohrs " << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> DensityEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //default constructor is sufficient
  return std::make_unique<DensityEstimator>(*this);
}

void DensityEstimator::resize()
{
  for (int i = 0; i < OHMMS_DIM; ++i)
  {
    delta_inv_[i] = 1.0 / delta_[i];
    num_grids_[i] = static_cast<int>(delta_inv_[i]);
    if (num_grids_[i] < 2)
    {
      APP_ABORT("DensityEstimator::resize invalid bin size");
    }
  }
  app_log() << " DensityEstimator bin_size= " << num_grids_ << " delta = " << delta_ << std::endl;
  num_grids_[OHMMS_DIM] = num_grids_[0] * num_grids_[1] * num_grids_[2];
}

int DensityEstimator::getGridIndex(int i, int j, int k) const noexcept
{
  return my_index_ + k + num_grids_[2] * (j + num_grids_[1] * i);
}

} // namespace qmcplusplus
