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

namespace qmcplusplus
{
using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
using GridType       = LRCoulombSingleton::GridType;
using RadFunctorType = LRCoulombSingleton::RadFunctorType;

DensityEstimator::DensityEstimator(ParticleSet& elns)
{
  update_mode_.set(COLLECTABLE, 1);
  Periodic = (elns.getLattice().SuperCellEnum != SUPERCELL_OPEN);
  for (int dim = 0; dim < OHMMS_DIM; ++dim)
  {
    density_max[dim] = elns.getLattice().Length[dim];
    ScaleFactor[dim] = 1.0 / elns.getLattice().Length[dim];
  }
}

void DensityEstimator::resetTargetParticleSet(ParticleSet& P) {}

DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
{
  RealType wgt = t_walker_->Weight;
  if (Periodic)
  {
    for (int iat = 0; iat < P.getTotalNum(); ++iat)
    {
      PosType ru;
      ru    = P.getLattice().toUnit(P.R[iat]);
      int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
      int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
      int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
      P.Collectables[getGridIndex(i, j, k)] += wgt; //1.0;
      //	P.Collectables[getGridIndexPotential(i,j,k)]-=1.0;
    }
  }
  else
  {
    for (int iat = 0; iat < P.getTotalNum(); ++iat)
    {
      PosType ru;
      for (int dim = 0; dim < OHMMS_DIM; dim++)
      {
        //ru[dim]=(P.R[iat][dim]-density_min[dim])/(density_max[dim]-density_min[dim]);
        ru[dim] = (P.R[iat][dim] - density_min[dim]) * ScaleFactor[dim];
      }
      if (ru[0] > 0.0 && ru[1] > 0.0 && ru[2] > 0.0 && ru[0] < 1.0 && ru[1] < 1.0 && ru[2] < 1.0)
      {
        int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
        int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
        int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
        P.Collectables[getGridIndex(i, j, k)] += wgt; //1.0;
        //	  P.Collectables[getGridIndexPotential(i,j,k)]-=1.0;
      }
    }
  }
  return 0.0;
}

void DensityEstimator::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  int nw = W.WalkerList.size();
  int N  = W.getTotalNum();
  if (Periodic)
  {
    for (int iw = 0; iw < nw; iw++)
    {
      Walker_t& w     = *W.WalkerList[iw];
      RealType weight = w.Weight / nw;
      for (int iat = 0; iat < N; iat++)
      {
        PosType ru;
        ru = W.getLattice().toUnit(w.R[iat]);
        // for (int dim=0; dim<OHMMS_DIM; dim++)
        // {
        //   ru[dim]=(w.R[iat][dim]-density_min[dim])*ScaleFactor[dim];
        // }
        int i = static_cast<int>(DeltaInv[0] * (ru[0] - std::floor(ru[0])));
        int j = static_cast<int>(DeltaInv[1] * (ru[1] - std::floor(ru[1])));
        int k = static_cast<int>(DeltaInv[2] * (ru[2] - std::floor(ru[2])));
        W.Collectables[getGridIndex(i, j, k)] += weight;
      }
    }
  }
}

void DensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  //current index
  my_index_ = collectables.current();
  std::vector<RealType> tmp(NumGrids[OHMMS_DIM]);
  collectables.add(tmp.begin(), tmp.end());
  //potentialIndex=collectables.current();
  //vector<RealType> tmp2(NumGrids[OHMMS_DIM]);
  //collectables.add(tmp2.begin(),tmp2.end());
}

void DensityEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  int loc = h5desc.size();
  std::vector<int> ng(OHMMS_DIM);
  for (int i = 0; i < OHMMS_DIM; ++i)
    ng[i] = NumGrids[i];
  h5desc.emplace_back(name_);
  auto& h5o = h5desc.back();
  h5o.set_dimensions(ng, my_index_);
  h5o.open(gid);
}

void DensityEstimator::setObservables(PropertySetType& plist)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex);
}

void DensityEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex+offset);
}

/** check xml elements
 *
 * <estimator name="density" debug="no" delta="0.1 0.1 0.1"/>
 */
bool DensityEstimator::put(xmlNodePtr cur)
{
  Delta = 0.1;
  std::vector<double> delta;
  std::string debug("no");
  std::string potential("no");
  OhmmsAttributeSet attrib;
  attrib.add(debug, "debug");
  attrib.add(potential, "potential");
  attrib.add(density_min[0], "x_min");
  attrib.add(density_min[1], "y_min");
  attrib.add(density_min[2], "z_min");
  attrib.add(density_max[0], "x_max");
  attrib.add(density_max[1], "y_max");
  attrib.add(density_max[2], "z_max");
  attrib.add(Delta, "delta");
  attrib.put(cur);
  if (!Periodic)
  {
    for (int dim = 0; dim < OHMMS_DIM; ++dim)
      ScaleFactor[dim] = 1.0 / (density_max[dim] - density_min[dim]);
  }
  resize();
  return true;
}

bool DensityEstimator::get(std::ostream& os) const
{
  os << name_ << " bin =" << Delta << " bohrs " << std::endl;
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
    DeltaInv[i] = 1.0 / Delta[i];
    NumGrids[i] = static_cast<int>(DeltaInv[i]);
    if (NumGrids[i] < 2)
    {
      APP_ABORT("DensityEstimator::resize invalid bin size");
    }
  }
  app_log() << " DensityEstimator bin_size= " << NumGrids << " delta = " << Delta << std::endl;
  NumGrids[OHMMS_DIM] = NumGrids[0] * NumGrids[1] * NumGrids[2];
}

} // namespace qmcplusplus
