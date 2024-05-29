//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "SelfHealingOverlapLegacy.h"
#include "TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h" 
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
SelfHealingOverlapLegacy::SelfHealingOverlapLegacy(TrialWaveFunction& wfn)
  : psi_ref(wfn)
{
  name_ = "SelfHealingOverlap";
  update_mode_.set(COLLECTABLE, 1);

  auto msd_refvec = wfn.findMSD();
  if (msd_refvec.size() != 1)
    throw std::runtime_error(
        "SelfHealingOverlap requires one and only one multi slater determinant component in the trial wavefunction.");

  const MultiSlaterDetTableMethod& msd = msd_refvec[0];

  ncoef = msd.getLinearExpansionCoefs().size();
}


std::unique_ptr<OperatorBase> SelfHealingOverlapLegacy::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return std::make_unique<SelfHealingOverlapLegacy>(psi);
}


bool SelfHealingOverlapLegacy::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(name_, "name");
  attrib.put(cur);
  return true;
}



void SelfHealingOverlapLegacy::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<RealType> tmp(ncoef);
  collectables.add(tmp.begin(), tmp.end());
}


void SelfHealingOverlapLegacy::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  std::vector<int> ng(1);
  ng[0] = ncoef;
  h5desc.push_back({{"sh_coeff"}});
  auto& h5o = h5desc.back();
  h5o.set_dimensions(ng, my_index_);
}


SelfHealingOverlapLegacy::Return_t SelfHealingOverlapLegacy::evaluate(ParticleSet& P)
{
  RealType weight = t_walker_->Weight;
  int offset = my_index_;
  auto& wcs  = psi_ref.getOrbitals();

  // separate jastrow and fermi wavefunction components
  std::vector<WaveFunctionComponent*> wcs_jastrow;
  std::vector<WaveFunctionComponent*> wcs_fermi;
  for (auto& wc : wcs)
    if (wc->isFermionic())
      wcs_fermi.push_back(wc.get());
    else
      wcs_jastrow.push_back(wc.get());
  auto msd_refvec = psi_ref.findMSD();
  MultiSlaterDetTableMethod& msd = msd_refvec[0];

  // fermionic must have only one component, and must be multideterminant
  assert(wcs_fermi.size() == 1);
  WaveFunctionComponent& wf = *wcs_fermi[0];
  if (!wf.isMultiDet())
    throw std::runtime_error("SelfHealingOverlap estimator requires use of multideterminant wavefunction");

  // collect parameter derivatives: (dpsi/dc_i)/psi
  msd.calcIndividualDetRatios(det_ratios);

  // collect jastrow prefactor
  WaveFunctionComponent::LogValue Jval = 0.0;
  for (auto& wc : wcs_jastrow)
    Jval += wc->get_log_value();
  auto Jprefactor = std::real(std::exp(-2. * Jval));

  // accumulate data
  assert(det_ratios.size() == ncoef);
  for (int ic = 0; ic < det_ratios.size(); ++ic)
    P.Collectables[offset+ic] += weight * Jprefactor * real(det_ratios[ic]); // only real supported for now

  return 0.0;
}


} // namespace qmcplusplus
