//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "SelfHealingOverlap.h"
#include "TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"

#include <iostream>
#include <numeric>


namespace qmcplusplus
{
SelfHealingOverlap::SelfHealingOverlap(SelfHealingOverlapInput&& inp_, const TrialWaveFunction& wfn, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(inp_))
{
  //my_name_ = input_.get_name();

  auto& inp = this->input_.input_section_;

  auto msd_refvec = wfn.findMSD();
  if (msd_refvec.size() != 1)
    throw std::runtime_error(
        "SelfHealingOverlap requires one and only one multi slater determinant component in the trial wavefunction.");

  const MultiSlaterDetTableMethod& msd = msd_refvec[0];
  const size_t data_size = msd.getLinearExpansionCoefs().size();
  data_.resize(data_size, 0.0);
}


SelfHealingOverlap::SelfHealingOverlap(const SelfHealingOverlap& sh, DataLocality dl) : SelfHealingOverlap(sh)
{
  data_locality_ = dl;
}

std::unique_ptr<OperatorEstBase> SelfHealingOverlap::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for SelfHealingOverlap");
  }

  auto spawn = std::make_unique<SelfHealingOverlap>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size);
  return spawn;
}

void SelfHealingOverlap::startBlock(int steps) {}

/** Gets called every step and writes to thread local data.
 *
 */
void SelfHealingOverlap::accumulate(const RefVector<MCPWalker>& walkers,
                                    const RefVector<ParticleSet>& psets,
                                    const RefVector<TrialWaveFunction>& wfns,
                                    const RefVector<QMCHamiltonian>& hams,
                                    RandomBase<FullPrecRealType>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker      = walkers[iw];
    ParticleSet& pset      = psets[iw];
    TrialWaveFunction& psi = wfns[iw];
    RealType weight        = walker.Weight;
    auto& wcs              = psi.getOrbitals();

    // separate jastrow and fermi wavefunction components
    std::vector<WaveFunctionComponent*> wcs_jastrow;
    std::vector<WaveFunctionComponent*> wcs_fermi;
    for (auto& wc : wcs)
      if (wc->isFermionic())
        wcs_fermi.push_back(wc.get());
      else
        wcs_jastrow.push_back(wc.get());

    // fermionic must have only one component, and must be multideterminant
    assert(wcs_fermi.size() == 1);
    WaveFunctionComponent& wf = *wcs_fermi[0];
    if (!wf.isMultiDet())
      throw std::runtime_error("SelfHealingOverlap estimator requires use of multideterminant wavefunction");
    auto msd_refvec = psi.findMSD();
    MultiSlaterDetTableMethod& msd = msd_refvec[0];

    // collect parameter derivatives: (dpsi/dc_i)/psi
    msd.calcIndividualDetRatios(det_ratios);

    // collect jastrow prefactor
    WaveFunctionComponent::LogValue Jval = 0.0;
    for (auto& wc : wcs_jastrow)
      Jval += wc->get_log_value();
    auto Jprefactor = std::real(std::exp(-2. * Jval));

    // accumulate weight (required by all estimators, otherwise inf results)
    walkers_weight_ += weight;

    // accumulate data
    assert(det_ratios.size() == data_.size());
    for (int ic = 0; ic < det_ratios.size(); ++ic)
      data_[ic] += weight * Jprefactor * real(det_ratios[ic]); // only real supported for now
  }
}


void SelfHealingOverlap::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a SelfHealingOverlap with this DataLocality");
  }
}


void SelfHealingOverlap::registerOperatorEstimator(hdf_archive& file)
{
  using namespace std::string_literals;
  ////descriptor for the data, 1-D data
  std::vector<int> ng(1);
  ng[0] = data_.size();
  h5desc_.push_back({{"sh_coeff"}});
  auto& h5o = h5desc_.back();
  h5o.set_dimensions(ng, 0); // JTK: doesn't seem right
}


} // namespace qmcplusplus
