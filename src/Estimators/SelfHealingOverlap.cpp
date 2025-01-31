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
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

#include <iostream>
#include <numeric>


namespace qmcplusplus
{
SelfHealingOverlap::SelfHealingOverlap(SelfHealingOverlapInput&& inp_, const TrialWaveFunction& wfn, DataLocality dl)
    : OperatorEstBase(dl),
      input_(std::move(inp_)),
      wf_type(no_wf),
      use_param_deriv(input_.input_section_.get<bool>("param_deriv"))
{
  auto msd_refvec = wfn.findMSD();
  auto sd_refvec  = wfn.findSD();

  auto nsd  = sd_refvec.size();
  auto nmsd = msd_refvec.size();

  size_t nparams;
  if (nmsd == 1 && nsd == 0)
  { // multi-slater-det wavefunction
    wf_type                              = msd_wf;
    const MultiSlaterDetTableMethod& msd = msd_refvec[0];
    if (!use_param_deriv)
      nparams = msd.getLinearExpansionCoefs().size();
    else
    {
      throw std::runtime_error("SelfHealingOverlap: use_param_deriv implementation incomplete, needs access to param "
                               "count from wavefunction component myVars");
    }
    if (nparams == 0)
      throw std::runtime_error("SelfHealingOverlap: multidet wavefunction has no parameters.");
  }
  else if (nmsd == 0 && nsd == 1)
  { // slater-det wavefunction
    throw std::runtime_error("SelfHealingOverlap: slaterdet wavefunction implementation incomplete");
  }
  else
  {
    throw std::runtime_error(
        "SelfHealingOverlap requires a single slater or multi-slater determinant component in the trial wavefunction.");
  }

#ifndef QMC_COMPLEX
  const size_t data_size = nparams;
#else
  const size_t data_size = 2 * nparams;
#endif
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

    // find jastrow wavefunction components
    std::vector<WaveFunctionComponent*> wcs_jastrow;
    for (auto& wc : wcs)
      if (!wc->isFermionic())
        wcs_jastrow.push_back(wc.get());

    if (wf_type == msd_wf)
    {
      auto msd_refvec                = psi.findMSD();
      MultiSlaterDetTableMethod& msd = msd_refvec[0];
      // collect parameter derivatives: (dpsi/dc_i)/psi
      if (!use_param_deriv)
        msd.calcIndividualDetRatios(det_ratios);
      else
      {
        throw std::runtime_error("SelfHealingOverlap: use_param_deriv implementation incomplete, needs call to "
                                 "msd.evaluateDerivatives with correct myVars");
      }
    }
    else if (wf_type == sd_rot_wf)
    {
      throw std::runtime_error("SelfHealingOverlap: slaterdet wavefunction implementation incomplete");
      auto sd_refvec = psi.findSD();
    }
    else
      throw std::runtime_error("SelfHealingOverlap: impossible branch reached, contact the developers");

    // collect jastrow prefactor
    WaveFunctionComponent::LogValue Jval = 0.0;
    for (auto& wc : wcs_jastrow)
      Jval += wc->get_log_value();
    RealType Jprefactor = std::exp(std::real(-2. * Jval));

    // accumulate weight (required by all estimators, otherwise inf results)
    walkers_weight_ += weight;

    // accumulate data
    assert(det_ratios.size() == data_.size());
    for (int ic = 0; ic < det_ratios.size(); ++ic)
    {
#ifndef QMC_COMPLEX
      data_[ic] += weight * Jprefactor * det_ratios[ic];
#else
      auto value = weight * Jprefactor * std::conj(det_ratios[ic]);
      data_[2 * ic] += std::real(value);
      data_[2 * ic + 1] += std::imag(value);
#endif
    }
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
