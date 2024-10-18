//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/{SkEstimator.cpp, SkAllEstimator.cpp}
//////////////////////////////////////////////////////////////////////////////////////

#include "StructureFactorEstimator.h"
#include "ParticleSet.h"
#include <LongRange/StructFact.h>

namespace qmcplusplus
{

StructureFactorEstimator::StructureFactorEstimator(StructureFactorInput& sfi,
                                                   ParticleSet& pset_ions,
                                                   ParticleSet& pset_elec,
                                                   DataLocality data_locality)
    : OperatorEstBase(data_locality), elns_(pset_elec), ions_(pset_ions)
{
  my_name_      = "StructureFactorEstimator";
  elec_species_ = elns_.getSpeciesSet().getTotalNum();
  ion_species_  = ions_.getSpeciesSet().getTotalNum();

  auto num_kpoints_ = pset_ions.getSimulationCell().getKLists().numk;
  kshell_offsets_   = pset_ions.getSimulationCell().getKLists().kshell;
  int max_kshell    = kshell_offsets_.size() - 1;

  rhok_tot_r_.resize(num_kpoints_);
  rhok_tot_i_.resize(num_kpoints_);

  // Legacy comment, but I don't trust it.
  //for values, we are including e-e structure factor, and e-Ion.  So a total of NumIonSpecies+1 structure factors.
  //+2 for the real and imaginary parts of rho_k^e
  //
  // skAll seems to be written for e-e sf + complex rho_k^e
  data_.resize(3 * num_kpoints_);
  kmags_.resize(max_kshell);
  one_over_degeneracy_kshell_.resize(max_kshell);
  for (int ks = 0; ks < max_kshell; ks++)
  {
    kmags_[ks]                      = std::sqrt(pset_elec.getSimulationCell().getKLists().ksq[kshell_offsets_[ks]]);
    one_over_degeneracy_kshell_[ks] = 1.0 / static_cast<Real>(kshell_offsets_[ks + 1] - kshell_offsets_[ks]);
  };
}

void StructureFactorEstimator::accumulate(const RefVector<MCPWalker>& walkers,
                                          const RefVector<ParticleSet>& psets,
                                          const RefVector<TrialWaveFunction>& wfns,
                                          const RefVector<QMCHamiltonian>& hams,
                                          RandomBase<FullPrecReal>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    Real weight = walkers[iw].get().Weight;

    //sum over species
    std::copy(psets[iw].get().getSK().rhok_r[0], psets[iw].get().getSK().rhok_r[0] + num_kpoints_, rhok_tot_r_.begin());
    std::copy(psets[iw].get().getSK().rhok_i[0], psets[iw].get().getSK().rhok_i[0] + num_kpoints_, rhok_tot_i_.begin());
    for (int i = 1; i < elec_species_; ++i)
      accumulate_elements(psets[iw].get().getSK().rhok_r[i], psets[iw].get().getSK().rhok_r[i] + num_kpoints_,
                          rhok_tot_r_.begin());
    for (int i = 1; i < elec_species_; ++i)
      accumulate_elements(psets[iw].get().getSK().rhok_i[i], psets[iw].get().getSK().rhok_i[i] + num_kpoints_,
                          rhok_tot_i_.begin());

    for (int k = 0; k < num_kpoints_; k++)
    {
      sfk_e_e_[k] += weight * (rhok_tot_r_[k] * rhok_tot_r_[k] + rhok_tot_i_[k] * rhok_tot_i_[k]);
      rhok_e_[k] += weight * std::complex<Real>{rhok_tot_r_[k], rhok_tot_i_[k]};
    }
  }
}

void StructureFactorEstimator::registerOperatorEstimator(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};
  hdf_path path_variables = hdf_name / std::string_view("kpoints");
  file.push(path_variables, true);
  file.write(const_cast<std::vector<KPt>&>(ions_.getSimulationCell().getKLists().kpts_cart), "value");
  file.pop();
}

void StructureFactorEstimator::write(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};
  file.push(hdf_name);
  // this is call rhok_e_e in the output of the legacy, but that is just wrong it is |rhok_e_e_|^2
  file.write(sfk_e_e_, "sfk_e_e");
  file.write(rhok_e_, "rhok_e_");
}

void StructureFactorEstimator::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  int num_crowds = type_erased_operator_estimators.size();
  for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
  {
    StructureFactorEstimator& crowd_sfe = dynamic_cast<StructureFactorEstimator&>(crowd_oeb);
    this->sfk_e_e_ += crowd_sfe.sfk_e_e_;
    this->rhok_e_ += crowd_sfe.rhok_e_;
  }
  OperatorEstBase::collect(type_erased_operator_estimators);
}

void StructureFactorEstimator::startBlock(int steps) {}

UPtr<OperatorEstBase> StructureFactorEstimator::spawnCrowdClone() const {}

} // namespace qmcplusplus
