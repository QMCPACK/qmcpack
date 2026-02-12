/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/EnergyDensityEstimator.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <OhmmsPETE/TinyVector.h>
#include "EnergyDensityEstimator.h"

namespace qmcplusplus
{

using Real = NEEnergyDensityEstimator::Real;

struct PosCharge
{
  NEEnergyDensityEstimator::ParticlePos r_ptcls;
  std::vector<NEEnergyDensityEstimator::Real> z_ptcls;
};

auto NEEnergyDensityEstimator::extractIonPositionsAndCharge(const ParticleSet& pset)
{
  const SpeciesSet& species(pset.getSpeciesSet());
  int charge_index = species.findAttribute("charge");
  int nspecies     = species.TotalNum;
  int nps          = pset.getTotalNum();
  std::vector<Real> z_spec;
  z_spec.resize(nspecies);
  std::vector<Real> z_ptcls;
  z_ptcls.resize(nps);
  for (int spec = 0; spec < nspecies; spec++)
    z_spec[spec] = species(charge_index, spec);
  for (int i = 0; i < nps; i++)
    z_ptcls[i] = z_spec[pset.GroupID[i]];
  ParticlePos r_ptcls;
  r_ptcls.resize(pset.R.size());
  for (int i = 0; i < pset.R.size(); i++)
    r_ptcls[i] = pset.R[i];
  if (pset.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    pset.applyMinimumImage(r_ptcls);
  return PosCharge{r_ptcls, z_ptcls};
}

NEEnergyDensityEstimator::NEEnergyDensityEstimator(const EnergyDensityInput& input,
                                                   const PSPool& pset_pool,
                                                   DataLocality data_locality)
    : OperatorEstBase(data_locality, input.get_name(), input.get_type()),
      input_(input),
      pset_dynamic_(getParticleSet(pset_pool, input.get_dynamic()))
{
  requires_listener_ = true;
  my_name_           = input.get_name();
  n_particles_       = pset_dynamic_.getTotalNum();
  if (!(input_.get_static().empty()))
  {
    pset_static_.emplace(getParticleSet(pset_pool, input.get_static()));
    if (!input.get_ion_points())
    {
      n_ions_ = pset_static_->getTotalNum();
      n_particles_ += n_ions_;
    }
    else
    {
      ed_ion_values_.resize(n_ions_, N_EDVALS);
      r_ion_work_.resize(n_ions_, OHMMS_DIM);
    }
  }

  constructToReferencePoints(pset_dynamic_, pset_static_);

  bool periodic = pset_dynamic_.getLattice().SuperCellEnum != SUPERCELL_OPEN;

  spacegrid_inputs_ = input.get_space_grid_inputs();
  spacegrids_.reserve(spacegrid_inputs_.size());
  for (int ig = 0; ig < spacegrid_inputs_.size(); ++ig)
  {
    spacegrids_.emplace_back(
        std::make_unique<NESpaceGrid<Real>>(spacegrid_inputs_[ig], ref_points_->get_points(), N_EDVALS, periodic));
  }
}

NEEnergyDensityEstimator::NEEnergyDensityEstimator(const NEEnergyDensityEstimator& ede, const DataLocality dl)
    : OperatorEstBase(dl, ede.input_.get_name(), ede.input_.get_type()),
      input_(ede.input_),
      pset_dynamic_(ede.pset_dynamic_),
      pset_static_(ede.pset_static_),
      n_particles_(ede.n_particles_),
      n_ions_(ede.n_ions_),
      spacegrid_inputs_(ede.spacegrid_inputs_)
{
  requires_listener_ = true;
  data_locality_     = dl;

  constructToReferencePoints(pset_dynamic_, pset_static_);

  for (const auto& space_grid : ede.spacegrids_)
    spacegrids_.emplace_back(std::make_unique<NESpaceGrid<Real>>(*space_grid));
}

void NEEnergyDensityEstimator::constructToReferencePoints(ParticleSet& pset_dynamic,
                                                          const std::optional<ParticleSet>& pset_static)
{
  // Bringing a bunch of adhoc setup from legacy.
  // removed redundant turnOnPerParticleSK this is handled via the EstimatorManagerNew if listeners are detected through
  // CoulombPBCAA{AB} which are the actual operators that need it.
  // Although it is possible that our copies of the particle sets will need per particle structure factors I don't think
  // they do.

  RefVector<ParticleSet> pset_refs;

  if (pset_static_)
  {
    dtable_index_ = pset_dynamic_.addTable(pset_static_.value());
    pset_refs.emplace_back(pset_static_.value());
  }

  r_work_.resize(n_particles_);
  ed_values_.resize(n_particles_, N_EDVALS);

  // right now this is only the case when pset_static_ && input.get_ion_points_ are true
  if (n_ions_ > 0)
  {
    ed_ion_values_.resize(n_ions_, N_EDVALS);
    r_ion_work_.resize(n_ions_, OHMMS_DIM);
    for (int i = 0; i < n_ions_; ++i)
      for (int d = 0; d < OHMMS_DIM; ++d)
        r_ion_work_[i][d] = pset_static_->R[i][d];
  }
  data_.resize(N_EDVALS);
  particles_outside_.resize(n_particles_, true);
  particles_outside_ions_.resize(n_ions_, true);
  ref_points_ = std::make_unique<NEReferencePoints>(input_.get_ref_points_input(), pset_dynamic_, pset_refs);
}

NEEnergyDensityEstimator::~NEEnergyDensityEstimator() {};

void NEEnergyDensityEstimator::registerListeners(QMCHamiltonian& ham_leader)
{
  ListenerVector<Real> kinetic_listener("kinetic", getListener(kinetic_values_));
  QMCHamiltonian::mw_registerKineticListener(ham_leader, kinetic_listener);
  ListenerVector<Real> potential_listener("potential", getListener(local_pot_values_));
  QMCHamiltonian::mw_registerLocalPotentialListener(ham_leader, potential_listener);
  ListenerVector<Real> ion_potential_listener("potential", getListener(local_ion_pot_values_));
  QMCHamiltonian::mw_registerLocalIonPotentialListener(ham_leader, ion_potential_listener);
}

/** This function collects the per particle energies.
 *  right now these are indentified by a string for each type.  This could be optimized but
 *  could also be an insiginificant cost versus the frequently large number of values handled.
 *  The values themselves are a vector of size particle_num.
 */
ListenerVector<QMCTraits::RealType>::ReportingFunction NEEnergyDensityEstimator::getListener(
    CrowdEnergyValues<Real>& local_values)
{
  return [&local_values](const int walker_index, const std::string& name, const Vector<Real>& inputV) {
    if (local_values.find(name) == local_values.end())
      local_values[name] = {};
    if (walker_index >= local_values[name].size())
      local_values[name].resize(walker_index + 1);
    local_values[name][walker_index] = inputV;
  };
}

const ParticleSet& NEEnergyDensityEstimator::getParticleSet(const PSPool& psetpool, const std::string& psname) const
{
  auto pset_iter(psetpool.find(psname));
  if (pset_iter == psetpool.end())
    throw UniformCommunicateError("Particle set pool does not contain \"" + psname +
                                  "\" so NEEnergyDensityEstimator::get_particleset fails!");
  return *(pset_iter->second.get());
}

std::size_t NEEnergyDensityEstimator::getFullDataSize() const
{
  std::size_t size = OperatorEstBase::get_data().size();
  for (const UPtr<NESpaceGrid<Real>>& grid : spacegrids_)
    size += grid->getDataVector().size();
  size += 1; // for nsamples;
  return size;
}

void NEEnergyDensityEstimator::packData(PooledData<Real>& buffer) const
{
  OperatorEstBase::packData(buffer);
  for (auto& grid : spacegrids_)
    buffer.add(grid->getDataVector().begin(), grid->getDataVector().end());
  buffer.add(nsamples_);
}

void NEEnergyDensityEstimator::unpackData(PooledData<Real>& buffer)
{
  OperatorEstBase::unpackData(buffer);
  for (auto& grid : spacegrids_)
    buffer.get(grid->getDataVector().begin(), grid->getDataVector().end());
  buffer.get(nsamples_);
}

// void NEEnergyDensityEstimator::unpackData(std::vector<Real>& operator_receive_buffer)
// {

// }

void NEEnergyDensityEstimator::accumulate(const RefVector<MCPWalker>& walkers,
                                          const RefVector<ParticleSet>& psets,
                                          const RefVector<TrialWaveFunction>& wfns,
                                          const RefVector<QMCHamiltonian>& hams,
                                          RandomBase<FullPrecReal>& rng)
{
  // Variable population is possible during DMC.
  reduced_local_kinetic_values_.resize(walkers.size());
  reduced_local_pot_values_.resize(walkers.size());
  reduced_local_ion_pot_values_.resize(walkers.size());

  // Depending on Hamiltonian setup one or more of these values could be absent.
  if (!local_pot_values_.empty())
    combinePerParticleEnergies(local_pot_values_, reduced_local_pot_values_);
  if (!kinetic_values_.empty())
    combinePerParticleEnergies(kinetic_values_, reduced_local_kinetic_values_);
  if (pset_static_)
  {
    if (!local_ion_pot_values_.empty())
      combinePerParticleEnergies(local_ion_pot_values_, reduced_local_ion_pot_values_);
  }
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    walkers_weight_ += walkers[iw].get().Weight;
    evaluate(psets[iw], walkers[iw], iw);
  }
  auto zero_reduced_values = [](auto& reduced_values) {
    for (auto& walker_values : reduced_values)
      std::fill(walker_values.begin(), walker_values.end(), 0.0);
  };
  zero_reduced_values(reduced_local_kinetic_values_);
  zero_reduced_values(reduced_local_pot_values_);
  zero_reduced_values(reduced_local_ion_pot_values_);
}

void NEEnergyDensityEstimator::evaluate(ParticleSet& pset, const MCPWalker& walker, const int walker_index)
{
  //Collect positions from ParticleSets
  int p_count = 0;
  {
    const ParticlePos& r_target = pset.R;
    for (int i = 0; i < r_target.size(); i++)
    {
      r_work_[p_count] = r_target[i];
      ++p_count;
    }
  }
  if (pset_static_ && !input_.get_ion_points())
  {
    const ParticlePos& r_static = pset_static_->R;
    for (int i = 0; i < r_static.size(); i++)
    {
      r_work_[p_count] = r_static[i];
      ++p_count;
    }
  }

  if (pset.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    pset.applyMinimumImage(r_work_);
  //Convert information accumulated in ParticleSets into EnergyDensity quantities
  Real weight = walker.Weight;
  p_count     = 0;
  {
    const auto& Ts = reduced_local_kinetic_values_[walker_index];
    const auto& Vd = reduced_local_pot_values_[walker_index];
    for (int i = 0; i < pset.getTotalNum(); i++)
    {
      ed_values_(p_count, W) = weight;
      ed_values_(p_count, T) = weight * Ts[i];
      ed_values_(p_count, V) = weight * Vd[i];
      ++p_count;
    }
  }
  if (pset_static_)
  {
    const ParticleSet& Ps = *pset_static_;
    auto& Vs              = reduced_local_ion_pot_values_[walker_index];
    if (!input_.get_ion_points())
    {
      // If there are no potentials involving ions only Vs will be
      // empty.
      for (auto ion_value : Vs)
      {
        ed_values_(p_count, W) = weight;
        ed_values_(p_count, T) = 0.0;
        ed_values_(p_count, V) = weight * ion_value;
        ++p_count;
      }
    }
    else
    {
      int ion_index{0};
      for (auto ion_value : Vs)
      {
        ed_ion_values_(ion_index, W) = weight;
        ed_ion_values_(ion_index, T) = 0.0;
        ed_ion_values_(ion_index, V) = weight * ion_value;
        ++ion_index;
      }
    }
  }
  //Accumulate energy density in spacegrids
  //const auto& dtab(pset.getDistTableAB(dtable_index_));
  fill(particles_outside_.begin(), particles_outside_.end(), true);
  for (int i = 0; i < spacegrids_.size(); i++)
  {
    NESpaceGrid<Real>& sg = *(spacegrids_[i]);
    sg.accumulate(r_work_, ed_values_, particles_outside_); //, dtab);
  }

  if (pset_static_ && input_.get_ion_points())
  {
    //
    fill(particles_outside_ions_.begin(), particles_outside_ions_.end(), true);
    for (int i = 0; i < spacegrids_.size(); i++)
    {
      NESpaceGrid<Real>& sg = *(spacegrids_[i]);
      sg.accumulate(r_ion_work_, ed_ion_values_, particles_outside_ions_);
    }
  }

  //Accumulate energy density of particles outside any spacegrid
  int bi, v;
  const int bimax = outside_buffer_offset + N_EDVALS;
  for (int p = 0; p < particles_outside_.size(); p++)
  {
    if (particles_outside_[p])
    {
      for (bi = outside_buffer_offset, v = 0; bi < bimax; bi++, v++)
      {
        data_[bi] += ed_values_(p, v);
      }
    }
  }
  if (input_.get_ion_points())
  {
    // Accumulate energy density for ions at a point field
    bi = outside_buffer_offset + N_EDVALS;
    for (int i = 0; i < n_ions_; i++)
      for (v = 0; v < N_EDVALS; v++, bi++)
      {
        data_[bi] += ed_ion_values_(i, v);
      }
  }
  ++nsamples_;
}

void NEEnergyDensityEstimator::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  int num_crowds = type_erased_operator_estimators.size();
  for (int ig = 0; ig < spacegrids_.size(); ++ig)
  {
    RefVector<NESpaceGrid<Real>> crowd_grids;
    crowd_grids.reserve(num_crowds);
    for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
    {
      NEEnergyDensityEstimator& crowd_ede = dynamic_cast<NEEnergyDensityEstimator&>(crowd_oeb);
      NESpaceGrid<Real>& grid_ref         = *(crowd_ede.spacegrids_[ig]);
      crowd_grids.push_back(grid_ref);
      //nsamples += grid_ref.samples;
    }
    NESpaceGrid<Real>::collect(*(spacegrids_[ig]), crowd_grids);
  }
  OperatorEstBase::collect(type_erased_operator_estimators);
}

void NEEnergyDensityEstimator::registerOperatorEstimator(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};

  hdf_path path_variables = hdf_name / std::string_view("variables");
  file.push(path_variables, true);
  file.write(n_particles_, "nparticles");
  auto nspacegrids = input_.get_space_grid_inputs().size();
  file.write(nspacegrids, "nspacegrids");
  file.write(nsamples_, "nsamples");

  if (input_.get_ion_points())
  {
    file.write(n_ions_, "nions");
  }
  file.pop();

  file.push(hdf_name);

  ref_points_->write(file);

  h5desc_.emplace_back(hdf_name / "outside");
  auto& ohOutside = h5desc_.back();
  std::vector<int> ng(1);
  ng[0] = N_EDVALS;
  ohOutside.set_dimensions(ng, outside_buffer_offset);
  for (int i = 0; i < spacegrids_.size(); i++)
    spacegrids_[i]->registerGrid(file, i);

  file.pop();
}

void NEEnergyDensityEstimator::write(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};
  file.push(hdf_name);
  file.write(particles_outside_, "outside");
  for (int i = 0; i < spacegrids_.size(); i++)
    spacegrids_[i]->write(file);
  if (input_.get_ion_points())
    file.write(ed_ion_values_, "ions");
  file.pop();
}

void NEEnergyDensityEstimator::zero()
{
  for (auto& spacegrid : spacegrids_)
    spacegrid->zero();
  walkers_weight_ = 0;
}

std::unique_ptr<OperatorEstBase> NEEnergyDensityEstimator::spawnCrowdClone() const
{
  auto spawn_data_locality = data_locality_;
  auto data_size           = this->data_.size();
  UPtr<NEEnergyDensityEstimator> spawn(std::make_unique<NEEnergyDensityEstimator>(*this, spawn_data_locality));
  spawn->get_data().resize(data_size);
  return spawn;
}

RefVector<NESpaceGrid<NEEnergyDensityEstimator::Real>> NEEnergyDensityEstimator::getSpaceGrids()
{
  return convertUPtrToRefVector(spacegrids_);
}

void NEEnergyDensityEstimator::normalize(QMCT::RealType invTotWgt)
{
  for (auto& spacegrid : spacegrids_)
    spacegrid->normalize(invTotWgt);
}

void NEEnergyDensityEstimator::startBlock(int steps) {}

// RefVector<std::vector<QMCT::RealType>>& NEEnergyDensityEstimator::getExtraData() {
//   RefVector<std::vector<QMCT::RealType>> refs;
//   for(auto& space_grid : space_grids_) {
//     refs.push_back(getDataVector());
//   }
//   return refs;
// }

} // namespace qmcplusplus
