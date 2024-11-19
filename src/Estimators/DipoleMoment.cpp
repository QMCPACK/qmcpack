//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "DipoleMoment.h"
#include "TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"

#include <iostream>
#include <numeric>


namespace qmcplusplus
{
DipoleMoment::DipoleMoment(DipoleMomentInput&& inp_, const PSPool& pset_pool, DataLocality dl)
    : OperatorEstBase(dl), input_(std::move(inp_))
{
  auto& inp = this->input_.input_section_;
  
  // pull data from inputs
  my_name_            = inp.get<std::string>("name");
  auto ion_pset_name  = inp.get<std::string>("ions");
  
  // extricate the ion particleset
  auto pit = pset_pool.find(ion_pset_name);
  if (pit == pset_pool.end())
    throw std::runtime_error("DipoleMoment estimator: ion particleset"+ion_pset_name+" does not exist");
  ParticleSet& pion  = *pit->second.get();
  
  // dipole moment is only defined (here) for open bc's
  if (pion.getLattice().SuperCellEnum != SUPERCELL_OPEN)
    throw std::runtime_error("DipoleMoment estimator is only compatible with open boundary conditions");
  
  // compute and store the ion contribution to the dipole moment
  //   resize and zero ion_dipole_moment_
  ion_dipole_moment_.resize(DIM);
  for (int d=0; d<DIM; d++)
    ion_dipole_moment_[d] = 0.0;
  //   heavily unpack ion charges from obtuse data structures
  SpeciesSet& species  = pion.getSpeciesSet();
  int ChargeAttribIndx = species.addAttribute("charge");
  int nspecies         = species.TotalNum;
  int nps              = pion.getTotalNum();
  assert(nps==P.R.size());
  std::vector<RealType> Zptcl;
  std::vector<RealType> Zspec;
  Zspec.resize(nspecies);
  Zptcl.resize(nps);
  for (int spec = 0; spec < nspecies; spec++)
    Zspec[spec] = species(ChargeAttribIndx, spec);
  for (int i = 0; i < nps; i++)
    Zptcl[i] = Zspec[pion.GroupID[i]];
  //   accumulate ion_dipole_moment_
  for (int i = 0; i < pion.R.size(); i++)
    for (int d=0; d<DIM; d++)
      ion_dipole_moment_[d] += Zptcl[i]*pion.R[i][d];
  
  // resize estimator data (effectively electron_dipole_moment_)
  const size_t data_size = DIM;
  data_.resize(data_size, 0.0);
}


DipoleMoment::DipoleMoment(const DipoleMoment& sh, DataLocality dl) : DipoleMoment(sh)
{
  data_locality_ = dl;
}

std::unique_ptr<OperatorEstBase> DipoleMoment::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for DipoleMoment");
  }

  auto spawn = std::make_unique<DipoleMoment>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size);
  return spawn;
}

void DipoleMoment::startBlock(int steps) {}

/** Gets called every step and writes to thread local data.
 *
 */
void DipoleMoment::accumulate(const RefVector<MCPWalker>& walkers,
                                    const RefVector<ParticleSet>& psets,
                                    const RefVector<TrialWaveFunction>& wfns,
                                    const RefVector<QMCHamiltonian>& hams,
                                    RandomBase<FullPrecRealType>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    MCPWalker& walker  = walkers[iw];
    RealType weight    = walker.Weight;
    
    // accumulate weight (required by all estimators, otherwise inf results)
    walkers_weight_ += weight;
    
    // accumulate data
    //   ion contribution
    for (int d = 0; d < DIM; ++d)
      data_[d] += weight * ion_dipole_moment_[d];
    //   electron contribution
    RealType elec_charge = -1.0;
    assert(det_ratios.size() == data_.size());
    for (int i = 0; i < walker.R.size(); i++)
      for (int d = 0; d < DIM; ++d)
        data_[d] += weight * elec_charge * walker.R[i][d];
  }
}


void DipoleMoment::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a DipoleMoment with this DataLocality");
  }
}


void DipoleMoment::registerOperatorEstimator(hdf_archive& file)
{
  using namespace std::string_literals;
  ////descriptor for the data, 1-D data
  std::vector<int> ng(1);
  ng[0] = data_.size();
  h5desc_.push_back({{my_name_}});
  auto& h5o = h5desc_.back();
  h5o.set_dimensions(ng, 0);
}


} // namespace qmcplusplus
