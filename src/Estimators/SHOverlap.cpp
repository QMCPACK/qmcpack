//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "SHOverlap.h"
#include "TrialWaveFunction.h"

#include <iostream>
#include <numeric>


namespace qmcplusplus
{
SHOverlap::SHOverlap(SHOverlapInput&& inp_,
                     DataLocality dl)
    : OperatorEstBase(dl),
      input_(std::move(inp_))
{

  //my_name_ = input_.get_name();

  auto& inp = this->input_.input_section_;

  //auto other = inp.get<std::string>("other");
  //app_log()<<"other: "<<other<<std::endl;


  // allocate data storage
  //size_t data_size = nofK.size();
  //data_.resize(data_size, 0.0);


}



SHOverlap::SHOverlap(const SHOverlap& sh, DataLocality dl) : SHOverlap(sh)
{
  data_locality_ = dl;
}

std::unique_ptr<OperatorEstBase> SHOverlap::spawnCrowdClone() const
{
  std::size_t data_size    = data_.size();
  auto spawn_data_locality = data_locality_;

  if (data_locality_ == DataLocality::rank)
  {
    // This is just a stub until a memory saving optimization is deemed necessary
    spawn_data_locality = DataLocality::queue;
    data_size           = 0;
    throw std::runtime_error("There is no memory savings implementation for SHOverlap");
  }

  auto spawn = std::make_unique<SHOverlap>(*this, spawn_data_locality);
  spawn->get_data().resize(data_size);
  return spawn;
}

void SHOverlap::startBlock(int steps)
{
}

/** Gets called every step and writes to thread local data.
 *
 */
void SHOverlap::accumulate(const RefVector<MCPWalker>& walkers,
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
    auto& wcs              = psi.getOrbitals();

    const int np = pset.getTotalNum();

    // separate jastrow and fermi wavefunction components
    std::vector<WaveFunctionComponent*> wcs_jastrow;
    std::vector<WaveFunctionComponent*> wcs_fermi;
    for(auto& wc: wcs)
      if(wc->isFermionic())
        wcs_fermi.push_back(wc.get());
      else
        wcs_jastrow.push_back(wc.get());

    // accumulate weight
    //  (required by all estimators, otherwise inf results)
    walkers_weight_ += weight;

    //// accumulate data
    //for (int ik = 0; ik < nofK.size(); ++ik)
    //  data_[ik] += weight * nofK[ik] * norm_nofK;


    
  }

  throw std::runtime_error("SHOverlap accumulate");
    
}


void SHOverlap::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  if (data_locality_ == DataLocality::crowd)
  {
    OperatorEstBase::collect(type_erased_operator_estimators);
  }
  else
  {
    throw std::runtime_error("You cannot call collect on a SHOverlap with this DataLocality");
  }
}


void SHOverlap::registerOperatorEstimator(hdf_archive& file)
{
  using namespace std::string_literals;
  ////descriptor for the data, 1-D data
  //std::vector<int> ng(1);
  ////add nofk
  //ng[0] = nofK.size();
  //h5desc_.push_back({{"nofk"s}});
  //auto& h5o = h5desc_.back();
  ////h5o.set_dimensions(ng, my_index_);
  //h5o.set_dimensions(ng, 0); // JTK: doesn't seem right
  //h5o.addProperty(const_cast<std::vector<PosType>&>(kPoints), "kpoints", file);
  //h5o.addProperty(const_cast<std::vector<int>&>(kWeights), "kweights", file);
}


} // namespace qmcplusplus
