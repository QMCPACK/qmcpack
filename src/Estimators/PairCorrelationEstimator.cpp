//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/PairCorrEstimator.h
//////////////////////////////////////////////////////////////////////////////////////

#include "PairCorrelationEstimator.h"
#include "OperatorEstBase.h"
#include <DistanceTable.h>

namespace qmcplusplus
{

PairCorrelationEstimator::PairCorrelationEstimator(const PairCorrelationInput& pci,
                                                   const PSPool& pset_pool,
                                                   ParticleSet& elecs,
                                                   DataLocality data_locality)
    : OperatorEstBase(data_locality, pci.get_name(), pci.get_type()), input_(pci)
{
  num_species_ = elecs.groups();
  n_vec_.resize(num_species_, 0);
  for (int i = 0; i < num_species_; i++)
    n_vec_[i] = elecs.last(i) - elecs.first(i);
  n_e_ = elecs.getTotalNum();

  d_aa_id_ = elecs.addTable(elecs, DTModes::NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP);

  // legacy behavior is that the input rmax wins these are overwritten by input tags for rmax
  if (elecs.getLattice().SuperCellEnum)
  {
    rmax_   = input_.get_explicit_set_rmax() ? input_.get_rmax() : elecs.getLattice().WignerSeitzRadius;
    volume_ = elecs.getLattice().Volume;
  }
  else // Open BC's
  {
    rmax_   = input_.get_rmax();
    volume_ = 1.0;
  }

  // This reproduces as well I as I could figure it out the behavior
  // of the legacy gofr wrt different combinations of input atrributes
  if (input_.get_explicit_set_nbins() && !(input_.get_explicit_set_delta()))
  {
    num_bins_ = input_.get_nbins();
    delta_    = rmax_ / static_cast<Real>(num_bins_);
  }
  else if (!(input_.get_explicit_set_nbins()) && input_.get_explicit_set_delta())
  {
    delta_    = input_.get_delta();
    num_bins_ = static_cast<int>(rmax_ / delta_);
  }
  else if (input_.get_explicit_set_nbins() && input_.get_explicit_set_delta())
  {
    // Here we must assume that certain checks were done in parsing.
    delta_               = input_.get_delta();
    num_bins_            = input_.get_nbins();
    Real calculated_rmax = delta_ * num_bins_;
    if (!(input_.get_explicit_set_rmax()))
      rmax_ = calculated_rmax;
    else
    {
      Real rmax_out_by = std::abs(rmax_ - calculated_rmax);
      if (rmax_out_by > delta_ / 2)
        throw UniformCommunicateError("Explicitly input PairCorrelation delta and num_bins product(" +
                                      std::to_string(calculated_rmax) + ") inconsistent with rmax (" +
                                      std::to_string(rmax_) + ")!");
    }
  }
  else // neither nbins nor delta set use the default num_bins_ to determine delta
  {
    num_bins_ = input_.get_nbins();
    delta_    = rmax_ / static_cast<Real>(num_bins_);
  }

  delta_inv_ = 1.0 / delta_;

  std::map<int, int> pair_map;
  int npairs = 0;
  for (int i = 0; i < num_species_; ++i)
    for (int j = i; j < num_species_; ++j)
    {
      std::ostringstream os;
      os << "gofr_" << elecs.getName() << "_" << i << "_" << j;
      gof_r_prefix_.push_back(os.str());
      pair_map[i * num_species_ + j] = npairs;
      ++npairs;
    }

  // source-target tables
  std::vector<std::string> slist, dlist;
  const int ntables = elecs.getNumDistTables();
  for (int k = 0; k < ntables; ++k)
    if (elecs.getName() != elecs.getDistTable(k).get_origin().getName())
      dlist.push_back(elecs.getDistTable(k).get_origin().getName());

  std::set<int> others_sorted;
  auto& sources = input_.get_sources();
  for (int i = 0; i < sources.size(); ++i)
  {
    int k = 0;
    while (k < dlist.size())
    {
      if (slist[i] == dlist[k])
      {
        others_sorted.insert(k + 1);
        break;
      }
      ++k;
    }
  }
  other_ids_.resize(others_sorted.size());
  other_offsets_.resize(others_sorted.size());
  copy(others_sorted.begin(), others_sorted.end(), other_ids_.begin());
  int toff = gof_r_prefix_.size();
  for (int k = 0; k < other_ids_.size(); ++k)
  {
    const DistanceTable& t(elecs.getDistTable(other_ids_[k]));
    app_log() << "  GOFR for " << t.getName() << " starts at " << toff << std::endl;
    other_offsets_[k] = toff;
    const SpeciesSet& species(t.get_origin().getSpeciesSet());
    int ng = species.size();
    for (int i = 0; i < ng; ++i)
    {
      std::ostringstream os;
      os << "gofr_" << t.getName() << "_" << species.speciesName[i];
      gof_r_prefix_.push_back(os.str());
    }
    toff += ng;
  }

  data_.resize(num_bins_ * toff);
  // Set normalization based on updated parameters & report status
  set_norm_factor();
  report();
}

PairCorrelationEstimator::PairCorrelationEstimator(const PairCorrelationEstimator& pce, DataLocality dl)
    : qmcplusplus::PairCorrelationEstimator(pce)
{
  data_locality_ = dl;
  data_.resize(pce.data_.size());
}

// The value should match the index to norm_factor in set_norm_factor
int PairCorrelationEstimator::gen_pair_id(const int ig, const int jg, const int ns)
{
  if (jg < ig)
    return ns * (ns - 1) / 2 - (ns - jg) * (ns - jg - 1) / 2 + ig;
  else
    return ns * (ns - 1) / 2 - (ns - ig) * (ns - ig - 1) / 2 + jg;
}

void PairCorrelationEstimator::accumulate(const RefVector<MCPWalker>& walkers,
                                          const RefVector<ParticleSet>& psets,
                                          const RefVector<TrialWaveFunction>& wfns,
                                          const RefVector<QMCHamiltonian>& hams,
                                          RandomBase<FullPrecReal>& rng)
{
  for (int iw = 0; iw < walkers.size(); ++iw)
  {
    Real weight       = walkers[iw].get().Weight;
    ParticleSet& pset = psets[iw];
    const auto& dii(pset.getDistTableAA(d_aa_id_));
    for (int iat = 1; iat < dii.centers(); ++iat)
    {
      const auto& dist = dii.getDistRow(iat);
      const int ig     = pset.GroupID[iat];
      for (int j = 0; j < iat; ++j)
      {
        const Real r = dist[j];
        if (r < rmax_)
        {
          const int loc     = static_cast<int>(delta_inv_ * r);
          const int jg      = pset.GroupID[j];
          const int pair_id = gen_pair_id(ig, jg, num_species_);
          data_[pair_id * num_bins_ + loc] += norm_factor_(pair_id + 1, loc);
        }
      }
    }
    for (int k = 0; k < other_ids_.size(); ++k)
    {
      const auto& d1(pset.getDistTableAB(other_ids_[k]));
      const ParticleSet::ParticleIndex& gid(d1.get_origin().GroupID);
      int koff    = other_offsets_[k];
      Real overNI = 1.0 / d1.centers();
      for (int iat = 0; iat < d1.targets(); ++iat)
      {
        const auto& dist = d1.getDistRow(iat);
        for (int j = 0; j < d1.centers(); ++j)
        {
          const Real r = dist[j];
          if (r < rmax_)
          {
            int toff = (gid[j] + koff) * num_bins_;
            int loc  = static_cast<int>(delta_inv_ * r);
            data_[toff + loc] += norm_factor_(0, loc) * overNI;
          }
        }
      }
    }
    walkers_weight_ += weight;
  }
}

ParticleSet& PairCorrelationEstimator::getParticleSet(PSPool& pset_pool, const std::string& psname) const
{
  auto pset_iter(pset_pool.find(psname));
  if (pset_iter == pset_pool.end())
  {
    throw UniformCommunicateError("Particle set pool does not contain \"" + psname +
                                  "\" so StructureFactorEstimator::get_particleset fails!");
  }
  return *(pset_iter->second.get());
}

void PairCorrelationEstimator::registerOperatorEstimator(hdf_archive& file) {}

void PairCorrelationEstimator::write(hdf_archive& file)
{
  hdf_path hdf_name{my_name_};
  file.push(hdf_name);
  std::vector<int> onedim(1, num_bins_);

  for (int i = 0; i < gof_r_prefix_.size(); ++i)
  {
    hdf_path prefix_path{gof_r_prefix_[i]};
    file.push(prefix_path);
    file.write(delta_, "delta");
    file.write(rmax_, "cutoff");
    Vector prefix_vector(data_.data() + num_bins_ * i, num_bins_);
    file.write(prefix_vector, "values");
    file.pop();
  }
  file.pop();
}

void PairCorrelationEstimator::resize(int nbins)
{
  num_bins_ = nbins;
  norm_factor_.resize((num_species_ * num_species_ - num_species_) / 2 + num_species_ + 1, num_bins_);
  Real r  = delta_ * 0.5;
  Real pf = volume_ * delta_inv_ / (4 * M_PI);
  for (int i = 0; i < num_bins_; ++i, r += delta_)
  {
    Real rm2           = pf / r / r;
    norm_factor_(0, i) = rm2 / n_e_;
    int indx(1);
    for (int m(0); m < num_species_; m++)
      for (int n(m); n < num_species_; n++, indx++)
        norm_factor_(indx, i) = rm2 * n_vec_[n] * n_vec_[m];
  }
}

void PairCorrelationEstimator::set_norm_factor()
{
  /*
     Number of species-pair-specific gofr's to compute
     E.g. "uu", "dd", "ud", & etc.
     Note the addition +1 bin, which is for the total gofr of the system
     (which seems not to get saved to h5 file. Why compute it then?)
  */
  const Real n_channels = num_species_ * (num_species_ - 1) / 2 + num_species_ + 1;
  norm_factor_.resize(n_channels, num_bins_);

  /*
     Compute the normalization V/Npairs/Nid, with
     V the volume of the system
     Npairs the number of (unique) pairs of particles of given types
     Nid the number of particles expected for a uniformly random distribution
     with the same number density
  */
  Real r_bin_start       = 0.;
  const Real ftpi        = 4. / 3 * M_PI;
  const Real n_tot_pairs = n_e_ * (n_e_ - 1) / 2;
  for (int i = 0; i < num_bins_; i++)
  {
    // Separation for current bin
    r_bin_start = static_cast<Real>(i) * delta_;

    // Number density of pairs
    Real rho = n_tot_pairs / volume_;

    // Volume of spherical shell of thickness Delta
    Real bin_volume = ftpi * (std::pow(r_bin_start + delta_, 3) - std::pow(r_bin_start, 3));

    // Expected number of pairs of particles separated by distance r assuming
    // they are uniformly randomly distributed (ideal gas-like)
    Real nid           = rho * bin_volume;
    norm_factor_(0, i) = 1. / nid;
    int indx{1};

    // Do same as above, but for each unique pair of species
    // e.g. uu, ud, dd...
    for (int m = 0; m < num_species_; m++)
    {
      const Real nm = n_vec_[m];
      for (int n = m; n < num_species_; n++)
      {
        const Real nn         = n_vec_[n];
        const Real npairs     = (m == n ? (nn * (nn - 1) / 2.) : (nn * nm));
        rho                   = npairs / volume_;
        nid                   = rho * bin_volume;
        norm_factor_(indx, i) = 1. / nid;
        indx++;
      }
    }
  }
}

void PairCorrelationEstimator::startBlock(int steps) {}

void PairCorrelationEstimator::report()
{
  app_log() << "PairCorrEstimator report" << std::endl;
  app_log() << "  num_species = " << num_species_ << std::endl;
  app_log() << "  Volume      = " << volume_ << std::endl;
  app_log() << "  Dmax        = " << rmax_ << std::endl;
  app_log() << "  NumBins     = " << num_bins_ << std::endl;
  app_log() << "  Delta       = " << delta_ << std::endl;
  app_log() << "  DeltaInv    = " << delta_inv_ << std::endl;
  //app_log()<<"  x = "<< x << std::endl;
  app_log() << "end PairCorrEstimator report" << std::endl;
}

UPtr<OperatorEstBase> PairCorrelationEstimator::spawnCrowdClone() const
{
  UPtr<PairCorrelationEstimator> spawn(std::make_unique<PairCorrelationEstimator>(*this, data_locality_));
  //default constructor is sufficient
  return spawn;
}


} // namespace qmcplusplus
