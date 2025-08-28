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

namespace qmcplusplus
{

PairCorrelationEstimator::PairCorrelationEstimator(const PairCorrelationInput& pci,
                                                   const PSPool& psp,
                                                   DataLocality data_locality)
    : input_(pci)
{
  update_mode_.set(COLLECTABLE, 1);
  const auto& elecs = *(psp[input_.getTarget()]);
  num_species_      = elecs.groups();
  n_vec_.resize(num_species_, 0);
  for (int i = 0; i < num_species_; i++)
    n_vec_[i] = elecs.last(i) - elecs.first(i);
  n_e_ = elecs_.getTotalNum();

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

  if (input_.get_explicit_set_nbins() && !(input_.get_explicit_set_delta()))
  {
    num_bins_ = input_.get_nbins();)
    delta = rmax_ / static_cast<Real>(num_bins_);
  }
  else if (!(input_.get_explicit_set_nbins()) && input_.get_explicit_set_delta())
  {
    delta_    = input_.get_delta();
    num_bins_ = static_cast<int>(rmax_ / delta_);
  }
  else if (input_.get_explicit_set_nbins() && input_.get_explicit_set_delta())
  {
    // Here we must assume that certain checks were done in parsing.
    delta     = input_.get_delta();
    num_bins_ = input_.get_nbins();
    if (rmax_ != delta * num_bins_)
      throw UniformCommunicateError("Explicitly input PairCorrelation delta and num_bins inconsistent with rmax!");
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
      gof_r_prefix.push_back(os.str());
      pair_map[i * num_species + j] = npairs;
      ++npairs;
    }

  // source-target tables
  std::vector<std::string> slist, dlist;
  const int ntables = elecs.getNumDistTables();
  for (int k = 0; k < ntables; ++k)
    if (elecs.getName() != elecs.getDistTable(k).get_origin().getName())
      dlist.push_back(elecs.getDistTable(k).get_origin().getName());

  std::set<int> others_sorted;
  for (int i = 0; i < sources_.size(); ++i)
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
    const DistanceTable& t(elns.getDistTable(other_ids_[k]));
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

  // Set normalization based on updated parameters & report status
  set_norm_factor();
  report();
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
  BufferType& collectables(P.Collectables);
  const auto& dii(P.getDistTableAA(d_aa_ID_));
  for (int iat = 1; iat < dii.centers(); ++iat)
  {
    const auto& dist = dii.getDistRow(iat);
    const int ig     = P.GroupID[iat];
    for (int j = 0; j < iat; ++j)
    {
      const RealType r = dist[j];
      if (r < Dmax)
      {
        const int loc     = static_cast<int>(DeltaInv * r);
        const int jg      = P.GroupID[j];
        const int pair_id = gen_pair_id(ig, jg, num_species);
        collectables[pair_id * NumBins + loc + my_index_] += norm_factor(pair_id + 1, loc);
      }
    }
  }
  for (int k = 0; k < other_ids.size(); ++k)
  {
    const auto& d1(P.getDistTableAB(other_ids[k]));
    const ParticleSet::ParticleIndex& gid(d1.get_origin().GroupID);
    int koff        = other_offsets[k];
    RealType overNI = 1.0 / d1.centers();
    for (int iat = 0; iat < d1.targets(); ++iat)
    {
      const auto& dist = d1.getDistRow(iat);
      for (int j = 0; j < d1.centers(); ++j)
      {
        const RealType r = dist[j];
        if (r < Dmax)
        {
          int toff = (gid[j] + koff) * NumBins;
          int loc  = static_cast<int>(DeltaInv * r);
          collectables[toff + loc + my_index_] += norm_factor(0, loc) * overNI;
        }
      }
    }
  }
  return 0.0;
}

void PairCorrEstimator::registerCollectables(std::vector<ObservableHelper>& h5list, hdf_archive& file) const
{
  std::vector<int> onedim(1, NumBins);
  int offset = my_index_;
  for (int i = 0; i < gof_r_prefix.size(); ++i)
  {
    h5list.emplace_back(hdf_path{gof_r_prefix[i]});
    auto& h5o = h5list.back();
    h5o.set_dimensions(onedim, offset);
    h5o.addProperty(const_cast<RealType&>(Delta), "delta", file);
    h5o.addProperty(const_cast<RealType&>(Dmax), "cutoff", file);
    offset += NumBins;
  }
}


void PairCorrEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.size();
  std::vector<RealType> g(gof_r_prefix.size() * NumBins, 0);
  collectables.add(g.begin(), g.end());
  ////only while debugging
  //if(gof_r.size())
  //{
  //  myDebugIndex=plist.size();
  //  for(int i=0; i<gof_r_prefix.size(); ++i)
  //  {
  //    for(int k=0; k<gof_r.cols(); ++k)
  //    {
  //      std::ostringstream h;
  //      h << gof_r_prefix[i]<< "_" << k;
  //      int dum=plist.add(h.str());
  //    }
  //  }
  //}
}

void PairCorrelationEstimator::set_norm_factor()
{
  /*
     Number of species-pair-specific gofr's to compute
     E.g. "uu", "dd", "ud", & etc.
     Note the addition +1 bin, which is for the total gofr of the system
     (which seems not to get saved to h5 file. Why compute it then?)
  */
  const RealType n_channels = num_species_ * (num_species_ - 1) / 2 + num_species_ + 1;
  norm_factor_.resize(n_channels, num_bins_);

  /*
     Compute the normalization V/Npairs/Nid, with
     V the volume of the system
     Npairs the number of (unique) pairs of particles of given types
     Nid the number of particles expected for a uniformly random distribution
     with the same number density
  */
  RealType r_bin_start   = 0.;
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
      for (int n = m; n < num_species; n++)
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


void PairCorrEstimator::report()
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

bool PairCorrEstimator::get(std::ostream& os) const
{
  os << name_ << " dmax=" << rmax_ << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> PairCorrEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //default constructor is sufficient
  return std::make_unique<PairCorrEstimator>(*this);
}


} // namespace qmcplusplus
