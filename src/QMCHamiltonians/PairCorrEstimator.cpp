//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "PairCorrEstimator.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/SimpleParser.h"
#include <set>

namespace qmcplusplus
{
PairCorrEstimator::PairCorrEstimator(ParticleSet& elns, std::string& sources)
    : Dmax(10.), Delta(0.5), num_species(2), d_aa_ID_(elns.addTable(elns))
{
  update_mode_.set(COLLECTABLE, 1);
  num_species = elns.groups();
  n_vec.resize(num_species, 0);
  for (int i = 0; i < num_species; i++)
    n_vec[i] = elns.last(i) - elns.first(i);
  N_e = elns.getTotalNum();

  // use the simulation cell radius if any direction is periodic
  if (elns.getLattice().SuperCellEnum)
  {
    Dmax   = elns.getLattice().WignerSeitzRadius;
    Volume = elns.getLattice().Volume;
  }
  else // Open BC's
    Volume = 1.0;
  NumBins  = static_cast<int>(Dmax / Delta);
  Delta    = Dmax / static_cast<RealType>(NumBins);
  DeltaInv = 1.0 / Delta;
  //ostringstream h;
  //h<<"gofr_" << elns.getName();
  //gof_r_prefix.push_back(h.str());
  std::map<int, int> pair_map;
  int npairs = 0;
  for (int i = 0; i < num_species; ++i)
    for (int j = i; j < num_species; ++j)
    {
      std::ostringstream os;
      os << "gofr_" << elns.getName() << "_" << i << "_" << j;
      gof_r_prefix.push_back(os.str());
      pair_map[i * num_species + j] = npairs;
      ++npairs;
    }

  // source-target tables
  std::vector<std::string> slist, dlist;
  const int ntables = elns.getNumDistTables();
  for (int k = 0; k < ntables; ++k)
    if (elns.getName() != elns.getDistTable(k).get_origin().getName())
      dlist.push_back(elns.getDistTable(k).get_origin().getName());
  parsewords(sources.c_str(), slist);
  std::set<int> others_sorted;
  for (int i = 0; i < slist.size(); ++i)
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
  other_ids.resize(others_sorted.size());
  other_offsets.resize(others_sorted.size());
  copy(others_sorted.begin(), others_sorted.end(), other_ids.begin());
  int toff = gof_r_prefix.size();
  for (int k = 0; k < other_ids.size(); ++k)
  {
    const DistanceTable& t(elns.getDistTable(other_ids[k]));
    app_log() << "  GOFR for " << t.getName() << " starts at " << toff << std::endl;
    other_offsets[k] = toff;
    const SpeciesSet& species(t.get_origin().getSpeciesSet());
    int ng = species.size();
    for (int i = 0; i < ng; ++i)
    {
      std::ostringstream os;
      os << "gofr_" << t.getName() << "_" << species.speciesName[i];
      gof_r_prefix.push_back(os.str());
    }
    toff += ng;
  }
}

void PairCorrEstimator::resetTargetParticleSet(ParticleSet& P) {}

// The value should match the index to norm_factor in set_norm_factor
int PairCorrEstimator::gen_pair_id(const int ig, const int jg, const int ns)
{
  if (jg < ig)
    return ns * (ns - 1) / 2 - (ns - jg) * (ns - jg - 1) / 2 + ig;
  else
    return ns * (ns - 1) / 2 - (ns - ig) * (ns - ig - 1) / 2 + jg;
}

PairCorrEstimator::Return_t PairCorrEstimator::evaluate(ParticleSet& P)
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

void PairCorrEstimator::registerCollectables(std::vector<ObservableHelper>& h5list, hid_t gid) const
{
  std::vector<int> onedim(1, NumBins);
  int offset = my_index_;
  for (int i = 0; i < gof_r_prefix.size(); ++i)
  {
    h5list.emplace_back(gof_r_prefix[i]);
    auto& h5o = h5list.back();
    h5o.set_dimensions(onedim, offset);
    h5o.open(gid);
    h5o.addProperty(const_cast<RealType&>(Delta), "delta");
    h5o.addProperty(const_cast<RealType&>(Dmax), "cutoff");
    //       h5o->addProperty(const_cast<std::vector<RealType>&>(norm_factor),"norm_factor");
    //       std::string blob("norm_factor[i]=1/r_m[i]^2 for r_m[i]=(r[i]+r[i+1])/2");
    //       h5o->addProperty(blob,"dictionary");
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


void PairCorrEstimator::setObservables(PropertySetType& plist)
{
  //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex);
}

void PairCorrEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myDebugIndex+offset);
}

bool PairCorrEstimator::put(xmlNodePtr cur)
{
  //set resolution
  int nbins = (int)std::ceil(Dmax * DeltaInv);
  std::string debug("no");
  OhmmsAttributeSet attrib;
  attrib.add(nbins, "num_bin");
  attrib.add(Dmax, "rmax");
  attrib.add(Delta, "dr");
  attrib.add(debug, "debug");
  attrib.put(cur);
  Delta    = Dmax / static_cast<RealType>(nbins);
  DeltaInv = 1.0 / Delta;
  NumBins  = nbins;

  // Set normalization based on updated parameters & report status
  set_norm_factor();
  report();

  return true;
}


// Called from inside put() after internals are set
// Sets the normalization, or norm_factor, for each channel and bin
void PairCorrEstimator::set_norm_factor()
{
  /*
     Number of species-pair-specific gofr's to compute
     E.g. "uu", "dd", "ud", & etc.
     Note the addition +1 bin, which is for the total gofr of the system
     (which seems not to get saved to h5 file. Why compute it then?)
  */
  const RealType n_channels = num_species * (num_species - 1) / 2 + num_species + 1;
  norm_factor.resize(n_channels, NumBins);

  /*
     Compute the normalization V/Npairs/Nid, with
     V the volume of the system
     Npairs the number of (unique) pairs of particles of given types
     Nid the number of particles expected for a uniformly random distribution
     with the same number density
  */
  RealType r                 = 0.;
  const RealType ftpi        = 4. / 3 * M_PI;
  const RealType N_tot_pairs = N_e * (N_e - 1) / 2;
  for (int i = 0; i < NumBins; i++)
  {
    // Separation for current bin
    r = static_cast<RealType>(i) * Delta;

    // Number density of pairs
    RealType rho = N_tot_pairs / Volume;

    // Volume of spherical shell of thickness Delta
    RealType bin_volume = ftpi * (std::pow(r + Delta, 3) - std::pow(r, 3));

    // Expected number of pairs of particles separated by distance r assuming
    // they are uniformly randomly distributed (ideal gas-like)
    RealType nid      = rho * bin_volume;
    norm_factor(0, i) = 1. / nid;
    int indx(1);

    // Do same as above, but for each unique pair of species
    // e.g. uu, ud, dd...
    for (int m = 0; m < num_species; m++)
    {
      const RealType nm = n_vec[m];
      for (int n = m; n < num_species; n++)
      {
        const RealType nn     = n_vec[n];
        const RealType npairs = (m == n ? (nn * (nn - 1) / 2.) : (nn * nm));
        rho                   = npairs / Volume;
        nid                   = rho * bin_volume;
        norm_factor(indx, i)  = 1. / nid;
        indx++;
      }
    }
  }

  // ***DEBUG*** print norm_factor
  /*
  std::cout << "norm_factor:\n";
  std::cout << std::fixed;
  for( int j=0; j<norm_factor.size2(); j++ )
    {
	std::cout << std::setw(4) << j;
	std::cout << std::setw(8) << std::setprecision(4) << j*Delta;
	for( int i=0; i<norm_factor.size1(); i++ )
	  {
	    std::cout << "  " << std::setw(10) << std::setprecision(4) << norm_factor(i,j);
	  }
	std::cout << std::endl;
    }
  std::cout << std::endl;
  */
}


void PairCorrEstimator::report()
{
  app_log() << "PairCorrEstimator report" << std::endl;
  app_log() << "  num_species = " << num_species << std::endl;
  app_log() << "  Volume      = " << Volume << std::endl;
  app_log() << "  Dmax        = " << Dmax << std::endl;
  app_log() << "  NumBins     = " << NumBins << std::endl;
  app_log() << "  Delta       = " << Delta << std::endl;
  app_log() << "  DeltaInv    = " << DeltaInv << std::endl;
  //app_log()<<"  x = "<< x << std::endl;
  app_log() << "end PairCorrEstimator report" << std::endl;
}

bool PairCorrEstimator::get(std::ostream& os) const
{
  os << name_ << " dmax=" << Dmax << std::endl;
  return true;
}

std::unique_ptr<OperatorBase> PairCorrEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //default constructor is sufficient
  return std::make_unique<PairCorrEstimator>(*this);
}

void PairCorrEstimator::resize(int nbins)
{
  NumBins = nbins;
  norm_factor.resize((num_species * num_species - num_species) / 2 + num_species + 1, NumBins);
  RealType r  = Delta * 0.5;
  RealType pf = Volume * DeltaInv / (4 * M_PI);
  for (int i = 0; i < NumBins; ++i, r += Delta)
  {
    RealType rm2      = pf / r / r;
    norm_factor(0, i) = rm2 / N_e;
    int indx(1);
    for (int m(0); m < num_species; m++)
      for (int n(m); n < num_species; n++, indx++)
        norm_factor(indx, i) = rm2 * n_vec[n] * n_vec[m];
  }
  //     norm_factor.resize(nbins);
  //     RealType r=Delta*0.5;
  //     for(int i=0; i<norm_factor.size(); ++i, r+=Delta)
  //     {
  //       norm_factor[i]=1.0/r/r;
  //     }
}
} // namespace qmcplusplus
