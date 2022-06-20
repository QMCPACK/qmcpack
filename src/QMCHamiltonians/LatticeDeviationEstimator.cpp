//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "LatticeDeviationEstimator.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
LatticeDeviationEstimator::LatticeDeviationEstimator(ParticleSet& P,
                                                     ParticleSet& sP,
                                                     const std::string& tgroup_in,
                                                     const std::string& sgroup_in)
    : tspecies(P.getSpeciesSet()),
      sspecies(sP.getSpeciesSet()),
      tpset(P),
      spset(sP),
      tgroup(tgroup_in),
      sgroup(sgroup_in),
      hdf5_out(false),
      per_xyz(false),
      myTableID_(P.addTable(sP))
{
  // calculate number of source particles to use as lattice sites
  int src_species_id = sspecies.findSpecies(sgroup);
  num_sites          = spset.last(src_species_id) - spset.first(src_species_id);
  int tar_species_id = tspecies.findSpecies(tgroup);
  int num_tars       = tpset.last(tar_species_id) - tpset.first(tar_species_id);
  if (num_tars != num_sites)
  {
    app_log() << "number of target particles = " << num_tars << std::endl;
    app_log() << "number of source particles = " << num_sites << std::endl;
    APP_ABORT("nsource != ntargets in LatticeDeviationEstimator");
  }
}

bool LatticeDeviationEstimator::put(xmlNodePtr cur)
{
  input_xml             = cur;
  std::string hdf5_flag = "no";
  std::string xyz_flag  = "no";
  OhmmsAttributeSet attrib;
  attrib.add(hdf5_flag, "hdf5");
  attrib.add(xyz_flag, "per_xyz");
  attrib.put(cur);

  if (hdf5_flag == "yes")
  {
    hdf5_out = true;
  }
  else if (hdf5_flag == "no")
  {
    hdf5_out = false;
  }
  else
  {
    APP_ABORT("LatticeDeviationEstimator::put() - Please choose \"yes\" or \"no\" for hdf5 flag");
  } // end if hdf5_flag
  if (hdf5_out)
  {
    // change the default behavior of addValue() called by addObservables()
    // YY: does this still matter if I have overwritten addObservables()?
    update_mode_.set(COLLECTABLE, 1);
  }

  if (xyz_flag == "yes")
  {
    per_xyz = true;
    xyz2.resize(OHMMS_DIM);
  }
  else if (xyz_flag == "no")
  {
    per_xyz = false;
  }
  else
  {
    APP_ABORT("LatticeDeviationEstimator::put() - Please choose \"yes\" or \"no\" for per_xyz flag");
  } // end if xyz_flag

  return true;
}

bool LatticeDeviationEstimator::get(std::ostream& os) const
{ // class description
  os << "LatticeDeviationEstimator: " << name_ << "lattice = " << spset.getName();
  return true;
}

LatticeDeviationEstimator::Return_t LatticeDeviationEstimator::evaluate(ParticleSet& P)
{ // calculate <r^2> averaged over lattice sites
  value_ = 0.0;
  std::fill(xyz2.begin(), xyz2.end(), 0.0);

  RealType wgt        = t_walker_->Weight;
  const auto& d_table = P.getDistTableAB(myTableID_);

  // temp variables
  RealType r, r2;
  PosType dr;

  int nsite(0);    // site index
  int cur_jat(-1); // target particle index
  for (int iat = 0; iat < spset.getTotalNum(); iat++)
  { // for each desired source particle
    if (sspecies.speciesName[spset.GroupID[iat]] == sgroup)
    { // find desired species
      for (int jat = cur_jat + 1; jat < tpset.getTotalNum(); jat++)
      { // find corresponding (!!!! assume next) target particle
        if (tspecies.speciesName[tpset.GroupID[jat]] == tgroup)
        {
          // distance between particle iat in source pset, and jat in target pset
          r  = d_table.getDistRow(jat)[iat];
          r2 = r * r;
          value_ += r2;

          if (hdf5_out & !per_xyz)
          { // store deviration for each lattice site if h5 file is available
            P.Collectables[h5_index + nsite] = wgt * r2;
          }

          if (per_xyz)
          {
            dr = d_table.getDisplRow(jat)[iat];
            for (int idir = 0; idir < OHMMS_DIM; idir++)
            {
              RealType dir2 = dr[idir] * dr[idir];
              xyz2[idir] += dir2;
              if (hdf5_out)
              {
                P.Collectables[h5_index + nsite * OHMMS_DIM + idir] = wgt * dir2;
              }
            }
          }

          cur_jat = jat;
          break;
        }
      }
      nsite += 1; // count the number of sites, for checking only
    }             // found desired species (source particle)
  }

  if (nsite != num_sites)
  {
    app_log() << "num_sites = " << num_sites << " nsites = " << nsite << std::endl;
    APP_ABORT("Mismatch in LatticeDeivationEstimator.");
  }

  // average per site
  value_ /= num_sites;
  if (per_xyz)
  {
    for (int idir = 0; idir < OHMMS_DIM; idir++)
    {
      xyz2[idir] /= num_sites;
    }
  }

  return value_;
}

void LatticeDeviationEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  // get myIndex for scalar.dat
  if (per_xyz)
  {
    my_index_ = plist.size();
    for (int idir = 0; idir < OHMMS_DIM; idir++)
    {
      std::stringstream ss;
      ss << idir;
      plist.add(name_ + "_dir" + ss.str());
    }
  }
  else
  {
    my_index_ = plist.add(name_); // same as OperatorBase::addObservables
  }

  // get h5_index for stat.h5
  if (hdf5_out)
  {
    h5_index = collectables.size();
    std::vector<RealType> tmp;
    if (per_xyz)
    {
      tmp.resize(num_sites * OHMMS_DIM);
    }
    else
    {
      tmp.resize(num_sites);
    }
    collectables.add(tmp.begin(), tmp.end());
  }
}

void LatticeDeviationEstimator::setObservables(PropertySetType& plist)
{
  if (per_xyz)
  {
    copy(xyz2.begin(), xyz2.end(), plist.begin() + my_index_);
  }
  else
  {
    plist[my_index_] = value_; // default behavior
  }
}

void LatticeDeviationEstimator::resetTargetParticleSet(ParticleSet& P) {}

std::unique_ptr<OperatorBase> LatticeDeviationEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  // default constructor does not work with threads
  //LatticeDeviationEstimator* myclone = new LatticeDeviationEstimator(*this);
  std::unique_ptr<LatticeDeviationEstimator> myclone =
      std::make_unique<LatticeDeviationEstimator>(qp, spset, tgroup, sgroup);
  myclone->put(input_xml);
  return myclone;
}

void LatticeDeviationEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  if (hdf5_out)
  {
    // one scalar per lattice site (i.e. source particle)
    std::vector<int> ndim(1, num_sites);
    if (per_xyz)
    {
      // one scalar per lattice site per dimension
      ndim[0] = num_sites * OHMMS_DIM;
    }

    // open hdf5 entry and resize
    h5desc.emplace_back(name_);
    auto& h5o = h5desc.back();
    h5o.set_dimensions(ndim, h5_index);
    h5o.open(gid);
  }
}


} // namespace qmcplusplus
