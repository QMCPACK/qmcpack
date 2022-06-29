//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SkEstimator.h"
#include "LongRange/StructFact.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
SkEstimator::SkEstimator(ParticleSet& source)
{
  sourcePtcl = &source;
  update_mode_.set(COLLECTABLE, 1);
  NumSpecies = source.getSpeciesSet().getTotalNum();
  NumK       = source.getSimulationCell().getKLists().numk;
  OneOverN   = 1.0 / static_cast<RealType>(source.getTotalNum());
  Kshell     = source.getSimulationCell().getKLists().kshell;
  MaxKshell  = Kshell.size() - 1;
  RhokTot_r.resize(NumK);
  RhokTot_i.resize(NumK);
  values.resize(NumK);
  Kmag.resize(MaxKshell);
  OneOverDnk.resize(MaxKshell);
  for (int ks = 0; ks < MaxKshell; ks++)
  {
    Kmag[ks]       = std::sqrt(source.getSimulationCell().getKLists().ksq[Kshell[ks]]);
    OneOverDnk[ks] = 1.0 / static_cast<RealType>(Kshell[ks + 1] - Kshell[ks]);
  }
  hdf5_out = true;
}

void SkEstimator::resetTargetParticleSet(ParticleSet& P) { sourcePtcl = &P; }

SkEstimator::Return_t SkEstimator::evaluate(ParticleSet& P)
{
  //sum over species
  std::copy(P.getSK().rhok_r[0], P.getSK().rhok_r[0] + NumK, RhokTot_r.begin());
  std::copy(P.getSK().rhok_i[0], P.getSK().rhok_i[0] + NumK, RhokTot_i.begin());
  for (int i = 1; i < NumSpecies; ++i)
    accumulate_elements(P.getSK().rhok_r[i], P.getSK().rhok_r[i] + NumK, RhokTot_r.begin());
  for (int i = 1; i < NumSpecies; ++i)
    accumulate_elements(P.getSK().rhok_i[i], P.getSK().rhok_i[i] + NumK, RhokTot_i.begin());
  if (hdf5_out)
  {
    Vector<RealType>::const_iterator iit_r(RhokTot_r.begin()), iit_r_end(RhokTot_r.end());
    Vector<RealType>::const_iterator iit_i(RhokTot_i.begin()), iit_i_end(RhokTot_i.end());
    for (int i = my_index_; iit_r != iit_r_end; ++iit_r, ++iit_i, ++i)
      P.Collectables[i] += OneOverN * ((*iit_r) * (*iit_r) + (*iit_i) * (*iit_i));
  }
  else
  {
    Vector<RealType>::const_iterator iit_r(RhokTot_r.begin()), iit_r_end(RhokTot_r.end());
    Vector<RealType>::const_iterator iit_i(RhokTot_i.begin()), iit_i_end(RhokTot_i.end());
    for (int i = 0; iit_r != iit_r_end; ++iit_r, ++iit_i, ++i)
      values[i] = OneOverN * ((*iit_r) * (*iit_r) + (*iit_i) * (*iit_i));
  }
  return 0.0;
}

void SkEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (hdf5_out)
  {
    my_index_ = collectables.size();
    std::vector<RealType> tmp(NumK);
    collectables.add(tmp.begin(), tmp.end());
  }
  else
  {
    my_index_ = plist.size();
    for (int i = 0; i < NumK; i++)
    {
      std::stringstream sstr;
      sstr << "sk_" << i;
      int id = plist.add(sstr.str());
    }
  }
}

void SkEstimator::addObservables(PropertySetType& plist)
{
  my_index_ = plist.size();
  for (int i = 0; i < NumK; i++)
  {
    std::stringstream sstr;
    sstr << "sk_" << i;
    int id = plist.add(sstr.str());
  }
}

void SkEstimator::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
    std::copy(values.begin(), values.end(), plist.begin() + my_index_);
}

void SkEstimator::setParticlePropertyList(PropertySetType& plist, int offset)
{
  if (!hdf5_out)
    std::copy(values.begin(), values.end(), plist.begin() + my_index_ + offset);
}


void SkEstimator::registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const
{
  if (hdf5_out)
  {
    std::vector<int> ndim(1, NumK);
    h5desc.emplace_back(name_);
    auto& h5o = h5desc.back();
    h5o.set_dimensions(ndim, my_index_);
    h5o.open(gid);

    hsize_t kdims[2];
    kdims[0]          = NumK;
    kdims[1]          = OHMMS_DIM;
    std::string kpath = name_ + "/kpoints";
    hid_t k_space     = H5Screate_simple(2, kdims, NULL);
    hid_t k_set       = H5Dcreate(gid, kpath.c_str(), H5T_NATIVE_DOUBLE, k_space, H5P_DEFAULT);
    hid_t mem_space   = H5Screate_simple(2, kdims, NULL);
    auto* ptr         = &(sourcePtcl->getSimulationCell().getKLists().kpts_cart[0][0]);
    herr_t ret        = H5Dwrite(k_set, H5T_NATIVE_DOUBLE, mem_space, k_space, H5P_DEFAULT, ptr);
    H5Dclose(k_set);
    H5Sclose(mem_space);
    H5Sclose(k_space);
    H5Fflush(gid, H5F_SCOPE_GLOBAL);
  }
}

bool SkEstimator::put(xmlNodePtr cur)
{
  OhmmsAttributeSet pAttrib;
  std::string hdf5_flag = "no";
  pAttrib.add(hdf5_flag, "hdf5");
  pAttrib.put(cur);
  if (hdf5_flag == "yes")
    hdf5_out = true;
  else
    hdf5_out = false;
  return true;
}

bool SkEstimator::get(std::ostream& os) const { return true; }

std::unique_ptr<OperatorBase> SkEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SkEstimator> myclone = std::make_unique<SkEstimator>(*this);
  myclone->hdf5_out                    = hdf5_out;
  myclone->my_index_                   = my_index_;
  return myclone;
}
} // namespace qmcplusplus
