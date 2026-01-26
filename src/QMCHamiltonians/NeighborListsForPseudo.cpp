//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "NeighborListsForPseudo.h"

namespace qmcplusplus
{
NeighborListsForPseudo::NeighborListsForPseudo(size_t num_elecs,
                                               size_t num_ions,
                                               const std::vector<NonLocalECPComponent*>& pp)
    : elec_neighbor_ions_(num_elecs), ion_neighbor_elecs_(num_ions), PP(pp)
{}

void NeighborListsForPseudo::markAffectedElecs(const DistanceTableAB& myTable,
                                               int iel,
                                               std::vector<bool>& elecTMAffected)
{
  std::vector<int>& neighbor_ions = elec_neighbor_ions_.getNeighborList(iel);
  for (int iat = 0; iat < PP.size(); iat++)
  {
    if (PP[iat] == nullptr)
      continue;
    auto old_distance = myTable.getDistRow(iel)[iat];
    auto new_distance = myTable.getTempDists()[iat];
    bool moved        = false;
    // move out
    if (old_distance < PP[iat]->getRmax() && new_distance >= PP[iat]->getRmax())
    {
      moved                = true;
      auto& neighbor_elecs = ion_neighbor_elecs_.getNeighborList(iat);
      auto iter_at         = std::find(neighbor_ions.begin(), neighbor_ions.end(), iat);
      auto iter_el         = std::find(neighbor_elecs.begin(), neighbor_elecs.end(), iel);
      if (iter_at == neighbor_ions.end())
        throw std::runtime_error("BUG! NeighborListsForPseudo::markAffectedElecs atom not found in the neighbor list.");
      if (iter_el == neighbor_elecs.end())
        throw std::runtime_error(
            "BUG! NeighborListsForPseudo::markAffectedElecs electron not found in the neighbor list.");
      *iter_at = neighbor_ions.back();
      *iter_el = neighbor_elecs.back();
      neighbor_ions.pop_back();
      neighbor_elecs.pop_back();
      elecTMAffected[iel] = true;
    }
    // move in
    if (old_distance >= PP[iat]->getRmax() && new_distance < PP[iat]->getRmax())
    {
      moved                = true;
      auto& neighbor_elecs = ion_neighbor_elecs_.getNeighborList(iat);
      neighbor_elecs.push_back(iel);
      neighbor_ions.push_back(iat);
    }
    // move around
    if (moved || (old_distance < PP[iat]->getRmax() && new_distance < PP[iat]->getRmax()))
    {
      auto& neighbor_elecs = ion_neighbor_elecs_.getNeighborList(iat);
      for (int jel = 0; jel < neighbor_elecs.size(); ++jel)
        elecTMAffected[neighbor_elecs[jel]] = true;
    }
  }
}

void NeighborListsForPseudo::clear()
{
  for (int iat = 0; iat < ion_neighbor_elecs_.size(); iat++)
    ion_neighbor_elecs_.getNeighborList(iat).clear();
  for (int jel = 0; jel < elec_neighbor_ions_.size(); jel++)
    elec_neighbor_ions_.getNeighborList(jel).clear();
}

void NeighborListsForPseudo::addElecIonPair(int jel, int iat)
{
  elec_neighbor_ions_.getNeighborList(jel).push_back(iat);
  ion_neighbor_elecs_.getNeighborList(iat).push_back(jel);
}
} // namespace qmcplusplus
