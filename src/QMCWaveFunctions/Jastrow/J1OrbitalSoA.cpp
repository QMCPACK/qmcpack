//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-


#include "J1OrbitalSoA.h"
#include "SoaDistanceTableABOMPTarget.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{

template<typename T>
struct J1OrbitalSoAMultiWalkerMem : public Resource
{
  // fused buffer for fast transfer
  Vector<char, OffloadPinnedAllocator<char>> transfer_buffer;
  // multi walker result
  Vector<T, OffloadPinnedAllocator<T>> mw_vals;
  // multi walker -1
  Vector<int, OffloadPinnedAllocator<int>> mw_minus_one;

  void resize_minus_one(size_t size)
  {
    if (mw_minus_one.size() < size)
    {
      mw_minus_one.resize(size, -1);
      mw_minus_one.updateTo();
    }
  }

  J1OrbitalSoAMultiWalkerMem() : Resource("J1OrbitalSoAMultiWalkerMem") {}

  J1OrbitalSoAMultiWalkerMem(const J1OrbitalSoAMultiWalkerMem&) : J1OrbitalSoAMultiWalkerMem() {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<J1OrbitalSoAMultiWalkerMem>(*this); }
};

template<typename FT>
J1OrbitalSoA<FT>::J1OrbitalSoA(const std::string& obj_name, const ParticleSet& ions, ParticleSet& els, bool use_offload)
    : WaveFunctionComponent(obj_name),
      use_offload_(use_offload),
      myTableID(els.addTable(ions, use_offload ? DTModes::ALL_OFF : DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
      Nions(ions.getTotalNum()),
      Nelec(els.getTotalNum()),
      NumGroups(ions.groups()),
      Ions(ions)
{
  if (my_name_.empty())
    throw std::runtime_error("J1OrbitalSoA object name cannot be empty!");

  if (use_offload_)
    assert(ions.getCoordinates().getKind() == DynamicCoordinateKind::DC_POS_OFFLOAD);

  initialize(els);

  // set up grp_ids
  grp_ids.resize(Nions);
  int count = 0;
  for (int ig = 0; ig < NumGroups; ig++)
    for (int j = ions.first(ig); j < ions.last(ig); j++)
      grp_ids[count++] = ig;
  assert(count == Nions);
  grp_ids.updateTo();
}

template<typename FT>
J1OrbitalSoA<FT>::~J1OrbitalSoA() = default;

template<typename FT>
void J1OrbitalSoA<FT>::checkSanity() const
{
  if (std::any_of(J1Functors.begin(), J1Functors.end(), [](auto* ptr) { return ptr == nullptr; }))
    app_warning() << "One-body Jastrow \"" << my_name_ << "\" doesn't cover all the particle pairs. "
                  << "Consider fusing multiple entries if they are of the same type for optimal code performance."
                  << std::endl;
}

template<typename FT>
void J1OrbitalSoA<FT>::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<J1OrbitalSoAMultiWalkerMem<RealType>>());
}

template<typename FT>
void J1OrbitalSoA<FT>::acquireResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<J1OrbitalSoA<FT>>();
  wfc_leader.mw_mem_handle_ = collection.lendResource<J1OrbitalSoAMultiWalkerMem<RealType>>();
}

template<typename FT>
void J1OrbitalSoA<FT>::releaseResource(ResourceCollection& collection,
                                       const RefVectorWithLeader<WaveFunctionComponent>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<J1OrbitalSoA<FT>>();
  collection.takebackResource(wfc_leader.mw_mem_handle_);
}

template<typename FT>
void J1OrbitalSoA<FT>::mw_evaluateRatios(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                         std::vector<std::vector<ValueType>>& ratios) const
{
  if (!use_offload_)
  {
    WaveFunctionComponent::mw_evaluateRatios(wfc_list, vp_list, ratios);
    return;
  }

  // add early return to prevent from accessing vp_list[0]
  if (wfc_list.size() == 0)
    return;
  auto& wfc_leader        = wfc_list.getCastedLeader<J1OrbitalSoA<FT>>();
  auto& vp_leader         = vp_list.getLeader();
  const auto& mw_refPctls = vp_leader.getMultiWalkerRefPctls();
  auto& mw_mem            = wfc_leader.mw_mem_handle_.getResource();
  auto& mw_vals           = mw_mem.mw_vals;
  auto& mw_minus_one      = mw_mem.mw_minus_one;
  const int nw            = wfc_list.size();

  const size_t nVPs = mw_refPctls.size();
  mw_vals.resize(nVPs);
  mw_mem.resize_minus_one(nVPs);

  const auto& dt_leader(vp_leader.getDistTableAB(wfc_leader.myTableID));

  FT::mw_evaluateV(NumGroups, GroupFunctors.data(), wfc_leader.Nions, grp_ids.data(), nVPs, mw_minus_one.data(),
                   dt_leader.getMultiWalkerDataPtr(), dt_leader.getPerTargetPctlStrideSize(), mw_vals.data(),
                   mw_mem.transfer_buffer);

  size_t ivp = 0;
  for (int iw = 0; iw < nw; ++iw)
  {
    const VirtualParticleSet& vp = vp_list[iw];
    auto& wfc                    = wfc_list.getCastedElement<J1OrbitalSoA<FT>>(iw);
    for (int k = 0; k < vp.getTotalNum(); ++k, ivp++)
      ratios[iw][k] = std::exp(wfc.Vat[mw_refPctls[ivp]] - mw_vals[ivp]);
  }
  assert(ivp == nVPs);
}

template class J1OrbitalSoA<BsplineFunctor<QMCTraits::RealType>>;
template class J1OrbitalSoA<
    CubicSplineSingle<QMCTraits::RealType, CubicBspline<QMCTraits::RealType, LINEAR_1DGRID, FIRSTDERIV_CONSTRAINTS>>>;
template class J1OrbitalSoA<UserFunctor<QMCTraits::RealType>>;
template class J1OrbitalSoA<ShortRangeCuspFunctor<QMCTraits::RealType>>;
template class J1OrbitalSoA<PadeFunctor<QMCTraits::RealType>>;
template class J1OrbitalSoA<Pade2ndOrderFunctor<QMCTraits::RealType>>;
} // namespace qmcplusplus
