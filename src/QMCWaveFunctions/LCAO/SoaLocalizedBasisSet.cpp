//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                   Anouar Benali, abenali@gmail.com, Qubit Pharmaceuticals
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#include <memory>
#include "SoaLocalizedBasisSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/RealSpacePositionsOMPTarget.h"
#include "SoaAtomicBasisSet.h"
#include "MultiQuinticSpline1D.h"
#include "MultiFunctorAdapter.h"
#include "Numerics/SoaCartesianTensor.h"
#include "Numerics/SoaSphericalTensor.h"

namespace qmcplusplus
{
template<class COT, typename ORBT>
struct SoaLocalizedBasisSet<COT, ORBT>::SoaLocalizedBSetMultiWalkerMem : public Resource
{
  SoaLocalizedBSetMultiWalkerMem() : Resource("SoaLocalizedBasisSet") {}

  SoaLocalizedBSetMultiWalkerMem(const SoaLocalizedBSetMultiWalkerMem&) : SoaLocalizedBSetMultiWalkerMem() {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<SoaLocalizedBSetMultiWalkerMem>(*this);
  }
  std::unique_ptr<DistanceTableABLCAO> lcao_distance_table;
  Vector<RealType, OffloadPinnedAllocator<RealType>> Tv_list;
  Vector<RealType, OffloadPinnedAllocator<RealType>> displ_list_tr;
  Vector<size_t, OffloadPinnedAllocator<size_t>> walker_offsets;
  Vector<int, OffloadPinnedAllocator<int>> active_particles;
  Vector<int, OffloadPinnedAllocator<int>> target_counts;
  Vector<RealType, OffloadPinnedAllocator<RealType>> ion_positions;
  Vector<RealType, OffloadPinnedAllocator<RealType>> electron_positions;
};

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::createResource(ResourceCollection& collection) const
{
  auto resource = std::make_unique<SoaLocalizedBSetMultiWalkerMem>();
  collection.addResource(std::make_unique<SoaLocalizedBSetMultiWalkerMem>());
  for (int i = 0; i < LOBasisSet.size(); i++)
    LOBasisSet[i]->createResource(collection);
}
template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::acquireResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basisset_list) const
{
  auto& loc_basis_leader = basisset_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();
  assert(this == &loc_basis_leader);
  loc_basis_leader.mw_mem_handle_ = collection.lendResource<SoaLocalizedBSetMultiWalkerMem>();

  auto& basisset_leader = loc_basis_leader.LOBasisSet;
  for (int i = 0; i < basisset_leader.size(); i++)
  {
    const auto one_species_basis_list(extractOneSpeciesBasisRefList(basisset_list, i));
    basisset_leader[i]->acquireResource(collection, one_species_basis_list);
  }
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::releaseResource(
    ResourceCollection& collection,
    const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basisset_list) const
{
  auto& loc_basis_leader = basisset_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();
  assert(this == &loc_basis_leader);

  collection.takebackResource(loc_basis_leader.mw_mem_handle_);

  auto& basisset_leader = loc_basis_leader.LOBasisSet;
  for (int i = 0; i < basisset_leader.size(); i++)
  {
    const auto one_species_basis_list(extractOneSpeciesBasisRefList(basisset_list, i));
    basisset_leader[i]->releaseResource(collection, one_species_basis_list);
  }
}
template<class COT, typename ORBT>
RefVectorWithLeader<COT> SoaLocalizedBasisSet<COT, ORBT>::extractOneSpeciesBasisRefList(
    const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basisset_list,
    int id)
{
  auto& loc_basis_leader = basisset_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();
  RefVectorWithLeader<COT> one_species_basis_list(*loc_basis_leader.LOBasisSet[id]);
  one_species_basis_list.reserve(basisset_list.size());
  for (size_t iw = 0; iw < basisset_list.size(); iw++)
    one_species_basis_list.push_back(
        *basisset_list.template getCastedElement<SoaLocalizedBasisSet<COT, ORBT>>(iw).LOBasisSet[id]);
  return one_species_basis_list;
}


template<class COT, typename ORBT>
SoaLocalizedBasisSet<COT, ORBT>::SoaLocalizedBasisSet(ParticleSet& ions, ParticleSet& els)
    : ions_(ions),
      myTableIndex(els.addTable(ions, DTModes::NEED_FULL_TABLE_ANYTIME | DTModes::NEED_VP_FULL_TABLE_ON_HOST)),
      SuperTwist(0.0),
      NumCenter_timer_(createGlobalTimer("SoaLocalizedBasisSet::mw_evaluateVGL_Numcenter", timer_level_fine))
{
  NumCenters = ions.getTotalNum();
  NumTargets = els.getTotalNum();
  LOBasisSet.resize(ions.getSpeciesSet().getTotalNum());
  BasisOffset.resize(NumCenters + 1);
  BasisSetSize = 0;
  initializeSpeciesOffsets();
}

template<class COT, typename ORBT>
SoaLocalizedBasisSet<COT, ORBT>::SoaLocalizedBasisSet(const SoaLocalizedBasisSet& a)
    : SoaBasisSetBase<ORBT>(a),
      NumCenters(a.NumCenters),
      NumTargets(a.NumTargets),
      ions_(a.ions_),
      myTableIndex(a.myTableIndex),
      SuperTwist(a.SuperTwist),
      BasisOffset(a.BasisOffset),
      NumCenter_timer_(createGlobalTimer("SoaLocalizedBasisSet::mw_evaluateVGL_Numcenter", timer_level_fine))
{
  LOBasisSet.reserve(a.LOBasisSet.size());
  for (auto& elem : a.LOBasisSet)
    LOBasisSet.push_back(std::make_unique<COT>(*elem));
  initializeSpeciesOffsets();
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::setPBCParams(
    const TinyVector<int, 3>& PBCImages,
    const TinyVector<double, 3> Sup_Twist,
    const Vector<ValueType, OffloadPinnedAllocator<ValueType>>& phase_factor,
    const Array<RealType, 2, OffloadPinnedAllocator<RealType>>& pbc_displacements)
{
  for (int i = 0; i < LOBasisSet.size(); ++i)
    LOBasisSet[i]->setPBCParams(PBCImages, Sup_Twist, phase_factor, pbc_displacements);

  SuperTwist = Sup_Twist;
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::setBasisSetSize(int nbs)
{
  const auto& IonID(ions_.GroupID);
  if (BasisSetSize > 0 && nbs == BasisSetSize)
    return;

  if (auto& mapping = ions_.get_map_storage_to_input(); mapping.empty())
  {
    //evaluate the total basis dimension and offset for each center
    BasisOffset[0] = 0;
    for (int c = 0; c < NumCenters; c++)
      BasisOffset[c + 1] = BasisOffset[c] + LOBasisSet[IonID[c]]->getBasisSetSize();
    BasisSetSize = BasisOffset[NumCenters];
  }
  else
  {
    // when particles are reordered due to grouping, AOs need to restore the input order to match MOs.
    std::vector<int> map_input_to_storage(mapping.size());
    for (int c = 0; c < NumCenters; c++)
      map_input_to_storage[mapping[c]] = c;

    std::vector<size_t> basis_offset_input_order(BasisOffset.size(), 0);
    for (int c = 0; c < NumCenters; c++)
      basis_offset_input_order[c + 1] =
          basis_offset_input_order[c] + LOBasisSet[IonID[map_input_to_storage[c]]]->getBasisSetSize();

    for (int c = 0; c < NumCenters; c++)
      BasisOffset[c] = basis_offset_input_order[mapping[c]];

    BasisSetSize = basis_offset_input_order[NumCenters];
  }
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::queryOrbitalsForSType(const std::vector<bool>& corrCenter,
                                                            std::vector<bool>& is_s_orbital) const
{
  const auto& IonID(ions_.GroupID);
  for (int c = 0; c < NumCenters; c++)
  {
    int idx = BasisOffset[c];
    int bss = LOBasisSet[IonID[c]]->getBasisSetSize();
    std::vector<bool> local_is_s_orbital(bss);
    LOBasisSet[IonID[c]]->queryOrbitalsForSType(local_is_s_orbital);
    for (int k = 0; k < bss; k++)
      is_s_orbital[idx++] = corrCenter[c] ? local_is_s_orbital[k] : false;
  }
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateVGL(const ParticleSet& P, int iat, vgl_type& vgl)
{
  const auto& IonID(ions_.GroupID);
  const auto& coordR  = P.activeR(iat);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

  PosType Tv;
  for (int c = 0; c < NumCenters; c++)
  {
    Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
    Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
    Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
    LOBasisSet[IonID[c]]->evaluateVGL(P.getLattice(), dist[c], displ[c], BasisOffset[c], vgl, Tv);
  }
}


template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::mw_evaluateVGL(
    const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basis_list,
    const RefVectorWithLeader<ParticleSet>& P_list,
    int iat,
    OffloadMWVGLArray& vgl_v)
{
  assert(this == &basis_list.getLeader());
  auto& basis_leader = basis_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();

  const int Nw = static_cast<int>(P_list.size());
  assert(vgl_v.size(0) == 5);
  assert(vgl_v.size(1) == static_cast<size_t>(Nw));
  assert(vgl_v.size(2) == BasisSetSize);

#ifndef NDEBUG
  for (int iw = 0; iw < Nw; ++iw)
    assert(P_list[iw].getActivePtcl() == iat);
#endif

  // Scratch buffers (same as before)
  auto& Tv_list       = basis_leader.mw_mem_handle_.getResource().Tv_list;        // pinned host + device
  auto& displ_list_tr = basis_leader.mw_mem_handle_.getResource().displ_list_tr;  // pinned host + device
  Tv_list.resize(3ULL * NumCenters * Nw);
  displ_list_tr.resize(3ULL * NumCenters * Nw);

  auto& lcao_dt = basis_leader.mw_mem_handle_.getResource().lcao_distance_table;
  if (!lcao_dt)
    lcao_dt = std::make_unique<DistanceTableABLCAO>(ions_, P_list.getLeader());

  const auto& lattice = ions_.getLattice();
  if (lattice.SuperCellEnum != SUPERCELL_OPEN)
	  lcao_dt->mw_evaluate_pbc(P_list, ions_, iat, NumCenters, displ_list_tr, Tv_list,myTableIndex);
  else
	  lcao_dt->mw_evaluate(P_list, ions_, iat, NumCenters, displ_list_tr, Tv_list);

  ScopedTimer NumCenter_Wrapper(NumCenter_timer_);
  const auto& species_names = ions_.getSpeciesSet().speciesName;
  const int num_species     = species_names.size();

  for (int s = 0; s < num_species; ++s)
    if (const auto& c_list = species_centers_[s]; c_list.size() > 0)
    {
      const auto& offs = species_center_coffsets_[s];
      auto basis_refs  = extractOneSpeciesBasisRefList(basis_list, s);

      LOBasisSet[s]->mw_evaluateVGL_multiCenter(basis_refs, P_list.getLeader().getLattice(),
                                                vgl_v, displ_list_tr, Tv_list,
                                                Nw, BasisSetSize, c_list, offs, NumCenters);
    }
}


template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateVGH(const ParticleSet& P, int iat, vgh_type& vgh)
{
  const auto& IonID(ions_.GroupID);
  const auto& coordR  = P.activeR(iat);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
  PosType Tv;
  for (int c = 0; c < NumCenters; c++)
  {
    Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
    Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
    Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
    LOBasisSet[IonID[c]]->evaluateVGH(P.getLattice(), dist[c], displ[c], BasisOffset[c], vgh, Tv);
  }
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateVGHGH(const ParticleSet& P, int iat, vghgh_type& vghgh)
{
  // APP_ABORT("SoaLocalizedBasisSet::evaluateVGH() not implemented\n");

  const auto& IonID(ions_.GroupID);
  const auto& coordR  = P.activeR(iat);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
  PosType Tv;
  for (int c = 0; c < NumCenters; c++)
  {
    Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
    Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
    Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
    LOBasisSet[IonID[c]]->evaluateVGHGH(P.getLattice(), dist[c], displ[c], BasisOffset[c], vghgh, Tv);
  }
}



template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::mw_evaluateValueVPs(
    const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basis_list,
    const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
    OffloadMWVArray& vp_basis_v)
{
  assert(this == &basis_list.getLeader());
  auto& basis_leader = basis_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();

  const size_t nVPs = vp_basis_v.size(0);
  assert(vp_basis_v.size(1) == BasisSetSize);

  auto& vps_leader = vp_list.getLeader();

  // Resource scratch
  auto& Tv_list       = basis_leader.mw_mem_handle_.getResource().Tv_list;        // host/device pinned
  auto& displ_list_tr = basis_leader.mw_mem_handle_.getResource().displ_list_tr;  // host/device pinned
  Tv_list.resize(3ULL * NumCenters * nVPs);
  displ_list_tr.resize(3ULL * NumCenters * nVPs);

  RealType* __restrict__ disp_out = displ_list_tr.device_data();
  RealType* __restrict__ tv_out   = Tv_list.device_data();

  // VP DT device buffer: [slabs] with per-target stride
  const auto& dt_leader  = vps_leader.getDistTableAB(myTableIndex);
  const RealType* dt_dev = dt_leader.getMultiWalkerDataDevicePtr();
  const size_t    stride = dt_leader.getPerTargetPctlStrideSize();
  const size_t    padded = getAlignedSize<RealType>(NumCenters);
  if (!dt_dev)
    throw std::runtime_error("mw_evaluateValueVPs: VP DT device buffer is null.");

  // --- 1) Pack displacements (device->device) ---
  // slabs: 0=dist, 1=dx, 2=dy, 3=dz
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                  is_device_ptr(dt_dev, disp_out)")
  for (int ivp = 0; ivp < static_cast<int>(nVPs); ++ivp)
    for (int c = 0; c < NumCenters; ++c)
    {
      const size_t baseDT = static_cast<size_t>(ivp) * stride;
      const RealType dx   = dt_dev[baseDT + 1 * padded + c];
      const RealType dy   = dt_dev[baseDT + 2 * padded + c];
      const RealType dz   = dt_dev[baseDT + 3 * padded + c];

      const size_t baseOut = 3ULL * (static_cast<size_t>(ivp) + static_cast<size_t>(c) * nVPs);
      disp_out[baseOut + 0] = dx;
      disp_out[baseOut + 1] = dy;
      disp_out[baseOut + 2] = dz;
    }

  // --- 2) Tv handling ---
  const auto& lat = vps_leader.getLattice();
  if (lat.SuperCellEnum == SUPERCELL_OPEN)
  {
    // OPEN cell: Tv = 0 (device-only)
    PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                    is_device_ptr(tv_out)")
    for (size_t i = 0; i < 3ULL * NumCenters * nVPs; ++i)
      tv_out[i] = RealType(0);
  }
  else
  {
    // PBC: compute Tv on host from host coords; one upload
    // host VP coords flattened in countVPs order
    const auto vp_host = vps_leader.extractVPCoords(vp_list); // size == nVPs
    // lattice matrices on host
    RealType G[9], Rm[9];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      { G[3*i + j] = lat.G(i, j); Rm[3*i + j] = lat.R(i, j); }

    for (size_t ivp = 0; ivp < nVPs; ++ivp)
    {
      const RealType vx = vp_host[ivp][0];
      const RealType vy = vp_host[ivp][1];
      const RealType vz = vp_host[ivp][2];

      for (int c = 0; c < NumCenters; ++c)
      {
        const RealType ix = ions_.R[c][0];
        const RealType iy = ions_.R[c][1];
        const RealType iz = ions_.R[c][2];

        const RealType drx = ix - vx;
        const RealType dry = iy - vy;
        const RealType drz = iz - vz;

        const RealType sx = G[0]*drx + G[1]*dry + G[2]*drz;
        const RealType sy = G[3]*drx + G[4]*dry + G[5]*drz;
        const RealType sz = G[6]*drx + G[7]*dry + G[8]*drz;

        const RealType nx = std::nearbyint(sx);
        const RealType ny = std::nearbyint(sy);
        const RealType nz = std::nearbyint(sz);

        // Tv = R * n
        const RealType tx = Rm[0]*nx + Rm[1]*ny + Rm[2]*nz;
        const RealType ty = Rm[3]*nx + Rm[4]*ny + Rm[5]*nz;
        const RealType tz = Rm[6]*nx + Rm[7]*ny + Rm[8]*nz;

        const size_t baseOut = 3ULL * (ivp + static_cast<size_t>(c) * nVPs);
        // write Tv into the pinned host view of Tv_list
        Tv_list[baseOut + 0] = tx;
        Tv_list[baseOut + 1] = ty;
        Tv_list[baseOut + 2] = tz;
      }
    }

    // single host->device upload for Tv_list
    Tv_list.updateTo();
  }

  // --- 3) Species-batched AO evaluation (unchanged) ---
  const auto& species_names = ions_.getSpeciesSet().speciesName;
  const int num_species     = species_names.size();
  for (int s = 0; s < num_species; ++s)
    if (auto& c_list_s = species_centers_[s]; c_list_s.size() > 0)
    {
      auto& offs_list_s = species_center_coffsets_[s];
      auto basis_refs   = extractOneSpeciesBasisRefList(basis_list, s);

      LOBasisSet[s]->mw_evaluateV_multiCenter(basis_refs, vps_leader.getLattice(), vp_basis_v,
                                              displ_list_tr, Tv_list,
                                              nVPs, BasisSetSize, c_list_s, offs_list_s, NumCenters);
    }
}
/*

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::mw_evaluateValueVPs(const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basis_list,
                                                          const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                                          OffloadMWVArray& vp_basis_v)
{
  assert(this == &basis_list.getLeader());
  auto& basis_leader = basis_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();

  const size_t nVPs = vp_basis_v.size(0);
  assert(vp_basis_v.size(1) == BasisSetSize); // shape [nVPs, BasisSetSize]
  const auto& IonID(ions_.GroupID);

  auto& vps_leader = vp_list.getLeader();

  const auto dt_list(vps_leader.extractDTRefList(vp_list, myTableIndex));
  const auto coordR_list(vps_leader.extractVPCoords(vp_list));

  // GPU arrays [NumCenters * nVPs, 3]
  // Each center c, vp index iVP => c*nVPs + iVP
  auto& Tv_list       = basis_leader.mw_mem_handle_.getResource().Tv_list;
  auto& displ_list_tr = basis_leader.mw_mem_handle_.getResource().displ_list_tr;
  Tv_list.resize(3ULL * NumCenters * nVPs);
  displ_list_tr.resize(3ULL * NumCenters * nVPs);

  auto* Tv_host    = Tv_list.data();
  auto* displ_host = displ_list_tr.data();

  // "index" is loop over each virtual point
  // i.e. we do "for each (iw, iat) => iVP" in [0..nVPs-1]
  size_t indexVP = 0;
  for (size_t iw = 0; iw < vp_list.size(); iw++)
  {
    // vp_list[iw].getTotalNum() = # of virtual points in this VPS
    int nVP_local = vp_list[iw].getTotalNum();
    for (int iat = 0; iat < nVP_local; iat++)
    {
      const auto& displ = dt_list[iw].getDisplRow(iat);
      // coords for this virtual point = coordR_list[indexVP]
      const auto& vpcoord = coordR_list[indexVP];

      // fill for all centers c in [0..NumCenters-1]
      for (int c = 0; c < NumCenters; c++)
      {
        for (int dim = 0; dim < 3; dim++)
        {
          size_t idx = dim + 3ULL * (indexVP + c * (size_t)nVPs);
          // Ion position is ions_.R[c]
          RealType val    = ions_.R[c][dim] - vpcoord[dim] - displ[c][dim];
          Tv_host[idx]    = val;
          displ_host[idx] = displ[c][dim];
        }
      }
      indexVP++;
    }
  }

#if defined(QMC_COMPLEX)
  Tv_list.updateTo();
#endif
  displ_list_tr.updateTo();

  const auto& species_names = ions_.getSpeciesSet().speciesName;
  const int num_species_    = species_names.size();
  std::vector<std::vector<size_t>> local_centers(num_species_);
  std::vector<std::vector<size_t>> local_offsets(num_species_);
  for (int c = 0; c < NumCenters; c++)
  {
    int s_id = IonID[c];
    local_centers[s_id].push_back(c);
    local_offsets[s_id].push_back(BasisOffset[c]);
  }

  for (int s = 0; s < num_species_; s++)
    if (auto& c_list_s = species_centers_[s]; c_list_s.size() > 0)
    {
      auto& offs_list_s = species_center_coffsets_[s];
      auto basis_refs   = extractOneSpeciesBasisRefList(basis_list, s);

      LOBasisSet[s]->mw_evaluateV_multiCenter(basis_refs, vps_leader.getLattice(), vp_basis_v, displ_list_tr, Tv_list,
                                              nVPs, BasisSetSize, c_list_s, offs_list_s, NumCenters);
    }
}
*/




template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateV(const ParticleSet& P, int iat, ORBT* restrict vals)
{
  const auto& IonID(ions_.GroupID);
  const auto& coordR  = P.activeR(iat);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

  PosType Tv;
  for (int c = 0; c < NumCenters; c++)
  {
    Tv[0] = (ions_.R[c][0] - coordR[0]) - displ[c][0];
    Tv[1] = (ions_.R[c][1] - coordR[1]) - displ[c][1];
    Tv[2] = (ions_.R[c][2] - coordR[2]) - displ[c][2];
    LOBasisSet[IonID[c]]->evaluateV(P.getLattice(), dist[c], displ[c], vals + BasisOffset[c], Tv);
  }
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::mw_evaluateValue(const RefVectorWithLeader<SoaBasisSetBase<ORBT>>& basis_list,
                                                       const RefVectorWithLeader<ParticleSet>& P_list,
                                                       int iat,
                                                       OffloadMWVArray& vals)
{
  assert(this == &basis_list.getLeader());
  auto& basis_leader = basis_list.template getCastedLeader<SoaLocalizedBasisSet<COT, ORBT>>();
  const auto& IonID(ions_.GroupID);
  auto& pset_leader = P_list.getLeader();

  size_t Nw = P_list.size();
  assert(vals.size(0) == Nw);
  assert(vals.size(1) == BasisSetSize);

  auto& Tv_list       = basis_leader.mw_mem_handle_.getResource().Tv_list;
  auto& displ_list_tr = basis_leader.mw_mem_handle_.getResource().displ_list_tr;
  Tv_list.resize(3 * NumCenters * Nw);
  displ_list_tr.resize(3 * NumCenters * Nw);

  for (size_t iw = 0; iw < P_list.size(); iw++)
  {
    const auto& coordR  = P_list[iw].activeR(iat);
    const auto& d_table = P_list[iw].getDistTableAB(myTableIndex);
    const auto& displ   = (P_list[iw].getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

    for (int c = 0; c < NumCenters; c++)
      for (size_t idim = 0; idim < 3; idim++)
      {
        Tv_list[idim + 3 * (iw + c * Nw)]       = (ions_.R[c][idim] - coordR[idim]) - displ[c][idim];
        displ_list_tr[idim + 3 * (iw + c * Nw)] = displ[c][idim];
      }
  }
#if defined(QMC_COMPLEX)
  Tv_list.updateTo();
#endif
  displ_list_tr.updateTo();

  for (int c = 0; c < NumCenters; c++)
  {
    auto one_species_basis_list = extractOneSpeciesBasisRefList(basis_list, IonID[c]);
    LOBasisSet[IonID[c]]->mw_evaluateV(one_species_basis_list, pset_leader.getLattice(), vals, displ_list_tr, Tv_list,
                                       Nw, BasisSetSize, c, BasisOffset[c], NumCenters);
  }
}


template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateGradSourceV(const ParticleSet& P,
                                                          int iat,
                                                          const ParticleSet& ions,
                                                          int jion,
                                                          vgl_type& vgl)
{
  //We need to zero out the temporary array vgl.
  auto* restrict gx = vgl.data(1);
  auto* restrict gy = vgl.data(2);
  auto* restrict gz = vgl.data(3);

  for (int ib = 0; ib < BasisSetSize; ib++)
  {
    gx[ib] = 0;
    gy[ib] = 0;
    gz[ib] = 0;
  }

  const auto& IonID(ions_.GroupID);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

  const auto& coordR = P.activeR(iat);

  PosType Tv;
  Tv[0] = (ions_.R[jion][0] - coordR[0]) - displ[jion][0];
  Tv[1] = (ions_.R[jion][1] - coordR[1]) - displ[jion][1];
  Tv[2] = (ions_.R[jion][2] - coordR[2]) - displ[jion][2];
  //PosType Tv;
  //Tv[0] = Tv[1] = Tv[2] = 0;
  //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
  //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
  //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.
  LOBasisSet[IonID[jion]]->evaluateVGL(P.getLattice(), dist[jion], displ[jion], BasisOffset[jion], vgl, Tv);
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::evaluateGradSourceVGL(const ParticleSet& P,
                                                            int iat,
                                                            const ParticleSet& ions,
                                                            int jion,
                                                            vghgh_type& vghgh)
{
  //We need to zero out the temporary array vghgh.
  auto* restrict gx = vghgh.data(1);
  auto* restrict gy = vghgh.data(2);
  auto* restrict gz = vghgh.data(3);

  auto* restrict hxx = vghgh.data(4);
  auto* restrict hxy = vghgh.data(5);
  auto* restrict hxz = vghgh.data(6);
  auto* restrict hyy = vghgh.data(7);
  auto* restrict hyz = vghgh.data(8);
  auto* restrict hzz = vghgh.data(9);

  auto* restrict gxxx = vghgh.data(10);
  auto* restrict gxxy = vghgh.data(11);
  auto* restrict gxxz = vghgh.data(12);
  auto* restrict gxyy = vghgh.data(13);
  auto* restrict gxyz = vghgh.data(14);
  auto* restrict gxzz = vghgh.data(15);
  auto* restrict gyyy = vghgh.data(16);
  auto* restrict gyyz = vghgh.data(17);
  auto* restrict gyzz = vghgh.data(18);
  auto* restrict gzzz = vghgh.data(19);


  for (int ib = 0; ib < BasisSetSize; ib++)
  {
    gx[ib] = 0;
    gy[ib] = 0;
    gz[ib] = 0;

    hxx[ib] = 0;
    hxy[ib] = 0;
    hxz[ib] = 0;
    hyy[ib] = 0;
    hyz[ib] = 0;
    hzz[ib] = 0;

    gxxx[ib] = 0;
    gxxy[ib] = 0;
    gxxz[ib] = 0;
    gxyy[ib] = 0;
    gxyz[ib] = 0;
    gxzz[ib] = 0;
    gyyy[ib] = 0;
    gyyz[ib] = 0;
    gyzz[ib] = 0;
    gzzz[ib] = 0;
  }

  // Since jion is indexed on the source ions not the ions_ the distinction between
  // ions_ and ions is extremely important.
  const auto& IonID(ions.GroupID);
  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.getActivePtcl() == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.getActivePtcl() == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);

  //Since LCAO's are written only in terms of (r-R), ionic derivatives only exist for the atomic center
  //that we wish to take derivatives of.  Moreover, we can obtain an ion derivative by multiplying an electron
  //derivative by -1.0.  Handling this sign is left to LCAOrbitalSet.  For now, just note this is the electron VGL function.

  const auto& coordR = P.activeR(iat);

  PosType Tv;
  Tv[0] = (ions_.R[jion][0] - coordR[0]) - displ[jion][0];
  Tv[1] = (ions_.R[jion][1] - coordR[1]) - displ[jion][1];
  Tv[2] = (ions_.R[jion][2] - coordR[2]) - displ[jion][2];
  LOBasisSet[IonID[jion]]->evaluateVGHGH(P.getLattice(), dist[jion], displ[jion], BasisOffset[jion], vghgh, Tv);
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::add(int icenter, std::unique_ptr<COT> aos)
{
  LOBasisSet[icenter] = std::move(aos);
}

template<class COT, typename ORBT>
void SoaLocalizedBasisSet<COT, ORBT>::initializeSpeciesOffsets()
{
  const auto& species_names = ions_.getSpeciesSet().speciesName;
  const size_t num_species  = species_names.size();
  const auto& IonID(ions_.GroupID);

  std::vector<std::vector<size_t>> local_centers(num_species);
  std::vector<std::vector<size_t>> local_offsets(num_species);
  for (int c = 0; c < NumCenters; ++c)
  {
    const int s = IonID[c];
    local_centers[s].push_back(c);
    local_offsets[s].push_back(BasisOffset[c]);
  }

  species_centers_.resize(num_species);
  species_center_coffsets_.resize(num_species);

  for (int s = 0; s < num_species; ++s)
  {
    auto& c_list = species_centers_[s];
    auto& o_list = species_center_coffsets_[s];
    c_list.resize(local_centers[s].size());
    o_list.resize(local_offsets[s].size());
    for (size_t i = 0; i < c_list.size(); ++i)
      c_list[i] = local_centers[s][i];
    for (size_t i = 0; i < o_list.size(); ++i)
      o_list[i] = local_offsets[s][i];
    c_list.updateTo();
    o_list.updateTo();
  }
}

template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiQuinticSpline1D<QMCTraits::RealType>, SoaCartesianTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiQuinticSpline1D<QMCTraits::RealType>, SoaSphericalTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<QMCTraits::RealType>>, SoaCartesianTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<QMCTraits::RealType>>, SoaSphericalTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<QMCTraits::RealType>>, SoaCartesianTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
template class SoaLocalizedBasisSet<
    SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<QMCTraits::RealType>>, SoaSphericalTensor<QMCTraits::RealType>>,
    QMCTraits::ValueType>;
} // namespace qmcplusplus
