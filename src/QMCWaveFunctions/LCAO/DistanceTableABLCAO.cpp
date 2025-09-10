//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Anouar Benali, abenali.sci@hotmail.com, Qubit Pharmaceuticals.
//
// File created by:  Anouar Benali, abenali.sci@hotmail.com, Qubit Pharmaceuticals.
//////////////////////////////////////////////////////////////////////////////////////

#include "DistanceTableABLCAO.h"
#include "Platforms/OMPTarget/OMPTargetMath.hpp"
#include "Particle/RealSpacePositionsOMPTarget.h"

namespace qmcplusplus
{

DistanceTableABLCAO::DistanceTableABLCAO(const ParticleSet& ions, const ParticleSet& elecs)
    : DistanceTable(ions, elecs, DTModes::ALL_OFF), num_ions_(ions.getTotalNum()), num_elec_(elecs.getTotalNum())
{}

// Required pure virtual method implementations (stubs for LCAO usage)
void DistanceTableABLCAO::evaluate(ParticleSet& P) {}
void DistanceTableABLCAO::move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) {}
void DistanceTableABLCAO::update(IndexType jat) {}
int DistanceTableABLCAO::get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const { return -1; }

/** GPU-optimized displacement computation for non-PBC molecular systems
 * @param elec_list multi-walker electron particle sets
 * @param ions ion particle set
 * @param iat active electron index
 * @param num_centers number of ion centers to evaluate
 * @param displ_list_tr [out] displacement vectors (ion - electron) in SoA layout [3*(iw + c*nw) + d]
 * @param Tv_list [out] lattice translation vectors (set to 0 for non-PBC)
 * 
 * Computes ion-electron displacements entirely on GPU using the fused new position buffer
 * that contains active electron positions for all walkers. This avoids CPU-GPU transfers
 * by leveraging data already resident on device memory.
 * 
 * Memory layout for output arrays:
 * - Both arrays have size 3 * num_centers * num_walkers
 * - Element [d + 3*(iw + c*nw)] contains dimension d of displacement from electron in walker iw to center c
 */

void DistanceTableABLCAO::mw_evaluate(const RefVectorWithLeader<ParticleSet>& elec_list,
                                      const ParticleSet& ions,
                                      int iat,
                                      int num_centers,
                                      OffloadPinnedVector<RealType>& displ_list_tr,
                                      OffloadPinnedVector<RealType>& Tv_list)
{
  const size_t nw = elec_list.size();

  // Ion coordinates on device
  const auto& ion_coords = dynamic_cast<const RealSpacePositionsOMPTarget&>(ions.getCoordinates());
  const RealType* ion_ptr = ion_coords.getDevicePtr();
  const size_t ion_stride = ion_coords.getAllParticlePos().capacity();

  // Get the fused buffer that's ALREADY ON DEVICE - no transfer needed!
  const auto& coords_leader = dynamic_cast<const RealSpacePositionsOMPTarget&>(
      elec_list.getLeader().getCoordinates());
  const auto& fused_new_pos = coords_leader.getFusedNewPosBuffer();
  const RealType* elec_ptr = fused_new_pos.device_data();  // Already on device!
  const size_t elec_stride = fused_new_pos.capacity();

  RealType* disp_out = displ_list_tr.device_data();
  RealType* tv_out = Tv_list.device_data();

  // Direct GPU kernel - no transfers!
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(3) \
                  is_device_ptr(elec_ptr, ion_ptr, disp_out, tv_out)")
  for (int iw = 0; iw < static_cast<int>(nw); ++iw)
    for (int c = 0; c < num_centers; ++c)
      for (int d = 0; d < 3; ++d)
      {
        // Fused buffer layout: [dim][walker]
        const RealType e_pos = elec_ptr[d * elec_stride + iw];
        const RealType ion_pos = ion_ptr[d * ion_stride + c];
        size_t out_idx = d + 3 * (iw + c * nw);

        disp_out[out_idx] = ion_pos - e_pos;
        tv_out[out_idx] = 0.0;
      }
}




void DistanceTableABLCAO::mw_evaluate_pbc(const RefVectorWithLeader<ParticleSet>& elec_list,
                                          const ParticleSet& ions,
                                          int iat,
                                          int num_centers,
                                          OffloadPinnedVector<RealType>& displ_list_tr,
                                          OffloadPinnedVector<RealType>& Tv_list,
                                          int table_index)
{
  const size_t nw = elec_list.size();
  
  auto* Tv_host = Tv_list.data();
  auto* displ_host = displ_list_tr.data();
  
  for (size_t iw = 0; iw < nw; iw++) {
    const auto& coordR = elec_list[iw].activeR(iat);
    const auto& d_table = elec_list[iw].getDistTableAB(table_index);
    const auto& displ = (elec_list[iw].getActivePtcl() == iat) ? 
                         d_table.getTempDispls() : d_table.getDisplRow(iat);
    
    for (int c = 0; c < num_centers; c++) {
      for (size_t idim = 0; idim < 3; idim++) {
        size_t idx = idim + 3 * (iw + c * nw);
        Tv_host[idx] = (ions.R[c][idim] - coordR[idim]) - displ[c][idim];
        displ_host[idx] = displ[c][idim];
      }
    }
  }
  
  // Transfer to device
  Tv_list.updateTo();
  displ_list_tr.updateTo();
}
} // namespace qmcplusplus
