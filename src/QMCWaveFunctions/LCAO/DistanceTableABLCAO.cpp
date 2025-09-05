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
  : DistanceTable(ions, elecs, DTModes::ALL_OFF),
    num_ions_(ions.getTotalNum()),
    num_elec_(elecs.getTotalNum())
	{
		ion_cached_.resize(3 * num_ions_);
  for (int i = 0; i < num_ions_; ++i)
    for (int d = 0; d < 3; ++d)
      ion_cached_[3*i + d] = ions.R[i][d];
  ion_cached_.updateTo();
	}

// Required pure virtual method implementations (stubs for LCAO usage)
void DistanceTableABLCAO::evaluate(ParticleSet& P) {}
void DistanceTableABLCAO::move(const ParticleSet& P, const PosType& rnew, const IndexType iat, bool prepare_old) {}
void DistanceTableABLCAO::update(IndexType jat) {}
int DistanceTableABLCAO::get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const { return -1; }

/** Compute displacements and lattice translation vectors for LCAO evaluation
 * @param elec_list multi-walker electron particle sets
 * @param ions ion particle set
 * @param iat target electron index
 * @param num_centers number of ion centers
 * @param displ_list_tr output displacement vectors [dim + 3*(walker + center*nw)]
 * @param Tv_list output lattice translation vectors (for PBC)
 * 
 * This function computes ion-electron displacements entirely on GPU using the
 * fused new position buffer. For PBC systems, it patches the results with
 * minimum image convention and lattice translations.
 */

void DistanceTableABLCAO::mw_evaluate(
    const RefVectorWithLeader<ParticleSet>& elec_list,
    const ParticleSet& ions,
    int iat,
    int num_centers,
    OffloadPinnedVector<RealType>& displ_list_tr,
    OffloadPinnedVector<RealType>& Tv_list)
{
 const size_t nw = elec_list.size();
  
  const auto* coords_leader = dynamic_cast<const RealSpacePositionsOMPTarget*>(
      &elec_list.getLeader().getCoordinates());
  if (!coords_leader)
    throw std::runtime_error("Electron coordinates not in OMPTarget format!");
  
  const RealType* ion_ptr = ion_cached_.device_data();  // Use cached ions
  RealType* disp_out = displ_list_tr.device_data();
  RealType* tv_out = Tv_list.device_data();
  
  const auto& fused_new_pos = coords_leader->getFusedNewPosBuffer();
  const RealType* elec_device_ptr = fused_new_pos.device_data();
  const size_t elec_stride = fused_new_pos.capacity();
  
  const auto& lattice = ions.getLattice();
  const bool is_pbc = (lattice.SuperCellEnum != SUPERCELL_OPEN);
  
  // GPU kernel for all cases
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(3) \
                  is_device_ptr(elec_device_ptr, ion_ptr, disp_out, tv_out)")
  for (int iw = 0; iw < static_cast<int>(nw); ++iw)
    for (int c = 0; c < num_centers; ++c)
      for (int d = 0; d < 3; ++d)
      {
        const RealType e_pos = elec_device_ptr[iw + d * elec_stride];
        const RealType ion_pos = ion_ptr[c * 3 + d];
        
        size_t out_idx = d + 3 * (iw + c * nw);
        disp_out[out_idx] = ion_pos - e_pos;
        tv_out[out_idx] = 0.0;
      } 
  // Apply PBC corrections on host if needed
  if (is_pbc) {
  // Need to get data back to apply minimum image
  displ_list_tr.updateFrom();
  
  for (size_t iw = 0; iw < nw; ++iw) {
    for (int c = 0; c < num_centers; ++c) {
      PosType disp;
      for (int d = 0; d < 3; ++d) {
        size_t idx = d + 3 * (iw + c * nw);
        disp[d] = displ_list_tr[idx];
      }
      
      // applyMinimumImage modifies disp in-place
      lattice.applyMinimumImage(disp);
      
      for (int d = 0; d < 3; ++d) {
        size_t idx = d + 3 * (iw + c * nw);
        Tv_list[idx] = displ_list_tr[idx] - disp[d];  // Original - wrapped
        displ_list_tr[idx] = disp[d];  // Store wrapped displacement
      }
    }
  }
  
  // Send corrected values back
  displ_list_tr.updateTo();
#if defined(QMC_COMPLEX)
  Tv_list.updateTo();
#endif
  }
}

} // namespace qmcplusplus
