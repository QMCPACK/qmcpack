#include "DistanceTableABLCAO.h"
#include "Platforms/OMPTarget/OMPTargetMath.hpp"

namespace qmcplusplus
{

DistanceTableABLCAO::DistanceTableABLCAO(
    const ParticleSet& ions, 
    const ParticleSet& elecs)
  : num_ions_(ions.getTotalNum()),
    num_elec_(elecs.getTotalNum()),
    num_centers_(0),
    nwalkers_(0)
{
  app_log() << "\n==============================================\n";
  app_log() << "LCAO Distance Table Created Successfully!\n";
  app_log() << "  Number of ions: " << num_ions_ << "\n";
  app_log() << "  Number of electrons: " << num_elec_ << "\n";
  app_log() << "==============================================\n" << std::endl;
}

void DistanceTableABLCAO::mw_evaluate(
    const RefVectorWithLeader<ParticleSet>& elec_list,
    const ParticleSet& ions,
    int iat,
    int jion)
{
  const size_t nw = elec_list.size();
  nwalkers_ = nw;
  num_centers_ = 1;
  
  // Resize storage
  mw_resource_.mw_distances.resize(nw);
  mw_resource_.mw_displacements.resize(3 * nw);
  
  // References for cleaner code
  auto& distances = mw_resource_.mw_distances;
  auto& displacements = mw_resource_.mw_displacements;
  
  // Prepare electron positions
  OffloadPinnedVector<RealType> elec_positions(3 * nw);
  for (size_t iw = 0; iw < nw; ++iw)
  {
    const auto& P = elec_list[iw];
    PosType elec_pos = (P.getActivePtcl() == iat) ? P.activeR(iat) : P.R[iat];
    for (int d = 0; d < 3; ++d)
      elec_positions[3*iw + d] = elec_pos[d];
  }
  
  // Transfer to GPU
  elec_positions.updateTo();
  
  // GPU kernel
  RealType* dist_ptr = distances.device_data();
  RealType* disp_ptr = displacements.device_data(); 
  const RealType* elec_ptr = elec_positions.device_data();
  const RealType ion_x = ions.R[jion][0];
  const RealType ion_y = ions.R[jion][1];
  const RealType ion_z = ions.R[jion][2];
  
  PRAGMA_OFFLOAD("omp target teams distribute parallel for \
                  is_device_ptr(dist_ptr, disp_ptr, elec_ptr)")
  for (int iw = 0; iw < static_cast<int>(nw); ++iw)
  {
    RealType dx = elec_ptr[3*iw + 0] - ion_x;
    RealType dy = elec_ptr[3*iw + 1] - ion_y;
    RealType dz = elec_ptr[3*iw + 2] - ion_z;
    
    dist_ptr[iw] = std::sqrt(dx*dx + dy*dy + dz*dz);
    disp_ptr[3*iw + 0] = dx;
    disp_ptr[3*iw + 1] = dy;
    disp_ptr[3*iw + 2] = dz;
  }
 static int call_count = 0;
  if (++call_count <= 3) {
    app_log() << ">>> GPU computed " << nw 
              << " distances for ion " << jion << "\n";
  }
}  
void DistanceTableABLCAO::mw_evaluate_all_centers(
    const RefVectorWithLeader<ParticleSet>& elec_list,
    const ParticleSet& ions,
    int iat,
    int num_centers)
{
  const size_t nw = elec_list.size();
  nwalkers_ = nw;
  num_centers_ = num_centers;
  
  // Resize for all centers
  mw_resource_.mw_distances.resize(num_centers * nw);
  mw_resource_.mw_displacements.resize(3 * num_centers * nw);
  
  // Prepare electron AND ion positions
  OffloadPinnedVector<RealType> elec_positions(3 * nw);
  OffloadPinnedVector<RealType> ion_positions(3 * num_centers);
  
  // Fill electron positions
  for (size_t iw = 0; iw < nw; ++iw)
  {
    const auto& P = elec_list[iw];
    PosType elec_pos = (P.getActivePtcl() == iat) ? P.activeR(iat) : P.R[iat];
    for (int d = 0; d < 3; ++d)
      elec_positions[3*iw + d] = elec_pos[d];
  }
  
  // Fill ion positions (FIX for the warning)
  for (int ic = 0; ic < num_centers; ++ic)
    for (int d = 0; d < 3; ++d)
      ion_positions[3*ic + d] = ions.R[ic][d];
  
  // Transfer both to GPU
  elec_positions.updateTo();
  ion_positions.updateTo();
  
  // GPU kernel
  RealType* dist_ptr = mw_resource_.mw_distances.device_data();
  RealType* disp_ptr = mw_resource_.mw_displacements.device_data();
  const RealType* elec_ptr = elec_positions.device_data();
  const RealType* ion_ptr = ion_positions.device_data();
  
  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
                  is_device_ptr(dist_ptr, disp_ptr, elec_ptr, ion_ptr)")
  for (int ic = 0; ic < num_centers; ++ic)
    for (int iw = 0; iw < static_cast<int>(nw); ++iw)
    {
      // Now use the copied ion positions
      RealType ion_x = ion_ptr[3*ic + 0];
      RealType ion_y = ion_ptr[3*ic + 1];
      RealType ion_z = ion_ptr[3*ic + 2];
      
      RealType elec_x = elec_ptr[3*iw + 0];
      RealType elec_y = elec_ptr[3*iw + 1];
      RealType elec_z = elec_ptr[3*iw + 2];
      
      // FIX: Use correct sign convention (ion - electron)
      RealType dx = ion_x - elec_x;  // Changed sign!
      RealType dy = ion_y - elec_y;  // Changed sign!
      RealType dz = ion_z - elec_z;  // Changed sign!
   
      // Store in layout [center * nw + walker]
      int idx = ic * nw + iw;
      dist_ptr[idx] = std::sqrt(dx*dx + dy*dy + dz*dz);
      disp_ptr[3*idx + 0] = dx;
      disp_ptr[3*idx + 1] = dy;
      disp_ptr[3*idx + 2] = dz;
    }
  
  static int call_count = 0;
  if (++call_count <= 3) {
    app_log() << ">>> GPU computed distances for " 
              << num_centers << " centers, " 
              << nw << " walkers\n";
  }
}
} // namespace qmcplusplus
