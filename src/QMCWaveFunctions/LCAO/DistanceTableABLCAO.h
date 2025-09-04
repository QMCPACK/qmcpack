#ifndef QMCPLUSPLUS_DISTANCETABLE_AB_LCAO_H
#define QMCPLUSPLUS_DISTANCETABLE_AB_LCAO_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{

class DistanceTableABLCAO
{
public:
  using RealType = OHMMS_PRECISION;
  using PosType = TinyVector<RealType, OHMMS_DIM>;
  
  template<typename DT>
  using OffloadPinnedVector = Vector<DT, OffloadPinnedAllocator<DT>>;

private:
  // Simple resource struct - no inheritance needed for now
  struct LCAODistanceResource
  {
    OffloadPinnedVector<RealType> mw_distances;
    OffloadPinnedVector<RealType> mw_displacements;
  };
  
  LCAODistanceResource mw_resource_;  // Direct member, not pointer
  
  // Problem sizes
  int num_ions_;
  int num_elec_;
  int num_centers_;
  size_t nwalkers_;
  
public:
  DistanceTableABLCAO(const ParticleSet& ions, const ParticleSet& elecs);
  ~DistanceTableABLCAO() = default;
  
  // Evaluation
  void mw_evaluate(
    const RefVectorWithLeader<ParticleSet>& elec_list,
    const ParticleSet& ions,
    int iat,
    int jion);
    
  // Clean accessors
    OffloadPinnedVector<RealType>& getDisplacementsVector()
  { return mw_resource_.mw_displacements; }

  const OffloadPinnedVector<RealType>& getDisplacementsVector() const
  { return mw_resource_.mw_displacements; }


  const RealType* getDistancesDevice() const 
  { return mw_resource_.mw_distances.device_data(); }
  
  const RealType* getDisplacementsDevice() const 
  { return mw_resource_.mw_displacements.device_data(); }
    // Add method to evaluate all centers
  void mw_evaluate_all_centers(
    const RefVectorWithLeader<ParticleSet>& elec_list,
    const ParticleSet& ions,
    int iat,
    int num_centers);

};

} // namespace qmcplusplus
#endif
