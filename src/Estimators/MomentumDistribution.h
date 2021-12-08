//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: MomentumEstimator.h
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H
#define QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H
#include <vector>

#include "Configuration.h"
#include "OperatorEstBase.h"
#include "Containers/OhmmsPETE/TinyVector.h"

#include "MomentumDistributionInput.h"


namespace qmcplusplus
{
namespace testing
{
class MomentumDistributionTests;
}
/** Class that collects momentum distribution of electrons
 *  
 */
class MomentumDistribution : public OperatorEstBase
{
public:
  using LatticeType = PtclOnLatticeTraits::ParticleLayout_t;
  using RealType    = QMCTraits::RealType;
  using ComplexType = QMCTraits::ComplexType;
  using ValueType   = QMCTraits::ValueType;
  using PosType     = QMCTraits::PosType;

  //data members set only during construction
  ///input values
  const MomentumDistributionInput input_;
  ///twist angle
  const PosType twist;
  ///lattice vector
  const LatticeType Lattice;
  ///number of samples
  const int M;
  ///normalization factor for n(k)
  const RealType norm_nofK;
  ///list of k-points in Cartesian Coordinates
  std::vector<PosType> kPoints;
  ///weight of k-points (make use of symmetry)
  std::vector<int> kWeights;

  /** @ingroup MomentumDistribution mutable data members
   */
  ///sample positions
  std::vector<PosType> vPos;
  ///wavefunction ratios
  std::vector<ValueType> psi_ratios;
  ///wavefunction ratios all samples
  Matrix<ValueType> psi_ratios_all;
  ///nofK internal
  Vector<RealType> kdotp;
  ///phases
  VectorSoaContainer<RealType, 2> phases;
  ///phases of vPos
  std::vector<VectorSoaContainer<RealType, 2>> phases_vPos;
  ///nofK
  aligned_vector<RealType> nofK;

public:
  /** Constructor for MomentumDistributionInput 
   */
  MomentumDistribution(MomentumDistributionInput&& mdi,
                       size_t np,
                       const PosType& twist,
                       const LatticeType& lattice,
                       DataLocality dl = DataLocality::crowd);

  /** Constructor used when spawing crowd clones
   *  needs to be public so std::make_unique can call it.
   *  Do not use directly unless you've really thought it through.
   */
  MomentumDistribution(const MomentumDistribution& md, DataLocality dl);

  /** This allows us to allocate the necessary data for the DataLocality::queue 
   */
  void startBlock(int steps) override;

  /** standard interface
   */
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;

  /** accumulate 1 or more walkers of MomentumDistribution samples
   */
  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator_t& rng) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of MomentumDistribution's data.
   *
   *  can use base class default until crowd level MomentumDistribution
   *  estimators don't have a copy of the density grid.
   */
  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;

  /** this allows the EstimatorManagerNew to reduce without needing to know the details
   *  of MomentumDistribution's data.
   *
   *  can use base class default until crowd level MomentumDistribution estimators don't have a copy of the density grid.
   */
  //void collect(const OperatorEstBase&  oeb);

  /** this gets us into the hdf5 file
   *
   *  Just parroting for now don't fully understand.
   *, needs to be unraveled and simplified the hdf5 output is another 
   *  big state big coupling design.
   */
  void registerOperatorEstimator(hid_t gid) override;

private:
  MomentumDistribution(const MomentumDistribution& md) = default;

  friend class testing::MomentumDistributionTests;
};

} // namespace qmcplusplus

#endif /* QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H */
