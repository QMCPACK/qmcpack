//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MAGNETIZATION_DENSITY_H
#define QMCPLUSPLUS_MAGNETIZATION_DENSITY_H

#include <vector>
#include <functional>

#include "Estimators/OperatorEstBase.h"
#include "type_traits/complex_help.hpp"
#include "ParticleBase/RandomSeqGenerator.h"
#include <SpeciesSet.h>
#include <StdRandom.h>
#include "Estimators/MagnetizationDensityInput.h"

namespace qmcplusplus
{

namespace testing
{
class MagnetizationDensityTests;
}

class MagnetizationDensity : public OperatorEstBase
{
public:
  using Value         = QMCTraits::ValueType;
  using FullPrecValue = QMCTraits::FullPrecValueType;
  using Real          = RealAlias<Value>;
  using FullPrecReal  = RealAlias<FullPrecValue>;
  using Grad          = TinyVector<Value, OHMMS_DIM>;
  using Lattice       = PtclOnLatticeTraits::ParticleLayout;
  using Position      = QMCTraits::PosType;
  using Integrator    = MagnetizationDensityInput::Integrator;
  static constexpr int DIM = QMCTraits::DIM;

  MagnetizationDensity(MagnetizationDensityInput&& min, const Lattice& lattice);
  MagnetizationDensity(const MagnetizationDensity& magdens, DataLocality dl);

  void startBlock(int steps) override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomGenerator& rng) override;

  size_t getFullDataSize();
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;
  void registerOperatorEstimator(hdf_archive& file) override;

  template<class VAL>
  VAL integrateBySimpsonsRule(const std::vector<VAL>& fgrid, Real gridDx);

  template<class VAL>
  VAL integrateMagnetizationDensity(const std::vector<VAL>& fgrid);
private:
  MagnetizationDensity(const MagnetizationDensity& magdens) = default;
  void generateSpinIntegrand(ParticleSet& pset, 
                             TrialWaveFunction& wfn, 
                             const int iat, 
                             std::vector<Value>& sx,
                             std::vector<Value>& sy,
                             std::vector<Value>& sz);
 
  void generateGrid(std::vector<Real>& sgrid);
  void generateUniformGrid(std::vector<Real>& sgrid, const Real start, const Real stop);

  template<class RAN_GEN> 
  void generateRandomGrid(std::vector<Real>& sgrid, RAN_GEN& rng, Real start, Real stop);
  
  size_t computeBin(const Position& r, const unsigned int component);
  MagnetizationDensityInput input_;
  Integrator integrator_;
  Lattice lattice_;
  Position rcorner_;
  Position center_;
  TinyVector<int, DIM> grid_;
  TinyVector<int, DIM> gdims_;
  size_t npoints_;
  int nsamples_;

  friend class testing::MagnetizationDensityTests;
};
} // namespace qmcplusplus
#endif  /* QMCPLUSPLUS_MAGNETIZATION_DENSITY_H */
