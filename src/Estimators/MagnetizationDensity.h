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

  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;
  size_t getFullDataSize();
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;
  void registerOperatorEstimator(hdf_archive& file) override;

  template<class VAL>
  VAL integrateMagnetizationDensity(const std::vector<VAL>& fgrid)
  {
    Real start=0.0;
    Real stop=2*M_PI;
    Real deltax=(stop-start)/(nsamples_ - 1.0);
    VAL val=0.0;
    switch(integrator_)
    {
      case Integrator::SIMPSONS:
        //There's a normalizing 2pi factor in this integral.  Take care of it here.  
        val=integrateBySimpsonsRule(fgrid,deltax)/Real(2.0*M_PI);
        
        break;
      case Integrator::MONTECARLO:
        //Integrand has form that can be monte carlo sampled.  This means the 2*PI
        //is taken care of if the integrand is uniformly sampled between [0,2PI),
        //which it is.  The following is just an average.  
        val=std::accumulate(fgrid.begin(),fgrid.end(),VAL(0));
        val/=Real(nsamples_);
        break;
    }  

    return val;
  };
  
private:
  MagnetizationDensity(const MagnetizationDensity& magdens) = default;
  void generateSpinIntegrand(ParticleSet& pset, 
                             TrialWaveFunction& wfn, 
                             const int iat, 
                             std::vector<Value>& sx,
                             std::vector<Value>& sy,
                             std::vector<Value>& sz);
 
  template<class VAL>
  VAL integrateBySimpsonsRule(const std::vector<VAL>& fgrid, Real gridDx) const
  {
    VAL sint(0.0);
    int gridsize = fgrid.size();
    for (int is = 1; is < gridsize - 1; is += 2)
    sint += Real(4. / 3.) * gridDx * fgrid[is];

    for (int is = 2; is < gridsize - 1; is += 2)
      sint += Real(2. / 3.) * gridDx * fgrid[is];
   
    sint += Real(1. / 3.) * gridDx * fgrid[0];
    sint += Real(1. / 3.) * gridDx * fgrid[gridsize - 1];
    
    return sint;
  };


  void generateGrid(std::vector<Real>& sgrid) const;
  void generateUniformGrid(std::vector<Real>& sgrid, const Real start, const Real stop) const;

  template<class RAN_GEN> 
  void generateRandomGrid(std::vector<Real>& sgrid, RAN_GEN& rng, Real start, Real stop) const
  {
    size_t npoints=sgrid.size();
    for(int i=0; i<npoints;i++)
      sgrid[i] = (stop-start)*rng();
  };

  
  size_t computeBin(const Position& r, const unsigned int component) const;
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
