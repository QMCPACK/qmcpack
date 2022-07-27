//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MAGDENSITY_HAMILTONIAN_H
#define QMCPLUSPLUS_MAGDENSITY_HAMILTONIAN_H
#include "QMCHamiltonians/OperatorBase.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "LongRange/LRCoulombSingleton.h"
namespace qmcplusplus
{
using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
using GridType       = LRCoulombSingleton::GridType;
using RadFunctorType = LRCoulombSingleton::RadFunctorType;

class MagDensityEstimator : public OperatorBase
{
public:
  enum MD_Integrator{MD_INT_SIMPSONS, MD_INT_TRAP, MD_INT_MC};

  MagDensityEstimator(ParticleSet& elns, TrialWaveFunction& psi);
  
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;
  // We'll let the base class throw an error if this gets touched.
  //void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override;

  void addObservables(PropertySetType& plist) {}
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  // index of grid point i=x, j=y, k=z, and vector component dim.
  inline int getMagGridIndex(int i, int j, int k, int dim) const {return my_index_ + dim + OHMMS_DIM*(k + NumGrids[2] * (j+NumGrids[1] * i));}

  
  ParticleSet::ParticleScalar generateUniformGrid(RealType start, RealType stop, int nGridPoints);
  ParticleSet::ParticleScalar generateRandomGrid(RealType start, RealType stop, int nSamples);
  RealType integrateBySimpsonsRule(const std::vector<RealType>& fgrid, RealType gridDx);
  RealType integrateByTrapzRule(const std::vector<RealType>& fgrid, RealType gridDx);
  RealType average(const std::vector<RealType>& fgrid);

private:
  ///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
  ///true if any direction of a supercell is periodic
  bool Periodic;
  ///normalization factor
  RealType Norm;
  ///number of grids
  TinyVector<int, OHMMS_DIM + 1> NumGrids;
  ///bin size
  TinyVector<RealType, OHMMS_DIM> Delta;
  ///inverse
  TinyVector<RealType, OHMMS_DIM> DeltaInv;
  ///scaling factor for conversion
  TinyVector<RealType, OHMMS_DIM> ScaleFactor;
  ///lower bound
  TinyVector<RealType, OHMMS_DIM> density_min;
  ///upper bound
  TinyVector<RealType, OHMMS_DIM> density_max;
  ///name of the density data
  
  std::string prefix;
  //Number of MC samples or grid points to perform spin integration.
  int nSamples_;
  //Enum to type of integration to be performed for spin integral.
  MD_Integrator integrator_;  

  void resize();

};

} // namespace qmcplusplus
#endif
