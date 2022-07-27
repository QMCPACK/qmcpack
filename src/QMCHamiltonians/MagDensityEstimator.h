//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
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
  enum MD_Integrator
  {
    MD_INT_SIMPSONS,
    MD_INT_TRAP,
    MD_INT_MC
  };

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

  /** Returns flattened array index of grid point i=x, j=y, k=z, and vector component dim.
   * \param[in] i x index.
   * \param[in] j y index.
   * \param[in] k z index.
   * \param[in] dim xyz vector component.
   * \return index in flattened array.
   */
  inline int getMagGridIndex(int i, int j, int k, int dim) const
  {
    return my_index_ + dim + OHMMS_DIM * (k + NumGrids[2] * (j + NumGrids[1] * i));
  }

  /** Returns a uniform grid with nGridPoints between start and stop.
   * \param[in] start start of interval.
   * \param[in] stop end of interval.
   * \param[in] nGridPoints number of grid points.  
   * \return ParticleScalar array of grid points.
   */
  ParticleSet::ParticleScalar generateUniformGrid(RealType start, RealType stop, int nGridPoints);

  /** Returns a random (uniformly distributed) grid with nGridPoints between start and stop.
   * \param[in] start start of interval.
   * \param[in] stop end of interval.
   * \param[in] nGridPoints number of grid points.  
   * \return ParticleScalar array of grid points.
   */
  ParticleSet::ParticleScalar generateRandomGrid(RealType start, RealType stop, int nSamples);

  /** Takes a discretized function fgrid with uniform grid spacing and integrates the function with
   *   Simpson's 3/8 rule.
   * \param[in] fgrid discretized function to integrate
   * \param[in] gridDx uniform grid spacing
   * \return  Value of int f(x).
   */
  RealType integrateBySimpsonsRule(const std::vector<RealType>& fgrid, RealType gridDx);

  /** Takes a discretized function fgrid with uniform grid spacing and integrates the function with
   *   Trapezoidal rule.
   * \param[in] fgrid discretized function to integrate
   * \param[in] gridDx uniform grid spacing
   * \return  Value of int f(x).
   */
  RealType integrateByTrapzRule(const std::vector<RealType>& fgrid, RealType gridDx);

  /** Returns the average of an array fgrid.  Used for MC integration. 
   * \param[in] fgrid discretized function to integrate
   * \return  Average of fgrid.
   */
  RealType average(const std::vector<RealType>& fgrid);

private:
  ///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
  ///true if any direction of a supercell is periodic
  bool Periodic;
  ///normalization factor
  RealType Norm;
  ///number of grids.  +1 is for the additional vector index.
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
