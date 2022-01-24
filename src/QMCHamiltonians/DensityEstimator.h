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


#ifndef QMCPLUSPLUS_DENSITY_HAMILTONIAN_H
#define QMCPLUSPLUS_DENSITY_HAMILTONIAN_H
#include "QMCHamiltonians/OperatorBase.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "LongRange/LRCoulombSingleton.h"
namespace qmcplusplus
{
using LRHandlerType  = LRCoulombSingleton::LRHandlerType;
using GridType       = LRCoulombSingleton::GridType;
using RadFunctorType = LRCoulombSingleton::RadFunctorType;

class DensityEstimator : public OperatorBase
{
public:
  DensityEstimator(ParticleSet& elns);
  int potentialIndex;
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;
  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override;

  void addObservables(PropertySetType& plist) {}
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;

  inline int getGridIndex(int i, int j, int k) const { return my_index_ + k + NumGrids[2] * (j + NumGrids[1] * i); }

  inline int getGridIndexPotential(int i, int j, int k) const
  {
    return potentialIndex + k + NumGrids[2] * (j + NumGrids[1] * i);
  }


private:
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
  ///density
  //Array<RealType,OHMMS_DIM> density, Vavg;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();
};

} // namespace qmcplusplus
#endif
