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

  std::string getClassName() const override;
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;


private:
  ///true if any direction of a supercell is periodic
  bool periodic_;
  ///number of grids
  TinyVector<int, OHMMS_DIM + 1> num_grids_;
  ///bin size
  TinyVector<RealType, OHMMS_DIM> delta_;
  ///inverse
  TinyVector<RealType, OHMMS_DIM> delta_inv_;
  ///scaling factor for conversion
  TinyVector<RealType, OHMMS_DIM> scale_factor_;
  ///lower bound
  TinyVector<RealType, OHMMS_DIM> density_min_;
  ///upper bound
  TinyVector<RealType, OHMMS_DIM> density_max_;
  ///name of the density data
  std::string prefix;

  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();

  /**
   * @brief Get the linearized grid Index object from 3D coordinates
   * @param i  x-index
   * @param j  y-index
   * @param k  k-index
   * @return int linearized index
   */
  int getGridIndex(int i, int j, int k) const noexcept;
};

} // namespace qmcplusplus
#endif
