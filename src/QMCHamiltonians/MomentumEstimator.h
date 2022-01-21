//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_MOMENTUM_HAMILTONIAN_H
#define QMCPLUSPLUS_MOMENTUM_HAMILTONIAN_H
#include "QMCHamiltonians/OperatorBase.h"
namespace qmcplusplus
{
class MomentumEstimator : public OperatorBase
{
public:
  MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi);
  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  void addObservables(PropertySetType& plist) {}
  void addObservables(PropertySetType& plist, BufferType& olist) override;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hid_t gid) const override;
  void setObservables(PropertySetType& plist) override;
  void setParticlePropertyList(PropertySetType& plist, int offset) override;
  bool putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode);
  bool put(xmlNodePtr cur) override { return false; };
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final;
  void setRandomGenerator(RandomGenerator* rng) override;
  //resize the internal data by input k-point list
  void resize(const std::vector<PosType>& kin, const int Min);
  ///number of samples
  int M;
  ///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
  ///lattice vector
  const ParticleSet::ParticleLayout_t& lattice_;
  ///normalization factor for n(k)
  RealType norm_nofK;
  ///random generator
  RandomGenerator myRNG;
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
  ///list of k-points in Cartesian Coordinates
  std::vector<PosType> kPoints;
  ///weight of k-points (make use of symmetry)
  std::vector<int> kWeights;
  ///nofK
  aligned_vector<RealType> nofK;
  /// print to hdf5 or scalar.dat
  bool hdf5_out;
  PosType twist;
};

} // namespace qmcplusplus
#endif
