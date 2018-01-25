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
#include <QMCHamiltonians/QMCHamiltonianBase.h>
namespace qmcplusplus
{

class MomentumEstimator: public QMCHamiltonianBase
{
public:

  MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi);
  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode);
  bool put(xmlNodePtr cur)
  {
    return false;
  };
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  void setRandomGenerator(RandomGenerator_t* rng);
  //resize the internal data by input k-point list
  void resize(const std::vector<PosType>& kin, const int Min);
  ///number of samples
  int M;
  ///normalization factor for n(k)
  RealType norm_nofK;
  ///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
  ///lattice vector
  ParticleSet::ParticleLayout_t Lattice;
  ///random generator
  RandomGenerator_t myRNG;
  ///sample positions
  std::vector<PosType> vPos;
  ///wavefunction ratios
  std::vector<ValueType> psi_ratios;
  ///wavefunction ratios all samples
  Matrix<ValueType> psi_ratios_all;
  ///nofK internal
  Vector<RealType> kdotp;
  ///phases
  VectorSoaContainer<RealType,2> phases;
  ///phases of vPos
  std::vector<VectorSoaContainer<RealType,2> > phases_vPos;
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

}
#endif

