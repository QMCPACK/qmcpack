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
  void resize(const std::vector<PosType>& kin,const std::vector<RealType>& qin);
  ///number of samples
  int M;
  ///normalization factor for n(k)
  RealType norm_nofK;
  ///normalization factor for the Compton profile
  RealType norm_compQ;
  ///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
  ///lattice vector
  ParticleSet::ParticleLayout_t Lattice;
  ///random generator
  RandomGenerator_t myRNG;
  ///wavefunction ratios
  std::vector<ValueType> psi_ratios;
  ///nofK internal
  Vector<RealType> kdotp;
  ///phases
  Vector<ComplexType> phases;
  ///list of k-points in Cartesian Coordinates
  std::vector<PosType> kPoints;
  ///weight of k-points (make use of symmetry)
  std::vector<int> kWeights;
  ///dims of a grid for k points
  int kgrid;
  ///nofK
  Vector<RealType> nofK;
  ///list of Q for the Compton profile
  std::vector<RealType> Q;
  ///compton profile at q
  Vector<RealType> compQ;
  /// print to hdf5 or scalar.dat
  bool hdf5_out;

  std::vector<std::vector<int> > mappedQtonofK;
//     std::vector<std::vector<int> > mappednofKtoK;
  std::vector<RealType> mappedQnorms;
  std::vector<RealType> mappedKnorms;
  PosType twist;
};

}
#endif

