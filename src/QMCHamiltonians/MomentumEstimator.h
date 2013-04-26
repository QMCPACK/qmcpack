//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(vector<observable_helper*>& h5desc, hid_t gid) const ;
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
  void resize(const vector<PosType>& kin,const vector<RealType>& qin);
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
  vector<ValueType> psi_ratios;
  ///nofK internal
  Vector<RealType> kdotp;
  ///phases
  Vector<ComplexType> phases;
  ///list of k-points in Cartesian Coordinates
  vector<PosType> kPoints;
  ///weight of k-points (make use of symmetry)
  vector<int> kWeights;
  ///dims of a grid for k points
  int kgrid;
  ///nofK
  Vector<RealType> nofK;
  ///list of Q for the Compton profile
  vector<RealType> Q;
  ///compton profile at q
  Vector<RealType> compQ;
  /// print to hdf5 or scalar.dat
  bool hdf5_out;

  vector<vector<int> > mappedQtonofK;
//     vector<vector<int> > mappednofKtoK;
  vector<RealType> mappedQnorms;
  vector<RealType> mappedKnorms;
  PosType twist;
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: MomentumEstimator.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/
