//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_NATURALORBITALS_HAMILTONIAN_H
#define QMCPLUSPLUS_NATURALORBITALS_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <QMCWaveFunctions/SPOSetBase.h>



namespace qmcplusplus
{

class NaturalOrbitals: public QMCHamiltonianBase
{
public:

  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t

  NaturalOrbitals(ParticleSet& elns, TrialWaveFunction& psi);


  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist,BufferType& olist);
  void moveRprime();

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

///reference to the trial wavefunction for ratio evaluations
  TrialWaveFunction& refPsi;
///use a sposet to compute the NO in.
  SPOSetBasePtr& Phi;

/// Stuff for generating new particle positions.
///random generator, timestep, and steps between samples. For generating $r^\prime$
  ParticleSet Rprime;
  RandomGenerator_t myRNG;
  ValueType tau, P_r, P_r_prime;
  int steps;

  ///wavefunction ratios
  std::vector<ValueType> psi_ratios;
  ValueVector_t phiV,temp_phiV;
  ValueMatrix_t phii_phij;
  GradVector_t dphiV;
  ///list of k-points in Cartesian Coordinates
  /// print to hdf5 or scalar.dat
  bool overwrite_hdf5;

};

}
#endif

