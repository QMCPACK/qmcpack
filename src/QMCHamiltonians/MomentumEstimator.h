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
    bool put(xmlNodePtr cur);
    bool get(std::ostream& os) const;
    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
    void setRandomGenerator(RandomGenerator_t* rng);
    //resize the internal data by input k-point list
    void resize(const vector<PosType>& kin);
    ///normalization factor
    RealType NormFactor;
    ///reference to the trial wavefunction for ratio evaluations
    TrialWaveFunction& refPsi;
    ///random generator
    RandomGenerator_t myRNG;
    ///nofK
    Vector<RealType> kdotp;
    ///list of k-points in Cartesian Coordinates
    vector<PosType> kPoints;
    ///nofK
    Vector<RealType> nofK;
    ///phases
    Matrix<ComplexType> phases;
    ///wavefunction ratios
    vector<RealType> psi_ratios;
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: MomentumEstimator.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
