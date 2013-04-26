//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
/** @file MCWalkerConfiguration.h
 * @brief Declaration of a MCWalkerConfiguration
 */
#ifndef QMCPLUSPLUS_BEAD_PARTICLESET_H
#define QMCPLUSPLUS_BEAD_PARTICLESET_H
#include "Particle/ParticleSet.h"
#include "Particle/Walker.h"
#include "Utilities/IteratorUtility.h"
#include "QMCDrivers/MultiChain.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"


namespace qmcplusplus
{

//Forward declaration
struct MultiChain;

/** A set of walkers that are to be advanced by Metropolis Monte Carlo.
 *
 *As a derived class from ParticleSet, MCWalkerConfiguration interacts with
 *QMCHamiltonian and TrialWaveFunction as a ParticleSet, while QMCDrivers
 *use it as multiple walkers whose configurations are advanced according
 to MC algorithms.
 *
 Each walker is represented by Walker<PosVector_t> and
 *MCWalkerConfiguration contains a list of
 *the walkers.  This class enables two possible moves:
 *<ul>
 *<li> move the entire active walkers, similarly to molecu. Suitable for
 *small and big systems with a small time step.
 *<li> move a particle for each walker. Suitable for large systems.

 *</ul>
 */
class Bead_ParticleSet: public ParticleSet
{

public:

  /**enumeration for update*/
  enum {Update_All = 0, ///move all the active walkers
        Update_Walker,  ///move a walker by walker
        Update_Particle ///move a particle by particle
       };

  // Need to typedef Walker_t first
  typedef Bead_ParticleSet::Walker_t Walker_t;
  ///container type of the Properties of a Walker
  typedef Walker_t::PropertyContainer_t  PropertyContainer_t;
  ///container type of Walkers
  typedef std::vector<Walker_t*>         WalkerList_t;
  ///iterator of Walker container
  typedef WalkerList_t::iterator         iterator;
  ///const_iterator of Walker container
  typedef WalkerList_t::const_iterator   const_iterator;
  typedef Walker_t::Buffer_t WBuffer_t;


  typedef Bead_ParticleSet::RealType RealType;
  typedef Bead_ParticleSet::ParticlePos_t ParticlePos_t;
  typedef Bead_ParticleSet::Scalar_t Scalar_t;
  typedef Bead_ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;

  vector<int> BeadSignWgt;
  ParticlePos_t Drift;
  //    ParticlePos_t R;
  vector<ParticleGradient_t*> Gradients;
  vector<ParticleGradient_t*> DriftVectors;
  vector<ParticleLaplacian_t*> Laplacians;
  vector<ParticleGradient_t*> dG_saved;
  vector<ParticleLaplacian_t*> dL_saved;
  Matrix<RealType> Action;
  RealType TransProb[2];


  //saved data!
  ParticlePos_t R_saved;
  vector<int> BeadSignWgt_saved;
  vector<ParticleGradient_t*> Gradients_saved;
  vector<ParticleGradient_t*> DriftVectors_saved;
  vector<ParticleLaplacian_t*> Laplacians_saved;
  Matrix<RealType> Action_saved;
  int nPsi;
  RealType TransProb_saved[2];
  Walker_t::PropertyContainer_t properties_saved;
  ParticlePos_t Drift_saved;

  ///true if using scaled drift
  bool ScaleDrift;
  vector<RealType> Tau_eff;

  void SaveOldData();
  void RestoreOldData();
  void SetGradientAndLaplacian(int ipsi);
  void CopyFromBead(Bead& b,vector<TrialWaveFunction*> &Psi);
  void CopyToBead(Bead& b,vector<TrialWaveFunction*> &Psi);







  ///default constructor
  Bead_ParticleSet();

  ///default constructor: copy only ParticleSet
  Bead_ParticleSet(const ParticleSet& mcw,int nPsi);

  ///default destructor
  ~Bead_ParticleSet();



  ///clean up the walker list and make a new list
  void resize(int numWalkers, int numPtcls);

  ///make random moves for all the walkers
  //void sample(iterator first, iterator last, value_type tauinv);
  ///make a random move for a walker
  void sample(iterator it, RealType tauinv);


  ///return the number of particles per walker
  inline int getParticleNum() const
  {
    return R.size();
  }
  void getDrift(vector<RealType>& LogNorm);

  void getScaledDrift(vector<RealType>& LogNorm,RealType Tau);
  void getScaledDriftSingle(vector<RealType>& LogNorm, RealType Tau, int ipsi);

  inline bool updatePbyP() const
  {
    return ReadyForPbyP;
  }




protected:

  ///true if the buffer is ready for particle-by-particle updates
  bool ReadyForPbyP;

  ///update-mode index
  int UpdateMode;



};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3051 $   $Date: 2008-08-26 13:23:42 -0500 (Tue, 26 Aug 2008) $
 * $Id: MCWalkerConfiguration.h 3051 2008-08-26 18:23:42Z jnkim $
 ***************************************************************************/
