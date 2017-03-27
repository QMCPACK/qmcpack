//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign  
//
//
// File created by:  Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign  
//////////////////////////////////////////////////////////////////////////////////////


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

  std::vector<int> BeadSignWgt;
  ParticlePos_t Drift;
  //    ParticlePos_t R;
  std::vector<ParticleGradient_t*> Gradients;
  std::vector<ParticleGradient_t*> DriftVectors;
  std::vector<ParticleLaplacian_t*> Laplacians;
  std::vector<ParticleGradient_t*> dG_saved;
  std::vector<ParticleLaplacian_t*> dL_saved;
  Matrix<RealType> Action;
  RealType TransProb[2];


  //saved data!
  ParticlePos_t R_saved;
  std::vector<int> BeadSignWgt_saved;
  std::vector<ParticleGradient_t*> Gradients_saved;
  std::vector<ParticleGradient_t*> DriftVectors_saved;
  std::vector<ParticleLaplacian_t*> Laplacians_saved;
  Matrix<RealType> Action_saved;
  int nPsi;
  RealType TransProb_saved[2];
  Walker_t::PropertyContainer_t properties_saved;
  ParticlePos_t Drift_saved;

  ///true if using scaled drift
  bool ScaleDrift;
  std::vector<RealType> Tau_eff;

  void SaveOldData();
  void RestoreOldData();
  void SetGradientAndLaplacian(int ipsi);
  void CopyFromBead(Bead& b,std::vector<TrialWaveFunction*> &Psi);
  void CopyToBead(Bead& b,std::vector<TrialWaveFunction*> &Psi);







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
  void getDrift(std::vector<RealType>& LogNorm);

  void getScaledDrift(std::vector<RealType>& LogNorm,RealType Tau);
  void getScaledDriftSingle(std::vector<RealType>& LogNorm, RealType Tau, int ipsi);

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
