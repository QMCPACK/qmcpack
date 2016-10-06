//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_RMCSINGLE_OMP_H
#define QMCPLUSPLUS_RMCSINGLE_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "Particle/Reptile.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a RMC using threaded execution.
 */
  class RMCSingleOMP:public QMCDriver, public CloneManager
  {
  public:
    /// Constructor.
    typedef PtclAttribTraits::ParticlePos_t ParticlePos_t;
    typedef Reptile::ReptileConfig_t ReptileConfig_t;

      RMCSingleOMP (MCWalkerConfiguration & w, TrialWaveFunction & psi,
		    QMCHamiltonian & h, HamiltonianPool & hpool,
		    WaveFunctionPool & ppool);
    bool run ();
    bool put (xmlNodePtr cur);
    //inline std::vector<RandomGenerator_t*>& getRng() { return Rng;}
  private:
    ///period for walker dump
    int myPeriod4WalkerDump;
    ///option to enable/disable drift equation for RMC
    std::string rescaleDrift;
    ///number of beads on the reptile, beta/tau
    int beads;
    //number of reptiles.  
    int nReptiles;
    ///rescale for time step studies. some int>2 and new beads are inserted in between the old ones.
    int resizeReptile;

    //Calculating the reptiles from scratch or from a previous VMC/DMC/RMC run.
    bool fromScratch;
    ///projection time of reptile
    RealType beta;
//       vector of indices for the action and transprob
      std::vector < int >Action;
      std::vector < int >TransProb;

    int prestepsVMC;
    ///check the run-time environments
    inline void resetVars ()
    {
      prestepsVMC = -1;
      beads = -1;
      beta = -1;
      nReptiles = -1;
    };
    void resetRun ();
    //This will resize the MCWalkerConfiguration and initialize the ReptileList.  Does not care for previous runs.  
    void resetReptiles (int nReptiles, int nbeads, RealType tau);
    //This will resize the MCWalkerConfiguration and initialize Reptile list.  It will then reinitialize the MCWC with a list of Reptile coordinates
    void resetReptiles (std::vector< ReptileConfig_t > &reptile_samps,
			RealType tau);
    //For # of walker samples, create that many reptiles with nbeads each.  Initialize each reptile to have the value of the walker "seed".
    void resetReptiles (std::vector< ParticlePos_t > &walker_samps, int nbeads,
			RealType tau);
    ///copy constructor
      RMCSingleOMP (const RMCSingleOMP & a):QMCDriver (a), CloneManager (a)
    {
    }
    /// Copy operator (disabled).
    RMCSingleOMP & operator= (const RMCSingleOMP &)
    {
      return *this;
    }


  };
}

#endif
