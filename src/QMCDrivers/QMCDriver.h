//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**
 * @file QMCDriver.h
 * @brief Declaration of QMCDriver
 */
#ifndef QMCPLUSPLUS_QMCDRIVER_H
#define QMCPLUSPLUS_QMCDRIVER_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/PooledData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/ScalarEstimatorManager.h"
#include "Utilities/OhmmsInfo.h"
//#include "QMCDrivers/MolecuFixedNodeBranch.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"

namespace qmcplusplus {

  /** @defgroup QMCDrivers QMC Driver group
   * QMC drivers that implement QMC algorithms
   */

  /** @defgroup WalkerByWalker QMC Drivers using walker-by-walker update
   * @brief Move all the particles for each walker
   */

  /** @defgroup ParticleByParticle QMC Drivers using particle-by-particle update
   * @brief Move particle by particle
   */

  /** @defgroup MultiplePsi QMC Drivers for energy differences
   * @brief Umbrella sampling over multiple H/Psi
   *
   * This class of QMC drivers are suitable to evaluate
   * the energy differences of multiple H-Psi pairs.
   */

  class MCWalkerConfiguration;
  /** @ingroup QMCDrivers
   * @{
   * @brief abstract base class for QMC engines 
   */
  class QMCDriver: public QMCTraits {

  public:

    /** enumeration coupled with QMCMode */
    enum {QMC_UPDATE_MODE, QMC_MULTIPLE, QMC_OPTIMIZE};

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef Walker_t::Buffer_t              Buffer_t;
    typedef SimpleFixedNodeBranch           BranchEngineType;
    //typedef MolecuFixedNodeBranch<RealType> BranchEngineType;

    /// Constructor.
    QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    virtual ~QMCDriver();

    ///return current step
    inline int current() const { return CurrentStep;}

    /** Set the status of the QMCDriver
     * @param aname the root file name
     * @param h5name root name of the master hdf5 file containing previous qmcrun
     * @param append if ture, the run is a continuation of the previous qmc
     *
     * All output files will be of
     * the form "aname.s00X.suffix", where "X" is number
     * of previous QMC runs for the simulation and "suffix"
     * is the suffix for the output file. 
     */
    void setStatus(const string& aname, const string& h5name, bool append);

    /** add QMCHamiltonian/TrialWaveFunction pair for multiple
     * @param h QMCHamiltonian
     * @param psi TrialWaveFunction
     *
     * *Multiple* drivers use multiple H/Psi pairs to perform correlated sampling
     * for energy difference evaluations.
     */
    void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi);

    /** initialize all the walkers
     */
    void initialize();

    void process(xmlNodePtr cur);

    /** return a xmlnode with update **/
    xmlNodePtr getQMCNode();

    void putWalkers(vector<xmlNodePtr>& wset);

    virtual bool run() = 0;
    
    virtual bool put(xmlNodePtr cur) = 0;


  protected:

    ///counts the number of qmc runs
    static int Counter;

    ///branch engine, declared to be static
    static BranchEngineType *branchEngine;

    ///randomize it
    bool ResetRandom;
    ///flag to append or restart the run
    bool AppendRun;

    ///the number of blocks to be rolled back
    int RollBackBlocks;

    /** period of dumping walker configurations and everything else for restart
     *
     * The unit is a block.
     */
    int Period4CheckPoint;

    /** period of recording walker configurations
     *
     * Default is 0 indicating that only the last configuration will be saved.
     */
    int Period4WalkerDump;

    /** bits to classify QMCDriver
     *
     * - QMCDriverMode[0]? particle-by-particle: walker-by-walker
     * - QMCDriverMode[1]? multiple H/Psi : single H/Psi
     * Other bits are unused.
     */
    bitset<4> QMCDriverMode;

    ///current step
    IndexType CurrentStep;

    ///maximum number of blocks
    IndexType nBlocks;

    ///maximum number of steps
    IndexType nSteps;

    ///counter for number of moves accepted
    IndexType nAccept;

    ///counter for number of moves /rejected
    IndexType nReject; 

    ///the number of walkers
    IndexType nTargetWalkers;

    ///timestep
    RealType Tau;

    ///Time-step factor \f$ 1/(2\Tau)\f$
    RealType m_oneover2tau;
    ///Time-step factor \f$ \sqrt{\Tau}\f$
    RealType m_sqrttau;

    ///Observables manager
    ScalarEstimatorManager* Estimators;

    ///pointer to qmc node in xml file
    xmlNodePtr qmcNode;

    ///type of qmc: assigned by subclasses
    string QMCType;
    ///the root of h5File
    string h5FileRoot;
    ///root of all the output files
    string RootName;

    ///store any parameter that has to be read from a file
    ParameterSet m_param;

    ///walker ensemble
    MCWalkerConfiguration& W;

    ///trial function
    TrialWaveFunction& Psi;

    ///Hamiltonian
    QMCHamiltonian& H;

    ///a list of TrialWaveFunctions for multiple method
    vector<TrialWaveFunction*> Psi1;

    ///a list of QMCHamiltonians for multiple method
    vector<QMCHamiltonian*> H1;

    ///a list of mcwalkerset element
    vector<xmlNodePtr> mcwalkerNodePtr;

    ///temporary storage for drift
    ParticleSet::ParticlePos_t drift;

    ///temporary storage for random displacement
    ParticleSet::ParticlePos_t deltaR;

    ///stream for the log file 
    //OhmmsInform *LogOut;

    ///temporary buffer to accumulate data
    //ostrstream log_buffer;

    //PooledData<RealType> HamPool;

    ///Copy Constructor (disabled).
    QMCDriver(const QMCDriver& a): W(a.W), Psi(a.Psi), H(a.H), Estimators(0){}
 
    bool putQMCInfo(xmlNodePtr cur);

    void addWalkers(int nwalkers);  

    void updateWalkers();

    /** record the state of the block
     * @param block current block
     *
     * virtual function with a default implementation
     */
    virtual void recordBlock(int block);

    /** finalize a qmc section
     * @param block current block
     *
     * Accumulate energy and weight is written to a hdf5 file.
     * Finialize the estimators
     */
    bool finalize(int block);

  };
  /**@}*/
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
