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
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
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
#include "Estimators/EstimatorManager.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"

class Communicate;

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
  class HDFWalkerOutput;
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

    /// Constructor.
    QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    virtual ~QMCDriver();

    ///return current step
    inline int current() const { return CurrentStep;}

    /** set the update mode
     * @param pbyp if true, use particle-by-particle update
     */
    inline void setUpdateMode(bool pbyp) 
    {
      QMCDriverMode[QMC_UPDATE_MODE]=pbyp;
    }

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

    /** set qmcComm 
     * @param c Communicate*
     */
    void setCommunicator(Communicate* c=0);
    /** add QMCHamiltonian/TrialWaveFunction pair for multiple
     * @param h QMCHamiltonian
     * @param psi TrialWaveFunction
     *
     * *Multiple* drivers use multiple H/Psi pairs to perform correlated sampling
     * for energy difference evaluations.
     */
    void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi);

    /** initialize with xmlNode
     */
    void process(xmlNodePtr cur);

    /** return a xmlnode with update **/
    xmlNodePtr getQMCNode();

    void putWalkers(vector<xmlNodePtr>& wset);

    virtual bool run() = 0;
    
    virtual bool put(xmlNodePtr cur) = 0;

    inline string getEngineName() const 
    { 
      return  QMCType;
    }

    template<class PDT>
    void setValue(const string& aname, PDT x) {
      m_param.setValue(aname,x);
    }

    ///set the BranchEngineType
    void setBranchEngine(BranchEngineType* be) {
      branchEngine=be;
    }

    ///return BranchEngineType*
    BranchEngineType* getBranchEngine() {
      return branchEngine;
    }

    int addObservable(const string& aname) {
      if(Estimators)
        return Estimators->addObservable(aname.c_str());
      else
        return -1;
    }

    RealType getObservable(int i) {
      return Estimators->getObservable(i);
    }

    ///Observables manager
    EstimatorManager* Estimators;

  protected:

    ///branch engine
    BranchEngineType *branchEngine;
    ///randomize it
    bool ResetRandom;
    ///flag to append or restart the run
    bool AppendRun;
    /** the number of times this QMCDriver is executed
     *
     * MyCounter is initialized to zero by the constructor and is incremented 
     * whenever a run is completed by calling finalize(int block) or 
     * using MyCounter++ as in RQMC.
     */
    int MyCounter;
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

    ///maximum cpu in secs
    RealType MaxCPUSecs;

    ///Time-step factor \f$ 1/(2\Tau)\f$
    RealType m_oneover2tau;
    ///Time-step factor \f$ \sqrt{\Tau}\f$
    RealType m_sqrttau;

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

    ///communicator
    Communicate* qmcComm;
    ///record engine for walkers
    HDFWalkerOutput* wOut;
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

    //void updateWalkers();

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
 * $RCSfile: QMCDriver.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
