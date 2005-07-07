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
#ifndef OHMMS_QMC_QMCDRIVER_H
#define OHMMS_QMC_QMCDRIVER_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Utilities/PooledData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Estimators/ScalarEstimatorManager.h"
#include <strstream>
class OhmmsInform;

namespace ohmmsqmc {

  class MCWalkerConfiguration;

  /** Base class to perform QMC simulations. */
  class QMCDriver: public QMCTraits {

  public:

    typedef MCWalkerConfiguration::Walker_t Walker_t;

    /// Constructor.
    QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

    virtual ~QMCDriver();

    void setFileRoot(const string& aname);

    void add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi) {
      H1.push_back(h);
      Psi1.push_back(psi);
    }

    void initialize();

    void process(xmlNodePtr cur);

    virtual bool run() = 0;
    
    virtual bool put(xmlNodePtr cur) = 0;

  protected:

    ///counts the number of qmc runs
    static int Counter;

    ///flag to print walker ensemble
    bool pStride;

    ///Index of the Acceptance Ratio
    int AcceptIndex;

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

    ///timestep to assign Walker::R at the start. Default = 0.0
    RealType FirstStep;

    ///reference energy, only used in DMC
    RealType e_ref;

    ///pointer to qmc node in xml file
    xmlNodePtr qmcNode;

    ParameterSet m_param;

    ///type of qmc: assigned by subclasses
    string QMCType;

    ///root of all the output files
    string RootName;

    ///walker ensemble
    MCWalkerConfiguration& W;

    ///trial function
    TrialWaveFunction& Psi;

    ///Hamiltonian
    QMCHamiltonian& H;

    ///Observables manager
    ScalarEstimatorManager* Estimators;

    vector<TrialWaveFunction*> Psi1;

    vector<QMCHamiltonian*> H1;

    ///temporary storage for drift
    ParticleSet::ParticlePos_t drift;

    ///temporary storage for random displacement
    ParticleSet::ParticlePos_t deltaR;

    ///stream for the log file 
    OhmmsInform *LogOut;

    ///temporary buffer to accumulate data
    ostrstream log_buffer;

    PooledData<RealType> HamPool;
   

    ///Copy Constructor (disabled).
    QMCDriver(const QMCDriver& a): W(a.W), Psi(a.Psi), H(a.H), Estimators(0){}
 
    bool putQMCInfo(xmlNodePtr cur);

    void addWalkers(int nwalkers);  

    //void getReady();
  };
  
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
