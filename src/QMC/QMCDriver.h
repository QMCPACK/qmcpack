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
#include "Estimators/ScalarEstimatorManager.h"
#include <strstream>
class OhmmsInform;

namespace ohmmsqmc {

  class MCWalkerConfiguration;
  class TrialWaveFunction;
  class QMCHamiltonian;

  /** Base class to perform QMC simulations. */
  class QMCDriver: public QMCTraits {

  public:
    /// Constructor.
    QMCDriver(MCWalkerConfiguration& w, 
	      TrialWaveFunction& psi, 
	      QMCHamiltonian& h, 
	      xmlNodePtr q);

    virtual ~QMCDriver();

    void setFileRoot(const string& aname);

    virtual bool run() = 0;
    
    virtual bool put(xmlNodePtr cur) = 0;

  protected:

    ///counts the number of qmc runs
    static int Counter;

    int AcceptIndex;

    ///timestep
    RealType Tau;

    RealType FirstStep;

    ///reference energy, only used in DMC
    RealType e_ref;

    ///maximum number of blocks
    IndexType nBlocks;

    ///maximum number of steps
    IndexType nSteps;

    ///flag to print walker ensemble
    bool pStride;

    ///counter for number of moves accepted
    int nAccept;

    ///counter for number of moves /rejected
    int nReject; 

    ///the number of walkers
    int nTargetWalkers;

    ParameterSet m_param;

    ///type of qmc: assigned by subclasses
    string QMCType;

    ///root of all the output files
    string RootName;

    ///Observables manager
    ScalarEstimatorManager Estimators;

    ///walker ensemble
    MCWalkerConfiguration& W;

    ///trial function
    TrialWaveFunction& Psi;

    ///Hamiltonian
    QMCHamiltonian& H;

    ///stream for the log file 
    OhmmsInform *LogOut;

    ///temporary buffer to accumulate data
    ostrstream log_buffer;

    PooledData<RealType> HamPool;
   
    ///pointer to qmc node in xml file
    xmlNodePtr qmc_node;

    ///Copy Constructor (disabled).
    QMCDriver(const QMCDriver& a): 
      W(a.W), Psi(a.Psi), H(a.H), qmc_node(a.qmc_node),Estimators(H),
      nAccept(a.nAccept), nReject(a.nReject), Tau(a.Tau), 
      e_ref(a.e_ref), nBlocks(a.nBlocks), nSteps(a.nSteps) { }
 
    bool putQMCInfo(xmlNodePtr cur);

    void addWalkers(int nwalkers);  

    void getReady();
  };
  
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
