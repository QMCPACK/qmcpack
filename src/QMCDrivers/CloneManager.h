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
/** @file CloneManager.h
 * @brief Manager class to handle multiple threads
 */
#ifndef QMCPLUSPLUS_CLONEMANAGER_H
#define QMCPLUSPLUS_CLONEMANAGER_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "CorrelatedSampling/CSUpdateBase.h"
// #include "QMCDrivers/EE/QMCRenyiUpdateBase.h"

namespace qmcplusplus
{

class HamiltonianPool;
//   class QMCRenyiUpdateBase;

/** Manager clones for threaded applications
 *
 * Clones for the ParticleSet, TrialWaveFunction and QMCHamiltonian
 * are static to ensure only one set of clones persist during a run.
 */
class CloneManager: public QMCTraits
{
public:
  /// Constructor.
  CloneManager(HamiltonianPool& hpool);
  ///virtual destructor
  virtual ~CloneManager();

  void makeClones(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& ham);
  void makeClones(MCWalkerConfiguration& w, std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& ham);
  void makeClones(std::vector<MCWalkerConfiguration*>& w, std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& ham);
  void makeClones_new(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& ham);
  void makeClones(MCWalkerConfiguration& wg, TrialWaveFunction& guide);
  void makeClones(TrialWaveFunction& guide);

  inline RealType acceptRatio() const
  {
    IndexType nAcceptTot=0;
    IndexType nRejectTot=0;
    for(int ip=0; ip<NumThreads; ip++)
    {
      nAcceptTot+=Movers[ip]->nAccept;
      nRejectTot+=Movers[ip]->nReject;
    }
    return static_cast<RealType>(nAcceptTot)/static_cast<RealType>(nAcceptTot+nRejectTot);
  }

protected:
  ///reference to HamiltonianPool to clone everything
  HamiltonianPool& cloneEngine;
  ///number of threads
  IndexType NumThreads;
  ///walkers
  static std::vector<MCWalkerConfiguration*> wClones;
  static std::vector<MCWalkerConfiguration*> wgClones;
  ///trial wavefunctions
  static std::vector<TrialWaveFunction*> psiClones;
  ///guide wavefunctions
  static std::vector<TrialWaveFunction*> guideClones;
  ///Hamiltonians
  static std::vector<QMCHamiltonian*> hClones;
  ///update engines
  std::vector<QMCUpdateBase*> Movers;
//     ///update engines
//     std::vector<QMCRenyiUpdateBase*> RenyiMovers;
  ///estimator managers
  std::vector<EstimatorManager*> estimatorClones;
  ///trace managers
  std::vector<TraceManager*> traceClones;
  ///Branch engines
  std::vector<SimpleFixedNodeBranch*> branchClones;
  
  //for correlated sampling.
  static std::vector<std::vector<MCWalkerConfiguration*> > WPoolClones; 
  static std::vector<std::vector<TrialWaveFunction*> > PsiPoolClones;
  static std::vector<std::vector<QMCHamiltonian*> > HPoolClones;
  std::vector<CSUpdateBase*> CSMovers;

  
  
  ///Walkers per node
  std::vector<int> wPerNode;
};
}
#endif
/***************************************************************************
 * $RCSfile: CloneManager.h,v $   $Author: jnkim $
 * $Revision: 1.2 $   $Date: 2006/02/26 17:41:10 $
 * $Id: CloneManager.h,v 1.2 2006/02/26 17:41:10 jnkim Exp $
 ***************************************************************************/
