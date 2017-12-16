//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file CloneManager.h
 * @brief Manager class to handle multiple threads
 */
#ifndef QMCPLUSPLUS_CLONEMANAGER_H
#define QMCPLUSPLUS_CLONEMANAGER_H
#include "QMCDrivers/QMCUpdateBase.h"
#include "CorrelatedSampling/CSUpdateBase.h"

namespace qmcplusplus
{

class HamiltonianPool;

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
  ///estimator managers
  std::vector<EstimatorManagerBase*> estimatorClones;
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
