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
class CloneManager : public QMCTraits
{
public:
  /// Constructor.
  CloneManager();
  ///virtual destructor
  virtual ~CloneManager();

  // Clear static arrays of clones owned by the manager
  static void clearClones();

  void makeClones(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& ham);
  void makeClones(MCWalkerConfiguration& w, std::vector<TrialWaveFunction*>& psi, std::vector<QMCHamiltonian*>& ham);
  void makeClones(MCWalkerConfiguration& wg, TrialWaveFunction& guide);
  void makeClones(TrialWaveFunction& guide);

  RealType acceptRatio() const;

protected:
  ///number of threads
  const IndexType NumThreads;
  /* A Non-owning pointer is passed into CloneManger, used by thread 0 and stored in XXX[0].
   * Clones for all other threads are owned by CloneManager. For thread N, XXX_uptr[N-1] stores 
   * a unique_ptr and a non-owning copy of the pointer is stored in XXX[N].
   */
  ///walkers
  static UPtrVector<MCWalkerConfiguration> wClones_uptr;
  static std::vector<MCWalkerConfiguration*> wClones;
  static UPtrVector<MCWalkerConfiguration> wgClones;
  ///trial wavefunctions
  static UPtrVector<TrialWaveFunction> psiClones_uptr;
  static std::vector<TrialWaveFunction*> psiClones;
  ///guide wavefunctions
  static UPtrVector<TrialWaveFunction> guideClones_uptr;
  static std::vector<TrialWaveFunction*> guideClones;
  ///Hamiltonians
  static UPtrVector<QMCHamiltonian> hClones_uptr;
  static std::vector<QMCHamiltonian*> hClones;
  ///update engines
  std::vector<QMCUpdateBase*> Movers;
  ///estimator managers
  std::vector<EstimatorManagerBase*> estimatorClones;
  ///trace managers
  std::vector<TraceManager*> traceClones;

  //for correlated sampling.
  static std::vector<UPtrVector<MCWalkerConfiguration>> WPoolClones_uptr;
  static std::vector<std::vector<MCWalkerConfiguration*>> WPoolClones;
  static std::vector<UPtrVector<TrialWaveFunction>> PsiPoolClones_uptr;
  static std::vector<std::vector<TrialWaveFunction*>> PsiPoolClones;
  static std::vector<UPtrVector<QMCHamiltonian>> HPoolClones_uptr;
  static std::vector<std::vector<QMCHamiltonian*>> HPoolClones;
  UPtrVector<CSUpdateBase> CSMovers;

  ///Walkers per MPI rank
  std::vector<int> wPerRank;
};
} // namespace qmcplusplus
#endif
