//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CloneManager.h"
#include "MemoryUsage.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Concurrency/OpenMP.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/qmc_common.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
using TraceManager = int;
#endif

//comment this out to use only method to clone
#define ENABLE_CLONE_PSI_AND_H

namespace qmcplusplus
{
//initialization of the static wClones
UPtrVector<MCWalkerConfiguration> CloneManager::wClones_uptr;
std::vector<MCWalkerConfiguration*> CloneManager::wClones;
//initialization of the static psiClones
UPtrVector<TrialWaveFunction> CloneManager::psiClones_uptr;
std::vector<TrialWaveFunction*> CloneManager::psiClones;
//initialization of the static guideClones
UPtrVector<TrialWaveFunction> CloneManager::guideClones_uptr;
std::vector<TrialWaveFunction*> CloneManager::guideClones;

UPtrVector<MCWalkerConfiguration> CloneManager::wgClones;
//initialization of the static hClones
UPtrVector<QMCHamiltonian> CloneManager::hClones_uptr;
std::vector<QMCHamiltonian*> CloneManager::hClones;

std::vector<UPtrVector<MCWalkerConfiguration>> CloneManager::WPoolClones_uptr;
std::vector<std::vector<MCWalkerConfiguration*>> CloneManager::WPoolClones;
std::vector<UPtrVector<TrialWaveFunction>> CloneManager::PsiPoolClones_uptr;
std::vector<std::vector<TrialWaveFunction*>> CloneManager::PsiPoolClones;
std::vector<UPtrVector<QMCHamiltonian>> CloneManager::HPoolClones_uptr;
std::vector<std::vector<QMCHamiltonian*>> CloneManager::HPoolClones;

void CloneManager::clearClones()
{
  HPoolClones.clear();
  HPoolClones_uptr.clear();
  PsiPoolClones.clear();
  PsiPoolClones_uptr.clear();
  WPoolClones.clear();
  WPoolClones_uptr.clear();

  hClones.clear();
  hClones_uptr.clear();
  guideClones.clear();
  guideClones_uptr.clear();
  psiClones.clear();
  psiClones_uptr.clear();
  wClones.clear();
  wClones_uptr.clear();
}

/// Constructor.
CloneManager::CloneManager() : NumThreads(omp_get_max_threads()) { wPerRank.resize(NumThreads + 1, 0); }

///cleanup non-static data members
CloneManager::~CloneManager()
{
  // delete_iter(CSMovers.begin(),CSMovers.end());
  delete_iter(Movers.begin(), Movers.end());
  delete_iter(estimatorClones.begin(), estimatorClones.end());

#if !defined(REMOVE_TRACEMANAGER)
  delete_iter(traceClones.begin(), traceClones.end());
#endif
}

void CloneManager::makeClones(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& ham)
{
  if (wClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  wClones.resize(NumThreads);
  psiClones.resize(NumThreads);
  hClones.resize(NumThreads);
  wClones[0]   = &w;
  psiClones[0] = &psi;
  hClones[0]   = &ham;
  if (NumThreads == 1)
    return;

  wClones_uptr.resize(NumThreads - 1);
  psiClones_uptr.resize(NumThreads - 1);
  hClones_uptr.resize(NumThreads - 1);

  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  print_mem("Memory Usage before cloning", app_log());
  outputManager.pause();
  // clang-format off
  #pragma omp parallel
  {
    // check sizes
    #pragma omp master
    if (NumThreads != omp_get_num_threads())
      throw std::runtime_error("CloneManager::makeClones Inconsist NumThreads and omp_get_num_threads()!\n");

    const int ip = omp_get_thread_num();
    if (ip > 0)
    {
      // all the [ip] objects must be created on the ip threads to have first touch accurate.
      wClones_uptr[ip - 1]   = std::make_unique<MCWalkerConfiguration>(w);
      wClones[ip]            = wClones_uptr[ip - 1].get();
      psiClones_uptr[ip - 1] = psi.makeClone(*wClones[ip]);
      psiClones[ip]          = psiClones_uptr[ip-1].get();
      hClones_uptr[ip - 1]   = ham.makeClone(*wClones[ip], *psiClones[ip]);
      hClones[ip]            = hClones_uptr[ip-1].get();
    }
  }
  // clang-format on
  infoLog.resume();
  infoSummary.resume();
  print_mem("Memory Usage after cloning", app_log());
}


void CloneManager::makeClones(MCWalkerConfiguration& w,
                              std::vector<TrialWaveFunction*>& psipool,
                              std::vector<QMCHamiltonian*>& hampool)
{
  if (WPoolClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  IndexType nPsi = psipool.size();

  wClones.resize(NumThreads);
  PsiPoolClones.resize(NumThreads);
  HPoolClones.resize(NumThreads);
  wClones[0]       = &w;
  PsiPoolClones[0] = psipool;
  HPoolClones[0]   = hampool;

  if (NumThreads == 1)
    return;

  wClones_uptr.resize(NumThreads - 1);
  PsiPoolClones_uptr.resize(NumThreads - 1);
  HPoolClones_uptr.resize(NumThreads - 1);

  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H Pools." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  outputManager.pause();

  bool io_node       = qmc_common.io_node;
  qmc_common.io_node = false;

  for (int ip = 1; ip < NumThreads; ++ip)
  {
    PsiPoolClones[ip].resize(nPsi);
    PsiPoolClones_uptr[ip - 1].resize(nPsi);
    HPoolClones[ip].resize(nPsi);
    HPoolClones_uptr[ip - 1].resize(nPsi);

    wClones_uptr[ip - 1] = std::make_unique<MCWalkerConfiguration>(w);
    wClones[ip]          = wClones_uptr[ip - 1].get();
    for (int ipsi = 0; ipsi < psipool.size(); ipsi++)
    {
      PsiPoolClones_uptr[ip - 1][ipsi] = psipool[ipsi]->makeClone(w);
      PsiPoolClones[ip][ipsi]          = PsiPoolClones_uptr[ip - 1][ipsi].get();
      HPoolClones_uptr[ip - 1][ipsi]   = hampool[ipsi]->makeClone(w, *psipool[ipsi]);
      HPoolClones[ip][ipsi]            = HPoolClones_uptr[ip - 1][ipsi].get();
    }
  }
  infoSummary.resume();
  infoLog.resume();
  qmc_common.io_node = io_node;
}


void CloneManager::makeClones(TrialWaveFunction& guide)
{
  if (guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  guideClones.resize(NumThreads);
  guideClones[0] = &guide;
  if (NumThreads == 1)
    return;
  guideClones_uptr.resize(NumThreads - 1);
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." << std::endl;
  outputManager.pause();
  for (int ip = 1; ip < NumThreads; ++ip)
  {
    guideClones_uptr[ip - 1] = guide.makeClone(*wClones[ip]);
    guideClones[ip]          = guideClones_uptr[ip - 1].get();
  }
  infoSummary.resume();
  infoLog.resume();
}


void CloneManager::makeClones(MCWalkerConfiguration& wg, TrialWaveFunction& guide)
{
  if (guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  guideClones.resize(NumThreads);
  wgClones.resize(NumThreads);
  guideClones[0] = &guide;
  wgClones[0]    = std::make_unique<MCWalkerConfiguration>(wg);
  if (NumThreads == 1)
    return;
  guideClones_uptr.resize(NumThreads - 1);
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." << std::endl;
  outputManager.pause();
  for (int ip = 1; ip < NumThreads; ++ip)
  {
    wgClones[ip]             = std::make_unique<MCWalkerConfiguration>(wg);
    guideClones_uptr[ip - 1] = guide.makeClone(*wgClones[ip]);
    guideClones[ip]          = guideClones_uptr[ip - 1].get();
  }
  infoSummary.resume();
  infoLog.resume();
}

CloneManager::RealType CloneManager::acceptRatio() const
{
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  for (int ip = 0; ip < NumThreads; ip++)
  {
    nAcceptTot += Movers[ip]->nAccept;
    nRejectTot += Movers[ip]->nReject;
  }
#if defined(__GNUC__) || !defined(NDEBUG)
  // Attempt to detect compiler vectorization errors by computing
  // acceptance ratio in a different way to the above loop
  IndexType nAcceptTot_debug = 0;
  IndexType nRejectTot_debug = 0;
  std::vector<int> vec(NumThreads);
  for (int ip = 0; ip < NumThreads; ip++)
  {
    vec[ip] = Movers[ip]->nAccept;
    nAcceptTot_debug += vec[ip];
    vec[ip] = Movers[ip]->nReject;
    nRejectTot_debug += vec[ip];
  }
  if (nAcceptTot != nAcceptTot_debug || nRejectTot != nRejectTot_debug)
  {
    app_warning() << " Potential compiler bug detected!"
                  << " Overwriting nAcceptTot wrong value " << nAcceptTot << " with correct value " << nAcceptTot_debug
                  << "."
                  << " Overwriting nRejectTot wrong value " << nRejectTot << " with correct value " << nRejectTot_debug
                  << "." << std::endl;
    nAcceptTot = nAcceptTot_debug;
    nRejectTot = nRejectTot_debug;
  }
#endif
  return static_cast<RealType>(nAcceptTot) / static_cast<RealType>(nAcceptTot + nRejectTot);
}

} // namespace qmcplusplus
