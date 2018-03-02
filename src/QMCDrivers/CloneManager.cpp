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
    
    

#include "Platforms/sysutil.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/IteratorUtility.h"
#include <qmc_common.h>
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

//comment this out to use only method to clone
#define ENABLE_CLONE_PSI_AND_H

namespace qmcplusplus
{

//initialization of the static wClones
std::vector<MCWalkerConfiguration*> CloneManager::wClones;
//initialization of the static psiClones
std::vector<TrialWaveFunction*> CloneManager::psiClones;
//initialization of the static guideClones
std::vector<TrialWaveFunction*> CloneManager::guideClones;

std::vector<MCWalkerConfiguration*> CloneManager::wgClones;
//initialization of the static hClones
std::vector<QMCHamiltonian*> CloneManager::hClones;

std::vector<std::vector<MCWalkerConfiguration*> > CloneManager::WPoolClones; 
std::vector<std::vector<TrialWaveFunction*> > CloneManager::PsiPoolClones;
std::vector<std::vector<QMCHamiltonian*> > CloneManager::HPoolClones;

/// Constructor.
CloneManager::CloneManager(HamiltonianPool& hpool): cloneEngine(hpool)
{
  NumThreads=omp_get_max_threads();
  wPerNode.resize(NumThreads+1,0);
}

///cleanup non-static data members
CloneManager::~CloneManager()
{
  // delete_iter(CSMovers.begin(),CSMovers.end());
  delete_iter(Movers.begin(),Movers.end());
  delete_iter(branchClones.begin(),branchClones.end());
  delete_iter(estimatorClones.begin(),estimatorClones.end());

#if !defined(REMOVE_TRACEMANAGER)
  delete_iter(traceClones.begin(),traceClones.end());
#endif
}

void CloneManager::makeClones(MCWalkerConfiguration& w,
                              TrialWaveFunction& psi, QMCHamiltonian& ham)
{
  if(wClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  wClones.resize(NumThreads);
  psiClones.resize(NumThreads);
  hClones.resize(NumThreads);
  wClones[0]=&w;
  psiClones[0]=&psi;
  hClones[0]=&ham;
  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  print_mem("Memory Usage before cloning", app_log());
  outputManager.pause();

  #pragma omp parallel for shared(w,psi,ham)
  for(int ip=1; ip<NumThreads; ++ip)
  {
#if defined(USE_PARTCILESET_CLONE)
    wClones[ip]=dynamic_cast<MCWalkerConfiguration*>(w.get_clone(ip));
#else
    wClones[ip]=new MCWalkerConfiguration(w);
#endif
    psiClones[ip]=psi.makeClone(*wClones[ip]);
    hClones[ip]=ham.makeClone(*wClones[ip],*psiClones[ip]);
  }
  infoLog.resume();
  infoSummary.resume();
  print_mem("Memory Usage after cloning", app_log());
}

void CloneManager::makeClones(std::vector<MCWalkerConfiguration*>& wpool,
                              std::vector<TrialWaveFunction*>& psipool, std::vector<QMCHamiltonian*>& hampool)
{
 /* if(WPoolClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  WPoolClones.resize(NumThreads);
  PsiPoolClones.resize(NumThreads);
  HPoolClones.resize(NumThreads);
  WPoolClones[0]=wpool;
  PsiPoolClones[0]=psipool;
  HPoolClones[0]=hampool;
  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H Pools." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();

  bool io_node=qmc_common.io_node;
  qmc_common.io_node=false;

  char pname[16];
  
  for(int ip=1; ip<NumThreads; ++ip)
  {
	IndexType nPsi=psipool.size();
    WPoolClones[ip].resize(nPsi,0);
    PsiPoolClones[ip].resize(nPsi,0);
    HPoolClones[ip].resize(nPsi,0);
    for(int ipsi=0; ipsi<psipool.size(); ipsi++)
    {

//#if defined(USE_PARTCILESET_CLONE)
//      wClones[ip]=dynamic_cast<MCWalkerConfiguration*>(w.get_clone(ip));
//#else
      WPoolClones[ip][ipsi]=new MCWalkerConfiguration(*wpool[ipsi]);
//#endif
      PsiPoolClones[ip][ipsi]=psipool[ipsi]->makeClone(*WPoolClones[ip][ipsi]);
      HPoolClones[ip][ipsi]=hampool[ipsi]->makeClone(*WPoolClones[ip][ipsi],*psiClones[ip]);
    }
  }
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
  qmc_common.io_node=io_node;*/
}

void CloneManager::makeClones(MCWalkerConfiguration& w,
                              std::vector<TrialWaveFunction*>& psipool, std::vector<QMCHamiltonian*>& hampool)
{
  if(WPoolClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  IndexType nPsi=psipool.size();
  
  wClones.resize(NumThreads);
  PsiPoolClones.resize(NumThreads);
  HPoolClones.resize(NumThreads);
  wClones[0]=&w;
  PsiPoolClones[0]=psipool;
  HPoolClones[0]=hampool;
  

  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H Pools." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  outputManager.pause();

  bool io_node=qmc_common.io_node;
  qmc_common.io_node=false;

  char pname[16];
 
  for(int ip=1; ip<NumThreads; ++ip)
  {
	
   // WPoolClones[ip].resize(nPsi,0);
    PsiPoolClones[ip].resize(nPsi);
    HPoolClones[ip].resize(nPsi);
    for(int ipsi=0; ipsi<psipool.size(); ipsi++)
    {

//#if defined(USE_PARTCILESET_CLONE)
//      wClones[ip]=dynamic_cast<MCWalkerConfiguration*>(w.get_clone(ip));
//#else
      wClones[ip]=new MCWalkerConfiguration(w);
//#endif
      PsiPoolClones[ip][ipsi]=psipool[ipsi]->makeClone(w);
      HPoolClones[ip][ipsi]=hampool[ipsi]->makeClone(w,*psipool[ipsi]);
    }
  }
  infoSummary.resume();
  infoLog.resume();
  qmc_common.io_node=io_node;
}


void CloneManager::makeClones_new(MCWalkerConfiguration& w,
                                  TrialWaveFunction& psi, QMCHamiltonian& ham)
{
  if(wClones.size())
  {
    delete_iter(wClones.begin()+1,wClones.end());
    delete_iter(psiClones.begin()+1,psiClones.end());
    delete_iter(hClones.begin()+1,hClones.end());
  }
  else
  {
    wClones.resize(NumThreads,0);
    psiClones.resize(NumThreads,0);
    hClones.resize(NumThreads,0);
  }
  wClones[0]=&w;
  psiClones[0]=&psi;
  hClones[0]=&ham;
  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H." << std::endl;
  app_log() << "  Cloning methods for both Psi and H are used" << std::endl;
  outputManager.pause();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    wClones[ip]=new MCWalkerConfiguration(w);
    psiClones[ip]=psi.makeClone(*wClones[ip]);
    hClones[ip]=ham.makeClone(*wClones[ip],*psiClones[ip]);
  }
  infoSummary.resume();
  infoLog.resume();
}

void CloneManager::makeClones(TrialWaveFunction& guide)
{
  if(guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  else
  {
    guideClones.resize(NumThreads,0);
  }
  guideClones.resize(NumThreads,0);
  guideClones[0]=&guide;
  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." << std::endl;
  outputManager.pause();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    guideClones[ip]=guide.makeClone(*wClones[ip]);
  }
  infoSummary.resume();
  infoLog.resume();
}
void CloneManager::makeClones(MCWalkerConfiguration& wg, TrialWaveFunction& guide)
{
  if(guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << std::endl;
    return;
  }
  else
  {
    guideClones.resize(NumThreads,0);
    wgClones.resize(NumThreads,0);
  }
  guideClones.resize(NumThreads,0);
  wgClones.resize(NumThreads,0);
  guideClones[0]=&guide;
  wgClones[0]=new MCWalkerConfiguration(wg);
  if(NumThreads==1)
    return;
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." << std::endl;
  outputManager.pause();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    wgClones[ip]=new MCWalkerConfiguration(wg);
    guideClones[ip]=guide.makeClone(*wgClones[ip]);
  }
  infoSummary.resume();
  infoLog.resume();
}

}

