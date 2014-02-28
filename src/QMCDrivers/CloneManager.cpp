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
#include "QMCDrivers/CloneManager.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/IteratorUtility.h"
#include <qmc_common.h>

//comment this out to use only method to clone
#define ENABLE_CLONE_PSI_AND_H

namespace qmcplusplus
{

//initialization of the static wClones
vector<MCWalkerConfiguration*> CloneManager::wClones;
//initialization of the static psiClones
vector<TrialWaveFunction*> CloneManager::psiClones;
//initialization of the static guideClones
vector<TrialWaveFunction*> CloneManager::guideClones;

vector<MCWalkerConfiguration*> CloneManager::wgClones;
//initialization of the static hClones
vector<QMCHamiltonian*> CloneManager::hClones;

vector<vector<MCWalkerConfiguration*> > CloneManager::WPoolClones; 
vector<vector<TrialWaveFunction*> > CloneManager::PsiPoolClones;
vector<vector<QMCHamiltonian*> > CloneManager::HPoolClones;

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
  delete_iter(traceClones.begin(),traceClones.end());
}

void CloneManager::makeClones(MCWalkerConfiguration& w,
                              TrialWaveFunction& psi, QMCHamiltonian& ham)
{
  if(wClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H." <<endl;
  app_log() << "  Cloning methods for both Psi and H are used" << endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();

  bool io_node=qmc_common.io_node;
  qmc_common.io_node=false;

  char pname[16];
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
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
  qmc_common.io_node=io_node;
}

void CloneManager::makeClones(vector<MCWalkerConfiguration*>& wpool,
                              vector<TrialWaveFunction*>& psipool, vector<QMCHamiltonian*>& hampool)
{
 /* if(WPoolClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H Pools." <<endl;
  app_log() << "  Cloning methods for both Psi and H are used" << endl;
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
                              vector<TrialWaveFunction*>& psipool, vector<QMCHamiltonian*>& hampool)
{
  if(WPoolClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H Pools." <<endl;
  app_log() << "  Cloning methods for both Psi and H are used" << endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();

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
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for W/Psi/H." <<endl;
  app_log() << "  Cloning methods for both Psi and H are used" << endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    wClones[ip]=new MCWalkerConfiguration(w);
    psiClones[ip]=psi.makeClone(*wClones[ip]);
    hClones[ip]=ham.makeClone(*wClones[ip],*psiClones[ip]);
  }
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
}

void CloneManager::makeClones(TrialWaveFunction& guide)
{
  if(guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." <<endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    guideClones[ip]=guide.makeClone(*wClones[ip]);
  }
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
}
void CloneManager::makeClones(MCWalkerConfiguration& wg, TrialWaveFunction& guide)
{
  if(guideClones.size())
  {
    app_log() << "  Cannot make clones again. Use existing " << NumThreads << " clones" << endl;
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
  app_log() << "  CloneManager::makeClones makes " << NumThreads << " clones for guide/wg." <<endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();
  char pname[16];
  for(int ip=1; ip<NumThreads; ++ip)
  {
    wgClones[ip]=new MCWalkerConfiguration(wg);
    guideClones[ip]=guide.makeClone(*wgClones[ip]);
  }
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
}

}

/***************************************************************************
 * $RCSfile: CloneManager.cpp,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/10/18 17:23:35 $
 * $Id:CloneManagerDMCPbyPOMP.cpp,v 1.5 2006/10/18 17:23:35 jnkim Exp $
 ***************************************************************************/
