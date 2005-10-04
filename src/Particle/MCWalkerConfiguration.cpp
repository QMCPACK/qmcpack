//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include <map>
using namespace ohmmsqmc;

MCWalkerConfiguration::MCWalkerConfiguration(): 
OwnWalkers(true),ReadyForPbyP(false),UpdateMode(Update_Walker) {
  initPropertyList();
}

MCWalkerConfiguration::MCWalkerConfiguration(const MCWalkerConfiguration& mcw)
: ParticleSet(mcw), OwnWalkers(true),
  UpdateMode(Update_Walker), ReadyForPbyP(false)
{
  initPropertyList();
}

///default destructor
MCWalkerConfiguration::~MCWalkerConfiguration(){
  DEBUGMSG("MCWalkerConfiguration::~MCWalkerConfiguration")
  if(OwnWalkers) destroyWalkers(WalkerList.begin(), WalkerList.end());
}


void MCWalkerConfiguration::createWalkers(int n) {
  while(n) {
    Walker_t* awalker=new Walker_t(GlobalNum);
    awalker->R = R;
    awalker->Drift = 0.0;
    WalkerList.push_back(awalker);
    --n;
  }
}

void MCWalkerConfiguration::resize(int numWalkers, int numPtcls) {

  WARNMSG("MCWalkerConfiguration::resize cleans up the walker list.")

  ParticleSet::resize(unsigned(numPtcls));

  int dn=numWalkers-WalkerList.size();
  if(dn>0) createWalkers(dn);

  if(dn<0) {
    int nw=-dn;
    if(nw<WalkerList.size())  {
      iterator it = WalkerList.begin();
      while(nw) {
        delete *it; ++it; --nw;
      }
      WalkerList.erase(WalkerList.begin(),WalkerList.begin()-dn);
    }
  }
  //iterator it = WalkerList.begin();
  //while(it != WalkerList.end()) {
  //  delete *it++;
  //}
  //WalkerList.erase(WalkerList.begin(),WalkerList.end());
  //R.resize(np);
  //GlobalNum = np;
  //createWalkers(nw);  
}

///returns the next valid iterator
MCWalkerConfiguration::iterator 
MCWalkerConfiguration::destroyWalkers(iterator first, iterator last) {
  if(OwnWalkers) {
    iterator it = first;
    while(it != last) { delete *it++;}
  }
  return WalkerList.erase(first,last);
}

void 
MCWalkerConfiguration::copyWalkerRefs(Walker_t* head, Walker_t* tail) {

  if(OwnWalkers) { //destroy the current walkers
    destroyWalkers(WalkerList.begin(), WalkerList.end());
    WalkerList.clear();
    OwnWalkers=false;//set to false to prevent deleting the Walkers
  }

  if(WalkerList.size()<2) {
    WalkerList.push_back(0);
    WalkerList.push_back(0);
  }

  WalkerList[0]=head;
  WalkerList[1]=tail;
}

/** Make Metropolis move to the walkers and save in a temporary array.
 * @param it the iterator of the first walker to work on
 * @param tauinv  inverse of the time step
 *
 * R + D + X
 */
void MCWalkerConfiguration::sample(iterator it, RealType tauinv) {
  makeGaussRandom(R);
  R *= tauinv;
  R += (*it)->R + (*it)->Drift;
}

void MCWalkerConfiguration::reset() {
  iterator it=WalkerList.begin();
  while(it != WalkerList.end()) {(*it)->reset();it++;}
}

bool MCWalkerConfiguration::createAuxDataSet(int nfield) {

  if(ReadyForPbyP) return false;

  ReadyForPbyP=true;
  UpdateMode=Update_Particle;
  iterator it(WalkerList.begin());
  iterator it_end(WalkerList.end());
  while(it!=it_end) {
    (*it)->DataSet.reserve(nfield); ++it;
  }

  return true;
}

void 
MCWalkerConfiguration::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
  R = awalker.R;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->evaluate(*this);
    DistTables[i]->registerData(buf);
  }
}

void 
MCWalkerConfiguration::updateBuffer(Walker_t& awalker, PooledData<RealType>& buf) {
  R = awalker.R;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->evaluate(*this);
    DistTables[i]->updateBuffer(buf);
  }
}
  
void 
MCWalkerConfiguration::copyToBuffer(PooledData<RealType>& buf) {
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->copyToBuffer(buf);
  }
}

void 
MCWalkerConfiguration::copyFromBuffer(PooledData<RealType>& buf) {
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->copyFromBuffer(buf);
  }
}

int MCWalkerConfiguration::branch(int maxcopy, int Nmax, int Nmin) {

  iterator it = WalkerList.begin();
  int iw=0, nw = WalkerList.size();

  vector<Walker_t*> good, bad;
  vector<int> ncopy;
  ncopy.reserve(nw);

  int num_walkers=0;
  while(it != WalkerList.end()) {
    int nc = std::min(static_cast<int>((*it)->Multiplicity),maxcopy);
    if(nc) {
      num_walkers += nc;
      good.push_back(*it);
      ncopy.push_back(nc-1);
    } else {
      bad.push_back(*it);
    }
    iw++;it++;
  }

  //remove bad walkers
  for(int i=0; i<bad.size(); i++) delete bad[i];

  if(good.empty()) {
    ERRORMSG("All the walkers have died. Abort. ")
    OHMMS::Controller->abort();
  }

  //check if the projected number of walkers is too small or too large
  if(num_walkers>Nmax) {
    int nsub=0;
    int nsub_target=num_walkers-static_cast<int>(0.9*Nmax);
    int i=0;
    while(i<ncopy.size() && nsub<nsub_target) {
      if(ncopy[i]) {ncopy[i]--; nsub++;}
      i++;
    }
    num_walkers -= nsub;
  } else  if(num_walkers < Nmin) {
    int nadd=0;
    int nadd_target = static_cast<int>(Nmin*1.1)-num_walkers;
    if(nadd_target>good.size()) {
      WARNMSG("The number of walkers is running low. Requested walkers " << nadd_target << " good walkers = " << good.size())
    }
    int i=0;
    while(i<ncopy.size() && nadd<nadd_target) {
      ncopy[i]++; nadd++;i++;
    }
    num_walkers +=  nadd;
  }

  //clear the WalkerList to populate them with the good walkers
  WalkerList.clear();
  WalkerList.insert(WalkerList.begin(), good.begin(), good.end());

  int cur_walker = good.size();
  for(int i=0; i<good.size(); i++) { //,ie+=ncols) {
    for(int j=0; j<ncopy[i]; j++, cur_walker++) {
      WalkerList.push_back(new Walker_t(*(good[i])));
    }
  }

  //DMC+MPI: swapWalkers() is disabled
  //int nw_tot = swapWalkers();
  int nw_tot = WalkerList.size();

  //set Weight and Multiplicity to default values
  iw=0;
  it=WalkerList.begin();
  while(it != WalkerList.end()) {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    it++;
  }

  return nw_tot;
}

void MCWalkerConfiguration::loadWalker(Walker_t& awalker) {
  R = awalker.R;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->evaluate(*this);
  }
}

/** reset the Property container of all the walkers
 */
void MCWalkerConfiguration::resetWalkerProperty(int ncopy) {
  int m(PropertyList.size());
  LOGMSG("Resetting Properties of the walkers " << ncopy << " x " << m)
  iterator it(WalkerList.begin()),it_end(WalkerList.end());
  while(it != it_end) {
    (*it)->resizeProperty(ncopy,m); ++it;
  }
}

void MCWalkerConfiguration::initPropertyList() {
  //Need to add the default Properties according to the enumeration
  PropertyList.add("LogPsi");
  PropertyList.add("SignPsi");
  PropertyList.add("UmbrellaWeight");
  PropertyList.add("LocalEnergy");
  PropertyList.add("LocalPotential");
}

#if defined(HAVE_MPI)

int MCWalkerConfiguration::swapWalkers() {

  //synchronize the nodes
  OHMMS::Controller->barrier();

  int mycontext=OHMMS::Controller->mycontext();
  int nrecv, ncontexts=OHMMS::Controller->ncontexts();
  vector<int> nsub(ncontexts,0), nsub_g(ncontexts,0);
  nsub[mycontext]=WalkerList.size();

  //get how many each node has
  int status=MPI_Allreduce(&nsub[0], &nsub_g[0], ncontexts, 
      MPI_INT, MPI_SUM,OHMMS::Controller->getID());

  //order it according to the number of particles
  multimap<int,int> nw_map;
  int nw_sum=0;
  for(int i=0; i<ncontexts; i++) {
    nw_sum+=nsub_g[i];
    nw_map.insert(pair<int,int>(nsub_g[i],i));
  }

  multimap<int,int>::iterator it(nw_map.begin());
  multimap<int,int>::reverse_iterator it_b(nw_map.end());
  bool notpaired=true;
  int target_context=-1;
  int half=ncontexts/2;
  int item=0;
  bool minorcontext;
  while(notpaired &&item<half) {
    int i=(*it).second;
    int j=(*it_b).second;
    if(i == mycontext) {
      target_context=j;
      notpaired=false;
      minorcontext=true;
    } else if(j == mycontext) {
      target_context= i;
      notpaired=false;
      minorcontext=false;
    } 
    ++it; ++it_b; ++item;
  }

  int nw_tot=nsub_g[mycontext]+nsub_g[target_context];
  int nw_L=nw_tot/2;
  int nw_R=nw_tot-nw_L;
  int dnw(0);
  if(minorcontext) {
    dnw=nw_R-nsub_g[mycontext];
  } else {
    dnw=nw_R-nsub_g[target_context];
  }

  //char fname[128];
  //sprintf(fname,"test.%d",mycontext);
  //ofstream fout(fname,ios::app);

  //if(minorcontext) 
  //  fout << mycontext << " recv from " << target_context <<  " " << dnw << endl;
  //else 
  //  fout << mycontext << " send to " << target_context << " " << dnw << endl;
  if(dnw) {//something to swap
    if(minorcontext) {//open recv buffer
      Walker_t& wRef(*WalkerList[0]);
      OOMPI_Packed recvBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      //To check if irecv is better than recv
      //OOMPI_COMM_WORLD[target_context].Recv(recvBuffer);
      OOMPI_Request recvRequest = OOMPI_COMM_WORLD[target_context].Irecv(recvBuffer, MPI_ANY_TAG);

      //create walkers
      for(int iw=0; iw<dnw; iw++) {
        WalkerList.push_back(new Walker_t(wRef));
      }
      recvRequest.Wait();
      int last=nsub_g[mycontext];
      while(dnw) {
        WalkerList[last++]->getMessage(recvBuffer);
        --dnw;
      }

    } else {
      Walker_t& wRef(*WalkerList[0]);
      OOMPI_Packed sendBuffer(dnw*wRef.byteSize(),OOMPI_COMM_WORLD);
      int last=WalkerList.size()-1;
      while(dnw) {
        WalkerList[last--]->putMessage(sendBuffer);
        --dnw; 
      }
      OOMPI_COMM_WORLD[target_context].Send(sendBuffer);
      destroyWalkers(WalkerList.begin()+nsub_g[mycontext], WalkerList.end());
    }
  }

  OHMMS::Controller->barrier();
  return nw_sum;
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
