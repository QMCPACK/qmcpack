//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <map>
namespace qmcplusplus {

MCWalkerConfiguration::MCWalkerConfiguration(): 
OwnWalkers(true),ReadyForPbyP(false),UpdateMode(Update_Walker),Polymer(0) {
  //move to ParticleSet
  //initPropertyList();
}

MCWalkerConfiguration::MCWalkerConfiguration(const MCWalkerConfiguration& mcw)
: ParticleSet(mcw), OwnWalkers(true),
  UpdateMode(Update_Walker), ReadyForPbyP(false), Polymer(0)
{
  //initPropertyList();
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
MCWalkerConfiguration::destroyWalkers(int nw) {
  if(WalkerList.size() == 1 || nw >= WalkerList.size()) {
    app_warning() << "  Cannot remove walkers. Current Walkers = " << WalkerList.size() << endl;
    return;
  }
  nw=WalkerList.size()-nw;
  iterator it(WalkerList.begin()+nw),it_end(WalkerList.end());
  while(it != it_end) {
    delete *it++;
  }
  WalkerList.erase(WalkerList.begin()+nw,WalkerList.end());
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
  iterator it(WalkerList.begin()), it_end(WalkerList.end());
  while(it != it_end) {//(*it)->reset();++it;}
    (*it)->Weight=1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
}

void MCWalkerConfiguration::clearAuxDataSet() {
  UpdateMode=Update_Particle;
  iterator it(WalkerList.begin());
  iterator it_end(WalkerList.end());
  while(it!=it_end) {
    (*it)->DataSet.clear(); ++it;
  }
  ReadyForPbyP = false;
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

//void 
//MCWalkerConfiguration::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
//  R = awalker.R;
//  for(int i=0; i< DistTables.size(); i++) {
//    DistTables[i]->evaluate(*this);
//    DistTables[i]->registerData(buf);
//  }
//}
//
//void 
//MCWalkerConfiguration::updateBuffer(Walker_t& awalker, PooledData<RealType>& buf) {
//  R = awalker.R;
//  for(int i=0; i< DistTables.size(); i++) {
//    DistTables[i]->evaluate(*this);
//    DistTables[i]->updateBuffer(buf);
//  }
//}
//  
//void 
//MCWalkerConfiguration::copyToBuffer(PooledData<RealType>& buf) {
//  for(int i=0; i< DistTables.size(); i++) {
//    DistTables[i]->copyToBuffer(buf);
//  }
//}
//
//void 
//MCWalkerConfiguration::copyFromBuffer(PooledData<RealType>& buf) {
//  for(int i=0; i< DistTables.size(); i++) {
//    DistTables[i]->copyFromBuffer(buf);
//  }
//}


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
  app_log() << "  Resetting Properties of the walkers " << ncopy << " x " << m << endl;
  iterator it(WalkerList.begin()),it_end(WalkerList.end());
  while(it != it_end) {
    (*it)->resizeProperty(ncopy,m); ++it;
  }
}

//void MCWalkerConfiguration::initPropertyList() {
//  //Need to add the default Properties according to the enumeration
//  PropertyList.add("LogPsi");
//  PropertyList.add("SignPsi");
//  PropertyList.add("UmbrellaWeight");
//  PropertyList.add("LocalEnergy");
//  PropertyList.add("LocalPotential");
//}

//int 
//MCWalkerConfiguration::branch(int maxcopy, int Nmax, int Nmin, bool swap) {
//
//  iterator it = WalkerList.begin();
//  int iw=0, nw = WalkerList.size();
//
//  vector<Walker_t*> good, bad;
//  vector<int> ncopy;
//  ncopy.reserve(nw);
//
//  int num_walkers=0;
//  while(it != WalkerList.end()) {
//    int nc = std::min(static_cast<int>((*it)->Multiplicity),maxcopy);
//    if(nc) {
//      num_walkers += nc;
//      good.push_back(*it);
//      ncopy.push_back(nc-1);
//    } else {
//      bad.push_back(*it);
//    }
//    iw++;it++;
//  }
//
//  //remove bad walkers
//  for(int i=0; i<bad.size(); i++) delete bad[i];
//
//  if(good.empty()) {
//    ERRORMSG("All the walkers have died. Abort. ")
//    OHMMS::Controller->abort();
//  }
//
//  //check if the projected number of walkers is too small or too large
//  if(num_walkers>Nmax) {
//    int nsub=0;
//    int nsub_target=num_walkers-static_cast<int>(0.9*Nmax);
//    int i=0;
//    while(i<ncopy.size() && nsub<nsub_target) {
//      if(ncopy[i]) {ncopy[i]--; nsub++;}
//      i++;
//    }
//    num_walkers -= nsub;
//  } else  if(num_walkers < Nmin) {
//    int nadd=0;
//    int nadd_target = static_cast<int>(Nmin*1.1)-num_walkers;
//    if(nadd_target>good.size()) {
//      WARNMSG("The number of walkers is running low. Requested walkers " << nadd_target << " good walkers = " << good.size())
//    }
//    int i=0;
//    while(i<ncopy.size() && nadd<nadd_target) {
//      ncopy[i]++; nadd++;i++;
//    }
//    num_walkers +=  nadd;
//  }
//
//  LOGMSG("Projected number of walkers " << num_walkers)
//
//  //WalkerControl
//  //MPI Send to the master, MPI Irecv by the master
//  //send the total number of walkers to the master
// 
//  //clear the WalkerList to populate them with the good walkers
//  WalkerList.clear();
//  WalkerList.insert(WalkerList.begin(), good.begin(), good.end());
//
//  int cur_walker = good.size();
//  for(int i=0; i<good.size(); i++) { //,ie+=ncols) {
//    for(int j=0; j<ncopy[i]; j++, cur_walker++) {
//      WalkerList.push_back(new Walker_t(*(good[i])));
//    }
//  }
//
//  int nw_tot = WalkerList.size();
//  LOGMSG("Real number of walkers " << nw_tot)
//
//  //WalkerControl
//  //Master check if criteria is met and send back 0/1, total walkers, max, min
//
//  if(swap) nw_tot= swapWalkers();
//  //if(swap) gsum(nw_tot,0);
//
//  //set Weight and Multiplicity to default values
//  iw=0;
//  it=WalkerList.begin();
//  while(it != WalkerList.end()) {
//    (*it)->Weight= 1.0;
//    (*it)->Multiplicity=1.0;
//    it++;
//  }
//
//  return nw_tot;
//}
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
