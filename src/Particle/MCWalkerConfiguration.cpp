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
using namespace ohmmsqmc;
#include "ParticleBase/RandomSeqGenerator.h"

MCWalkerConfiguration::MCWalkerConfiguration(): 
UpdateMode(Update_Walker) {

}


///default destructor
MCWalkerConfiguration::~MCWalkerConfiguration(){

  destroyWalker(WalkerList.begin(), WalkerList.end());
  for(int iw=0; iw<DataSet.size(); iw++) delete DataSet[iw];

  DEBUGMSG("MCWalkerConfiguration::~MCWalkerConfiguration")
 //WalkerList.clear();
}

void MCWalkerConfiguration::createWalkers(int n) {
  for(int i=0; i<n; i++) {
    WalkerList.push_back(new Walker_t(GlobalNum));
    WalkerList.back()->R = R;
    WalkerList.back()->Drift = 0.0;
  }
}

void MCWalkerConfiguration::resize(int nw, int np) {
  WARNMSG("MCWalkerConfiguration::resize cleans up the walker list.")
  iterator it = WalkerList.begin();
  while(it != WalkerList.end()) {
    delete *it++;
  }
  WalkerList.erase(WalkerList.begin(),WalkerList.end());
  R.resize(np);
  GlobalNum = np;
  createWalkers(nw);  
}


void MCWalkerConfiguration::addWalkers(vector<int>& Copies,
				       vector<Walker_t*>& InactiveList) {
  for(int i=0; i<Copies.size(); i++) {
    for(int j=0; j<Copies[i]; j++) {
      WalkerList.push_back(new Walker_t(*InactiveList[i]));
    }
  }
}


void MCWalkerConfiguration::copyWalker(iterator walkerit, int n) {

  for(int i=0; i<n; i++)
    WalkerList.push_back(new Walker_t(**walkerit));
    //WalkerList.insert(WalkerList.end(),new Walker_t(**walkerit));
}

///returns the next valid iterator
MCWalkerConfiguration::iterator 
MCWalkerConfiguration::destroyWalker(iterator walkerit) {
  delete *walkerit;
  return WalkerList.erase(walkerit);
}

///returns the next valid iterator
MCWalkerConfiguration::iterator 
MCWalkerConfiguration::destroyWalker(iterator itstart, iterator itend) {
  iterator it = itstart;
  while(it != itend) { delete *it++;}
  return WalkerList.erase(itstart,itend);
}

/**@param first the iterator of the first walker to work on
 *@param last the iterator of the last walker to work on
 *@brief Make Metropolis move to the walkers and save in a temporary array.
 */
void MCWalkerConfiguration::sample(iterator it, RealType tauinv) {
  /// R + D + X
  makeGaussRandom(R);
  R *= tauinv;
  R += (*it)->R + (*it)->Drift;
}

void MCWalkerConfiguration::clear() {
  WalkerList.clear();
}


void MCWalkerConfiguration::copy(iterator first, iterator last){
  while(first != last) {
    WalkerList.push_back(*first); first++;
  }
}

void MCWalkerConfiguration::reset() {
  iterator it=WalkerList.begin();
  while(it != WalkerList.end()) {(*it)->reset();it++;}
}

bool MCWalkerConfiguration::createAuxDataSet(int nfield) {

  if(DataSet.size()) return false;

  for(int iw=0; iw<WalkerList.size(); iw++) {
    DataSet.push_back(new WalkerData_t);
    DataSet[iw]->reserve(nfield);
  }

  return true;
}

void MCWalkerConfiguration::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
  R = awalker.R;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->evaluate(*this);
    DistTables[i]->registerData(buf);
  }
}
  
void MCWalkerConfiguration::copyToBuffer(PooledData<RealType>& buf) {
  //cout << "ParticleSet::copyToBuffer " << endl;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->copyToBuffer(buf);
    //DistTables[i]->print(cout);
  }
}

void MCWalkerConfiguration::copyFromBuffer(PooledData<RealType>& buf) {
  //cout << "ParticleSet::copyFromBuffer " << endl;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->copyFromBuffer(buf);
    //DistTables[i]->print(cout);
  }
}


int MCWalkerConfiguration::branch(int maxcopy, int Nmax, int Nmin) {

  iterator it = WalkerList.begin();
  int iw=0, nw = WalkerList.size();

  vector<Walker_t*> good, bad;
  vector<int> ncopy;
  ncopy.reserve(nw);

  //Temporary for vectorization: any line with energy or Energy is commented out
  //since the performance is not improved by using the vectroized container.
  // vector<RealType> energy;
  // energy.reserve(Energy.size());
  // int ncols = Energy.cols();
  int num_walkers=0;
  while(it != WalkerList.end()) {
    int nc = std::min(static_cast<int>((*it)->Properties(MULTIPLICITY)),maxcopy);
    if(nc) {
      num_walkers += nc;
      good.push_back(*it);
      ncopy.push_back(nc-1);
      //  energy.insert(energy.end(), Energy[iw], Energy[iw]+ncols);
    } else {
      bad.push_back(*it);
    }
    iw++;it++;
  }


  //remove bad walkers
  for(int i=0; i<bad.size(); i++) delete bad[i];

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
      cerr << "Too few walkers to copy! Abort." << endl;
      exit(-1);
    } else {
      int i=0;
      while(i<ncopy.size() && nadd<nadd_target) {
	ncopy[i]++; nadd++;i++;
      }
    }
    num_walkers +=  nadd;
  }

  //clear the WalkerList to populate them with the good walkers
  WalkerList.clear();
  WalkerList.insert(WalkerList.begin(), good.begin(), good.end());

  //  Energy.resize(num_walkers,ncols);
  //  std::copy(energy.begin(), energy.end(), Energy.data());
  //  vector<RealType>::iterator ie=energy.begin();
  //  RealType *e_start = Energy.data()+energy.size();
  int cur_walker = good.size();
  for(int i=0; i<good.size(); i++) { //,ie+=ncols) {
    for(int j=0; j<ncopy[i]; j++, cur_walker++) {
      WalkerList.push_back(new Walker_t(*(good[i])));
      //std::copy(ie,ie+ncols,e_start);
      //e_start += ncols;
    }
  }

  iw=0;
  it=WalkerList.begin();
  while(it != WalkerList.end()) {
    (*it)->Properties(WEIGHT) = 1.0;
    (*it)->Properties(MULTIPLICITY) = 1.0;
    it++;
  }

  return iw;
}

void MCWalkerConfiguration::setUpdateMode(int updatemode) { 

  ///get the tables for which this particle set is the visitor or source
  if(DistTables.empty())  DistanceTable::getTables(ObjectTag,DistTables);

  UpdateMode = updatemode;
  if(UpdateMode == Update_All) {
    DistanceTable::create(WalkerList.size());
  } else {
    DistanceTable::create(1);
  }
}

void MCWalkerConfiguration::loadWalker(Walker_t& awalker) {
  R = awalker.R;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->evaluate(*this);
  }
}

/** move a particle iat
 *@param iat the index of the particle to be moved
 *@param newpos the new position
 *
 *Update activePtcl index and activePos position for the proposed move.
 *Evaluate the related distance table data DistanceTableData::Temp.
 */
void MCWalkerConfiguration::makeMove(int iat, const SingleParticlePos_t& newpos) {
  activePtcl=iat;
  activePos=R[iat]+newpos;
  for(int i=0; i< DistTables.size(); i++) {
    DistTables[i]->move(*this,activePos,iat);
  }
}

/** update the particle attribute by the proposed move
 *@param iat the particle index
 *
 *When the activePtcl is equal to iat, overwrite the position and update the
 *content of the distance tables.
 */
void MCWalkerConfiguration::acceptMove(int iat) {
  if(iat == activePtcl) {
    R[iat]=activePos; 
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->update(iat);
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
