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
using namespace ohmmsqmc;
#include "ParticleBase/RandomSeqGenerator.h"

MCWalkerConfiguration::MCWalkerConfiguration(): 
UpdateMode(Update_Walker) {

}


///default destructor
MCWalkerConfiguration::~MCWalkerConfiguration(){

  destroyWalker(WalkerList.begin(), WalkerList.end());
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
  while(it != WalkerList.end()) {(*it)->Properties(Weight) = 1.0;it++;}
}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
