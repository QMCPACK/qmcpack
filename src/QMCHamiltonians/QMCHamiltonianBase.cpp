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
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
using namespace ohmmsqmc;

QMCHamiltonian::QMCHamiltonian() { }

QMCHamiltonian::~QMCHamiltonian() {
  
  DEBUGMSG("QMCHamiltonian::~QMCHamiltonian")
    
    }

/** add a new Hamiltonian the the list of Hamiltonians.
 *@param h the Hamiltonian
 *@param aname the name of the Hamiltonian
 */
void 
QMCHamiltonian::add(QMCHamiltonianBase* h, const string& aname) {
  //check if already added, if not add at the end
  map<string,int>::iterator it = Hmap.find(aname);
  if(it == Hmap.end()) {
    Hmap[aname] = H.size();
    Hname.push_back(aname);
    H.push_back(h);
  }
  Hvalue.resize(H.size()+1,RealType());
}

/** remove a named Hamiltonian from the list
 *@param aname the name of the Hamiltonian
 *@param return true, if the request hamiltonian exists and is removed.
 */
bool 
QMCHamiltonian::remove(const string& aname) {

  map<string,int>::iterator it = Hmap.find(aname);
  if(it != Hmap.end()) {
    int n = (*it).second;
    it++;
    while(it != Hmap.end()) {
      (*it).second--; it++;
    }
    Hmap.erase(aname); 
    delete H[n];
    for(int i=n+1; i<H.size(); i++) H[i-1]=H[i];
    H.pop_back();
    Hvalue.resize(H.size()+1,RealType());
    return true;
  }

  return false;
}

/**
 *@param P input configuration containing N particles
 *@return the local energy
 *@brief Evaluate all the Hamiltonians for the N-particle
 *configuration
 */
QMCHamiltonian::ValueType 
QMCHamiltonian::evaluate(ParticleSet& P) {

  ValueType esum = 0.0;
  register int i=0;
  for(; i<H.size(); i++)  {
    esum += H[i]->evaluate(P,Hvalue[i]);
  }
  return Hvalue[i]=esum;
}

/**
 *@param aname the name of Hamiltonian
 *@return the pointer to the named term.
 *@brief Get a QMCHamiltonian with the name.
 *@note If not found, return NULL
 */
QMCHamiltonianBase* 
QMCHamiltonian::getHamiltonian(const string& aname) {

  //check if already added, if not add at the end
  map<string,int>::iterator it = Hmap.find(aname);
  if(it == Hmap.end()) 
    return NULL;
  else 
    return H[(*it).second];
}

QMCHamiltonian::RealType 
QMCHamiltonian::getLocalPotential() {
  
  return Hvalue.back()-Hvalue[0];
}

/**
 *@param W a set of walkers (N-particle configurations)
 *@param LE return a vector containing the local 
 energy for each walker
 *@brief Evaluate all the Hamiltonians for a set of N-particle
 *configurations
 */
void QMCHamiltonian::evaluate(WalkerSetRef& W, ValueVectorType& LE) {
  LE = 0.0;
  for(int i=0; i<H.size(); i++) H[i]->evaluate(W, LE);
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

