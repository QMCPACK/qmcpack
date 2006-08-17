/////////////////////////////////////////////////////////////////
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
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"
using namespace qmcplusplus;

/** constructor
 */
QMCHamiltonian::QMCHamiltonian(){ }

/// copy constructor is distable by declaring it as private
QMCHamiltonian::QMCHamiltonian(const QMCHamiltonian& qh) {}

/** destructor
 */
QMCHamiltonian::~QMCHamiltonian() {
  
  DEBUGMSG("QMCHamiltonian::~QMCHamiltonian")
    
}

bool QMCHamiltonian::get(std::ostream& os) const {
  for(int i=0; i<H.size(); i++) {
    os << "    "; H[i]->get(os); os << "\n";
  }
  return true;
}

/** add a new Hamiltonian the the list of Hamiltonians.
 *@param h the Hamiltonian
 *@param aname the name of the Hamiltonian
 */
void 
QMCHamiltonian::addOperator(QMCHamiltonianBase* h, const string& aname) {
  //check if already added, if not add at the end
  map<string,int>::iterator it = Hmap.find(aname);
  if(it == Hmap.end()) {
    Hmap[aname] = H.size();
    Hname.push_back(aname);
    H.push_back(h);
  }
  Hvalue.resize(H.size(),RealType());
}

/** remove a named Hamiltonian from the list
 *@param aname the name of the Hamiltonian
 *@return true, if the request hamiltonian exists and is removed.
 */
bool 
QMCHamiltonian::remove(const string& aname) {
  map<string,int>::iterator it = Hmap.find(aname);
  if(it != Hmap.end()) {
    int n = (*it).second;
    Hmap.erase(aname); 
    map<string,int>::iterator jt = Hmap.begin();
    while(jt != Hmap.end()) {
      if((*jt).second > n) (*jt).second -= 1;
      jt++;
    }
    delete H[n];
    for(int i=n+1; i<H.size(); i++) H[i-1]=H[i];
    H.pop_back();
    Hvalue.resize(H.size(),RealType());
    return true;
  }
  return false;
}

/** add a number of properties to the ParticleSet
 * @param P ParticleSet to which multiple columns to be added
 * 
 * QMCHamiltonian can add any number of properties to a ParticleSet.
 * Hindex contains the index map to the ParticleSet::PropertyList.
 * This enables assigning the properties evaluated by each QMCHamiltonianBase
 * object to the correct property column.
 */
void 
QMCHamiltonian::add2WalkerProperty(ParticleSet& P) {
  Hindex.resize(Hname.size());
  for(int i=0; i<Hname.size(); i++) Hindex[i]=P.addProperty(Hname[i]);
}

/** Evaluate all the Hamiltonians for the N-particle  configuration
 *@param P input configuration containing N particles
 *@return the local energy
 */
QMCHamiltonian::Return_t 
QMCHamiltonian::evaluate(ParticleSet& P) {
  LocalEnergy = 0.0;
  vector<QMCHamiltonianBase*>::iterator hit(H.begin()),hit_end(H.end());
  int i(0);
  while(hit != hit_end) {
    LocalEnergy += (*hit)->evaluate(P);
    Hvalue[i]=(*hit)->Value; 
    ++hit;++i;
  }
  return LocalEnergy;
}

QMCHamiltonian::Return_t 
QMCHamiltonian::evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
  LocalEnergy = 0.0;
  vector<QMCHamiltonianBase*>::iterator hit(H.begin()),hit_end(H.end());
  int i(0);
  while(hit != hit_end) {
    LocalEnergy += (*hit)->evaluate(P,Txy);
    Hvalue[i]=(*hit)->Value; 
    ++hit;++i;
  }
  return LocalEnergy;
}


QMCHamiltonian::Return_t 
QMCHamiltonian::getEnsembleAverage() {
  Return_t sum=0.0;
  for(int i=0; i<H.size(); i++) sum += H[i]->getEnsembleAverage();
  return sum;
}
/** return pointer to the QMCHamtiltonian with the name
 *@param aname the name of Hamiltonian
 *@return the pointer to the named term.
 *
 * If not found, return 0
 */
QMCHamiltonianBase* 
QMCHamiltonian::getHamiltonian(const string& aname) {

  //check if already added, if not add at the end
  map<string,int>::iterator it = Hmap.find(aname);
  if(it == Hmap.end()) 
    return 0;
  else 
    return H[(*it).second];
}

void
QMCHamiltonian::resetTargetParticleSet(ParticleSet& P) {
  for(int i=0; i<H.size(); i++) H[i]->resetTargetParticleSet(P);
}

//Meant for vectorized operators. Not so helpful.
///**
// *@param W a set of walkers (N-particle configurations)
// *@param LE return a vector containing the local 
// energy for each walker
// *@brief Evaluate all the Hamiltonians for a set of N-particle
// *configurations
// */
//void QMCHamiltonian::evaluate(WalkerSetRef& W, ValueVectorType& LE) {
//  LE = 0.0;
//  for(int i=0; i<H.size(); i++) H[i]->evaluate(W, LE);
//}

//bool QMCHamiltonian::put(std::istream& ) {
//  return true;
//}
//
//bool QMCHamiltonian::put(xmlNodePtr cur) {
//  return true;
//}
//
//void QMCHamiltonian::reset() { }

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

