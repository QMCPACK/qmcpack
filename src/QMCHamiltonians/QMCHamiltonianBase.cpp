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

/*!
 *\param h the Hamiltonian
 *\param aname the name of the Hamiltonian
 *\brief add a new Hamiltonian the the list of Hamiltonians.
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

/*!
 *\param P input configuration containing N particles
 *\return the local energy
 *\brief Evaluate all the Hamiltonians for the N-particle
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

/*!
 *\param W a set of walkers (N-particle configurations)
 *\param LE return a vector containing the local 
 energy for each walker
 *\brief Evaluate all the Hamiltonians for a set of N-particle
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

