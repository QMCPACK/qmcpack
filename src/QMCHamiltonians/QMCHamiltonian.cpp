/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/NewTimer.h"
namespace qmcplusplus
{

  /** constructor
  */
  QMCHamiltonian::QMCHamiltonian():myIndex(0)
  { }

///// copy constructor is distable by declaring it as private
//QMCHamiltonian::QMCHamiltonian(const QMCHamiltonian& qh) {}

/** destructor
 */
QMCHamiltonian::~QMCHamiltonian() 
{
  //@todo clean up H and auxH
}

bool QMCHamiltonian::get(std::ostream& os) const 
{
  for(int i=0; i<H.size(); i++) {
    os.setf(ios::left);
    os << "  " << setw(16) << H[i]->myName; 
    H[i]->get(os); os << "\n";
  }
  return true;
}

/** add a new Hamiltonian the the list of Hamiltonians.
 * @param h an operator
 * @param aname name of h
 * @param physical if true, a physical operator
 */
void 
QMCHamiltonian::addOperator(QMCHamiltonianBase* h, const string& aname, bool physical) 
{
  if(physical)
  {
    for(int i=0; i<H.size(); ++i)
    {
      if(H[i]->myName == aname) 
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to H, physical Hamiltonian " << endl;
    h->myName=aname;
    H.push_back(h);
    string tname="Hamiltonian::"+aname;
    NewTimer *atimer=new NewTimer(tname);
    myTimers.push_back(atimer);
    TimerManager.addTimer(atimer);
  }
  else
  {//ignore timers for now
    for(int i=0; i<auxH.size(); ++i)
    {
      if(auxH[i]->myName == aname) 
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << endl;
        return;
      }
    }

    app_log() << "  QMCHamiltonian::addOperator " << aname << " to auxH " << endl;
    h->myName=aname;
    auxH.push_back(h);
  }
}

///** remove a named Hamiltonian from the list
// *@param aname the name of the Hamiltonian
// *@return true, if the request hamiltonian exists and is removed.
// */
//bool 
//QMCHamiltonian::remove(const string& aname) 
//{
//  return false;
//}

/** add a number of properties to the ParticleSet
 * @param P ParticleSet to which multiple columns to be added
 * 
 * QMCHamiltonian can add any number of properties to a ParticleSet.
 * Hindex contains the index map to the ParticleSet::PropertyList.
 * This enables assigning the properties evaluated by each QMCHamiltonianBase
 * object to the correct property column.
 */
void 
QMCHamiltonian::addObservables(PropertySetType& plist) 
{
  //first add properties to Observables
  Observables.clear();
  for(int i=0; i<H.size(); ++i) H[i]->addObservables(Observables);
  for(int i=0; i<auxH.size(); ++i) auxH[i]->addObservables(Observables);

  myIndex=plist.add(Observables.Names[0]);
  for(int i=1; i<Observables.size(); ++i) plist.add(Observables.Names[i]);

  app_log() << "  QMCHamiltonian::add2WalkerProperty added " << Observables.size() << " data to PropertyList" << endl;
  app_log() << "    starting Index = " << myIndex << endl;
}

/** Evaluate all the Hamiltonians for the N-particle  configuration
 *@param P input configuration containing N particles
 *@return the local energy
 */
QMCHamiltonian::Return_t 
QMCHamiltonian::evaluate(ParticleSet& P) 
{
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    LocalEnergy += H[i]->evaluate(P);
    H[i]->setObservables(Observables);
    myTimers[i]->stop();
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for(int i=0; i<auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
  }

  return LocalEnergy;
}

QMCHamiltonian::Return_t 
QMCHamiltonian::evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
{
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    LocalEnergy += H[i]->evaluate(P,Txy);
    H[i]->setObservables(Observables);
    myTimers[i]->stop();
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for(int i=0; i<auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
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

  for(int i=0; i<H.size(); ++i)
    if(H[i]->myName == aname) return H[i];

  for(int i=0; i<auxH.size(); ++i)
    if(auxH[i]->myName == aname) return auxH[i];

  return 0;
}

void
QMCHamiltonian::resetTargetParticleSet(ParticleSet& P) 
{
  for(int i=0; i<H.size(); i++) H[i]->resetTargetParticleSet(P);
  for(int i=0; i<auxH.size(); i++) auxH[i]->resetTargetParticleSet(P);
}

void
QMCHamiltonian::setRandomGenerator(RandomGenerator_t* rng) 
{
  for(int i=0; i<H.size(); i++) H[i]->setRandomGenerator(rng);
  for(int i=0; i<auxH.size(); i++) auxH[i]->setRandomGenerator(rng);
}

QMCHamiltonian* QMCHamiltonian::makeClone(ParticleSet& qp, TrialWaveFunction& psi) 
{
  QMCHamiltonian* myclone=new QMCHamiltonian;
  for(int i=0; i<H.size(); ++i)
    myclone->addOperator(H[i]->makeClone(qp,psi),H[i]->myName,true);
  for(int i=0; i<auxH.size(); ++i)
    myclone->addOperator(auxH[i]->makeClone(qp,psi),auxH[i]->myName,false);
  //myclone->addObservables(qp.PropertyList);
  myclone->resetObservables(myIndex);
  return myclone;
}

void 
QMCHamiltonian::resetObservables(int start)
{
  //first add properties to Observables
  Observables.clear();
  for(int i=0; i<H.size(); ++i) H[i]->addObservables(Observables);
  for(int i=0; i<auxH.size(); ++i) auxH[i]->addObservables(Observables);
  myIndex=start;
}
QMCHamiltonian::Return_t QMCHamiltonian::registerData(ParticleSet& P, BufferType& buffer)
{
  LocalEnergy=0.0;
  for(int i=0; i<H.size(); ++i) 
  {
    LocalEnergy+=H[i]->registerData(P,buffer);
    H[i]->setObservables(Observables);
  }
  buffer.add(LocalEnergy);
  return LocalEnergy;
}

void QMCHamiltonian::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  for(int i=0; i<H.size(); ++i) H[i]->copyFromBuffer(P,buffer);
  buffer.get(LocalEnergy);
}

QMCHamiltonian::Return_t QMCHamiltonian::evaluate(ParticleSet& P, BufferType& buffer)
{
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    LocalEnergy += H[i]->Value;
    H[i]->setObservables(Observables);
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for(int i=0; i<auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
  }

  for(int i=0; i<H.size(); ++i) H[i]->copyToBuffer(P,buffer);
  buffer.put(LocalEnergy);
  return LocalEnergy;
}

QMCHamiltonian::Return_t QMCHamiltonian::evaluatePbyP(ParticleSet& P, int active)
{
  NewLocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i) NewLocalEnergy +=H[i]->evaluatePbyP(P,active);
  return NewLocalEnergy;
}
void QMCHamiltonian::acceptMove(int active)
{
  for(int i=0; i<H.size(); ++i) H[i]->acceptMove(active);
  LocalEnergy=NewLocalEnergy;
}

void QMCHamiltonian::rejectMove(int active)
{
  for(int i=0; i<H.size(); ++i) H[i]->rejectMove(active);
}
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

