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
/**@file QMCHamiltonian.h
 *@brief Declaration of QMCHamiltonian
 */
#ifndef QMCPLUSPLUS_HAMILTONIAN_H
#define QMCPLUSPLUS_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {

  class NewTimer;
  /**  Collection of Local Energy Operators 
   *
   * Note that QMCHamiltonian is not derived from QMCHmailtonianBase.
   */
  class QMCHamiltonian {

  public:
  
    typedef QMCHamiltonianBase::RealType  RealType;
    typedef QMCHamiltonianBase::ValueType ValueType;
    typedef QMCHamiltonianBase::Return_t  Return_t;
    
    ///constructor
    QMCHamiltonian();

    ///copy constructor

    ///destructor
    ~QMCHamiltonian();

    void addOperator(QMCHamiltonianBase* h, const string& aname);
    bool remove(const string& aname);

    ///add each term to the PropertyList for averages
    void add2WalkerProperty(ParticleSet& P);

    ///retrun the starting index
    inline int startIndex() const { return Hindex[0];}

    ///return the name of ith Hamiltonian 
    inline string getName(int i) const { return Hname[i];}

    ///return the value of Hamiltonian i
    inline Return_t operator[](int i) const { return Hvalue[i];}

    ///return the number of Hamiltonians
    inline int size() const { return H.size();}

    ///save the values of Hamiltonian elements to the Properties
    template<class IT>
    inline 
    void saveProperty(IT first) {
      first[LOCALPOTENTIAL]= LocalEnergy-Hvalue[0];
      std::copy(Hvalue.begin(),Hvalue.end(),first+Hindex[0]);
    }

    /** return QMCHamiltonianBase with the name aname
     * @param aname name of a QMCHamiltonianBase
     * @return 0 if aname is not found.
     */
    QMCHamiltonianBase* getHamiltonian(const string& aname);

    /** return i-th QMCHamiltonianBase
     * @param i index of the QMCHamiltonianBase
     * @return H[i]
     */
    QMCHamiltonianBase* getHamiltonian(int i) {
      return H[i];
    }

    ////return the LocalEnergy \f$=\sum_i H^{qmc}_{i}\f$
    inline Return_t getLocalEnergy() { return LocalEnergy;}
    ////return the LocalPotential \f$=\sum_i H^{qmc}_{i} - KE\f$
    inline Return_t getLocalPotential() { return LocalEnergy-Hvalue[0];}
    ///return the energy that does not depend on variational parameters
    inline Return_t getInvariantEnergy() const { 
      Return_t s=0;
      for(int i=0; i<H.size(); i++) {
        if(!H[i]->UpdateMode[QMCHamiltonianBase::OPTIMIZABLE]) s+=Hvalue[i];
      }
      return s;
    }

    /** set Tau for each Hamiltonian
     */
    inline void setTau(RealType tau) {
      for(int i=0; i< H.size();i++)H[i]->setTau(tau); 
    }

    /** set Tau for each Hamiltonian
     */
    inline void setPrimary(bool primary) {
      for(int i=0; i< H.size();i++) 
        H[i]->UpdateMode.set(QMCHamiltonianBase::PRIMARY,1);
    }
    
    /** return if WaveFunction Ratio needs to be evaluated
     *
     * This is added to handle orbital-dependent QMCHamiltonianBase during
     * orbital optimizations.
     */
    inline bool needRatio() {
      bool dependOnOrbital=false;
      for(int i=0; i< H.size();i++)  
        if(H[i]->UpdateMode[QMCHamiltonianBase::RATIOUPDATE]) dependOnOrbital=true;
      return dependOnOrbital;
    }

    /** evaluate Local Energy
     * @param P ParticleSet
     * @return the local energy
     *
     * P.R, P.G and P.L are used to evaluate the LocalEnergy.
     */
    Return_t evaluate(ParticleSet& P);


    /** evaluate Local and NonLocal energies
     * @param P ParticleSEt
     * @param Txy transition matrix of nonlocal Hamiltonians
     * @return Local energy
     */
    Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy); 

    /** return an average value of the LocalEnergy 
     *
     * Introduced to get a collective value
     */ 
    Return_t getEnsembleAverage();

    void resetTargetParticleSet(ParticleSet& P);

    /** By mistake, QMCHamiltonian::getName(int i) is used
     * and this is in conflict with the declaration of OhmmsElementBase.
     * For the moment, QMCHamiltonian is not inherited from OhmmsElementBase.
     */
    void setName(const string& aname) { 
      myName=aname;
    }

    
    string getName() const { return myName;}

    bool get(std::ostream& os) const;

    void setRandomGenerator(RandomGenerator_t* rng);

    /** return a clone */
    QMCHamiltonian* makeClone(ParticleSet& qp, TrialWaveFunction& psi); 

   private:

    ///Current Local Energy
    Return_t LocalEnergy;
    ///getName is in the way
    string myName;
    ///vector of Hamiltonians
    std::vector<QMCHamiltonianBase*> H;
    ///timers
    std::vector<NewTimer*> myTimers;
    ///vector containing the index of the Hamiltonians
    std::vector<int> Hindex;
    ///vector containing the values of the Hamiltonians
    std::vector<Return_t> Hvalue;
    ///vector containing the names of the Hamiltonians
    std::vector<string> Hname;
    ///map the name to an index
    std::map<string,int> Hmap;

    /////disable copy constructor
    //QMCHamiltonian(const QMCHamiltonian& qh);
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

