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
#ifndef OHMMS_QMC_HAMILTONIAN_H
#define OHMMS_QMC_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  /**  Collection of Local Energy Operators 
   *
   * Note that QMCHamiltonian is not derived from QMCHmailtonianBase.
   */
  class QMCHamiltonian  {

  public:
  
    typedef QMCHamiltonianBase::RealType  RealType;
    typedef QMCHamiltonianBase::ValueType ValueType;
    typedef QMCHamiltonianBase::Return_t  Return_t;
    
    ///constructor
    QMCHamiltonian();
    ///destructor
    ~QMCHamiltonian();

    void add(QMCHamiltonianBase* h, const string& aname);
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

    /** returnn QMCHamiltonian with the name aname
     * @param aname name of a QMCHamiltonian
     * @return 0 if aname is not found.
     */
    QMCHamiltonianBase* getHamiltonian(const string& aname);
    ////return the LocalEnergy \f$=\sum_i H^{qmc}_{i}\f$
    inline Return_t getLocalEnergy() { return LocalEnergy;}
    ////return the LocalPotential \f$=\sum_i H^{qmc}_{i} - KE\f$
    inline Return_t getLocalPotential() { return LocalEnergy-Hvalue[0];}

    /** set Tau for each Hamiltonian
     */
    inline void setTau(RealType tau) {
      for(int i=0; i< H.size();++i)H[i]->setTau(tau); 
    }
    
    /** evaluate Local Energy
     * @param P ParticleSet
     * @return the local energy
     *
     * P.R, P.G and P.L are used to evaluate the LocalEnergy.
     */
    Return_t evaluate(ParticleSet& P);

   private:

    ///Current Local Energy
    Return_t LocalEnergy;
    ///vector of Hamiltonians
    std::vector<QMCHamiltonianBase*> H;
    ///vector containing the index of the Hamiltonians
    std::vector<int> Hindex;
    ///vector containing the values of the Hamiltonians
    std::vector<Return_t> Hvalue;
    ///vector containing the names of the Hamiltonians
    std::vector<string> Hname;
    ///map the name to an index
    std::map<string,int> Hmap;
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

