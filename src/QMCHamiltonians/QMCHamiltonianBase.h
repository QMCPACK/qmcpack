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
#ifndef OHMMS_QMC_HAMILTONIANBASE_H
#define OHMMS_QMC_HAMILTONIANBASE_H
#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

/**@file QMCHamiltonianBase.h
 *@brief Declaration of QMCHamiltonianBase and QMCHamiltonian
 */
namespace ohmmsqmc {

  class ParticeBase;
  class WalkerSetRef;
  class DistanceTableData;

  /** An abstract class for Local Energy operators 
   *
   * Use of ValueType is questionable. The types should be checked when using
   * complex wave functions.
   */ 
  struct QMCHamiltonianBase: public QMCTraits {

    RealType Tau;
    RealType Value;

    typedef ParticleAttrib<ValueType>  ValueVectorType;

    ///constructor
    QMCHamiltonianBase():Tau(0.0),Value(0.0){}

    ///virtual destructor
    virtual ~QMCHamiltonianBase() { }

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@return the value of the Hamiltonian
     */
    virtual ValueType evaluate(ParticleSet& P) = 0; 

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@param x the sum of local energies
     *@return the value the Local Energy
    */
    virtual ValueType evaluate(ParticleSet& P, RealType& x) = 0;

    /** Evaluate the local energies of the entire walkers
     *@param W a set of walkers (N-particle configurations)
     *@param LE return a vector containing the value
     */
    virtual 
    void evaluate(WalkerSetRef& W, ValueVectorType& LE) = 0;

    inline void setTau(RealType tau) { Tau = tau;}
  };

  /**  Collection of Local Energy Operators 
   *
   * Note that QMCHamiltonian is not derived from QMCHmailtonianBase.
   */
  class QMCHamiltonian  {

  public:

    typedef QMCHamiltonianBase::RealType        RealType;
    typedef QMCHamiltonianBase::ValueType       ValueType;
    typedef QMCHamiltonianBase::ValueVectorType ValueVectorType;

    ///constructor
    QMCHamiltonian();
    ///destructor
    ~QMCHamiltonian();

    void add(QMCHamiltonianBase* h, const string& aname);
    bool remove(const string& aname);

    ///return the name of ith Hamiltonian 
    inline string getName(int i) const { return Hname[i];}

    ///return the value of Hamiltonian i
    inline RealType operator[](int i) const { return Hvalue[i];}

    ///return the number of Hamiltonians
    inline int size() const { return H.size();}

    ///copy the local content to iterator first
    template<class IT>
    inline void copy(IT first) { std::copy(Hvalue.begin(), Hvalue.end(), first);}

    template<class IT>
    inline void copy(IT first, RealType wgt) { 
      vector<RealType>::iterator it(Hvalue.begin()),it_end(Hvalue.end());
      while(it != it_end) {*first++ = wgt*(*it++);}
      //for(int i=0; i<Hvalue.size(); i++,first++) *first = wgt*Hvalue[i];
    }

    QMCHamiltonianBase* getHamiltonian(const string& aname);
    
    inline RealType getLocalEnergy() { return LocalEnergy;}
    inline RealType getLocalPotential() { return LocalEnergy-Hvalue[0];}

    inline void setTau(RealType tau) {for(int i=0; i< H.size();++i)H[i]->setTau(tau); }
    
    ValueType evaluate(ParticleSet& P);

    void evaluate(WalkerSetRef& W, ValueVectorType& LE);

   private:
    ///Current Local Energy
    RealType LocalEnergy;
    ///vector of Hamiltonians
    std::vector<QMCHamiltonianBase*> H;
    ///map the name to an index
    std::map<string,int> Hmap;
    ///vector containing the names of the Hamiltonians
    std::vector<string> Hname;
    ///vector containing the values of the Hamiltonians
    std::vector<RealType> Hvalue;
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

