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

  /** An abstract class for Local Energy operators */

  struct QMCHamiltonianBase: public QMCTraits {

    typedef ParticleAttrib<ValueType>  ValueVectorType;

    ///constructor
    QMCHamiltonianBase(){}

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
  };

  /**  Collection of Local Energy Operators */

  class QMCHamiltonian {

  public:

    typedef QMCHamiltonianBase::RealType        RealType;
    typedef QMCHamiltonianBase::ValueType       ValueType;
    typedef QMCHamiltonianBase::ValueVectorType ValueVectorType;

    ///constructor
    QMCHamiltonian();
    ///destructor
    ~QMCHamiltonian();

    void add(QMCHamiltonianBase* h, const string& aname);

    ///return the name of ith Hamiltonian 
    inline string getName(int i) const { return Hname[i];}

    ///return the value of Hamiltonian i
    inline RealType operator[](int i) const { return Hvalue[i];}

    ///return the number of Hamiltonians
    inline int size() const { return H.size();}

    ///assign the Hamiltonian values to a vector a
    template<class V>
    inline void get(V& a) const {
      std::copy(Hvalue.begin(), Hvalue.end(), a.begin());
    }

    QMCHamiltonianBase* getHamiltonian(const string& aname);
    
    RealType getLocalPotential();
    
    ValueType evaluate(ParticleSet& P);

    void evaluate(WalkerSetRef& W, ValueVectorType& LE);

   private:
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

