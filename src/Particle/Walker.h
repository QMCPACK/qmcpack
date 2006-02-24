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
#ifndef QMCPLUSPLUS_WALKER_H
#define QMCPLUSPLUS_WALKER_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus {


  /** an enum denoting index of physical properties */
  enum {LOGPSI=0,       /*!< log(fabs(psi)) instead of square of the many-body wavefunction \f$|\Psi|^2\f$ */
	SIGN,           /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
        UMBRELLAWEIGHT, /*!< sum of wavefunction ratios for multiple H and Psi */
	LOCALENERGY,    /*!< local energy, the sum of all the components */
	LOCALPOTENTIAL, /*!< local potential energy = local energy - kinetic energy */
	NUMPROPERTIES   /*!< the number of properties */
       };
  
  /** A container class to represent a walker.
   *
   * A walker stores the particle configurations {R}  and a property container.
   * The template (P)articleSet(A)ttribute is a generic container  of position types.
   * The template (G)radient(A)ttribute is a generic container of gradients types.
   * Data members for each walker
   * - ID : identity for a walker. default is 0. 
   * - Age : generation after a move is accepted.
   * - Weight : weight to take the ensemble averages
   * - Multiplicity : multiplicity for branching. Probably can be removed.
   * - Properties  : 2D container. The first index corresponds to the H/Psi index and second index >=NUMPROPERTIES.
   * - DataSet : anonymous container. 
   */
  template<class T, class PA, class GA=PA>
  struct Walker {
    
    ///typedef for the property container, fixed size
    typedef Matrix<T>      PropertyContainer_t;
    typedef PooledData<T>  Buffer_t;

    ///id reserved for forward walking
    int ID;
    ///Age of this walker.
    int Age;

    ///Weight of the walker
    T Weight;

    /** Number of copies for branching
     *
     * When Multiplicity = 0, this walker will be destroyed.
     */
    T Multiplicity;

    ///scalar properties of a walker
    PropertyContainer_t  Properties;

    /**the configuration vector (3N-dimensional vector to store
       the positions of all the particles for a single walker)*/
    PA R;
    
    ///drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
    GA Drift;

    ///buffer for the data for particle-by-particle update
    Buffer_t DataSet;

    ///default constructor
    inline Walker() : Age(0),Weight(1.0e0),Multiplicity(1.0e0) {
      Properties.resize(1,NUMPROPERTIES);
      reset();
    }

    ///create a walker for n-particles
    inline explicit Walker(int nptcl) : Age(0), Weight(1.0e0),Multiplicity(1.0e0){  
      Properties.resize(1,NUMPROPERTIES);
      resize(nptcl);
      reset();
    }

    ///copy constructor
    inline Walker(const Walker& a):Age(0),Weight(1.0e0), Multiplicity(1.0e0){
      makeCopy(a);
    }

    inline ~Walker() { }
    
    ///assignment operator
    inline Walker& operator=(const Walker& a) {
      makeCopy(a);
      return *this;
    }

    inline void assign(const Walker& a) {
      Age=a.Age;
      Weight=a.Weight;
      Multiplicity=a.Multiplicity;
      R = a.R;
      Drift = a.Drift;
      Properties=a.Properties;
      if(a.DataSet.size()) {
        DataSet=a.DataSet;
      }
    }

    ///return the number of particles per walker
    inline int size() const { return R.size(); }

    ///resize for n-particles
    inline void resize(int nptcl) {
      R.resize(nptcl); Drift.resize(nptcl); 
    }

    ///copy the content of a walker
    inline void makeCopy(const Walker& a) {    
      resize(a.R.size());
      Age=a.Age;
      Multiplicity=a.Multiplicity;
      R = a.R;
      Drift = a.Drift;
      Properties.copy(a.Properties);
      if(a.DataSet.size()) {
        DataSet=a.DataSet;
      }
    }

    //return the address of the values of Hamiltonian terms
    inline T* restrict getPropertyBase() {
      return Properties.data();
    }

    //return the address of the values of Hamiltonian terms
    inline const T* restrict getPropertyBase() const {
      return Properties.data();
    }

    ///return the address of the i-th properties
    inline T* restrict getPropertyBase(int i) {
      return Properties[i];
    }

    ///return the address of the i-th properties
    inline const T* restrict getPropertyBase(int i) const {
      return Properties[i];
    }

    /** reset the property of a walker
     *@param logpsi \f$\log |\Psi|\f$
     *@param sigN  sign of the trial wavefunction
     *@param ene the local energy
     *
     *Assign the values and reset the age
     * but leave the weight and multiplicity 
     */
    inline void resetProperty(T logpsi, T sigN, T ene) {
      Age=0;
      //Weight=1.0;
      Properties(LOGPSI)=logpsi;
      Properties(SIGN)=sigN;
      Properties(LOCALENERGY) = ene;
    }

    /** marked to die
     *
     * Multiplicity and weight are set to zero.
     */
    inline void willDie() {
      Multiplicity=0;
      Weight=0.0;
    }

    /** reset the walker weight, multiplicity and age */
    inline void reset() {
      Age=0;
      Multiplicity=1.0e0;
      Weight=1.0e0;
    }

    inline void resizeProperty(int n, int m) {
      Properties.resize(n,m);
    }


    /** byte size for a packed message 
     *
     * ID, Age, Properties, R, Drift, DataSet is packed
     */
    inline int byteSize() {
      return 2*sizeof(int)+(Properties.size()+OHMMS_DIM*2*R.size()+DataSet.size())*sizeof(T);
    }

    template<class Msg>
    inline Msg& putMessage(Msg& m) {
      int nat=R.size();
      m << ID << Age;
      for(int iat=0; iat<nat;iat++) R[iat].putMessage(m);
      for(int iat=0; iat<nat;iat++) Drift[iat].putMessage(m);
      Properties.putMessage(m);
      DataSet.putMessage(m);
      return m;
    }

    template<class Msg>
    inline Msg& getMessage(Msg& m) {
      int nat=R.size();
      m >> ID >> Age;
      for(int iat=0; iat<nat;iat++) R[iat].getMessage(m);
      for(int iat=0; iat<nat;iat++) Drift[iat].getMessage(m);
      Properties.getMessage(m);
      DataSet.getMessage(m);
      return m;
    }

  };

  template<class T, class PA>
    ostream& operator<<(ostream& out, const Walker<T,PA>& rhs)
    {
      copy(rhs.Properties.begin(), rhs.Properties.end(), 
      	   ostream_iterator<double>(out," "));
      out << endl;
      out << rhs.R;
      return out;
    }
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
