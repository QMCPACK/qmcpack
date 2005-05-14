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
#ifndef OHMMS_QMC_WALKER_H
#define OHMMS_QMC_WALKER_H

#include "OhmmsPETE/OhmmsMatrix.h"

namespace ohmmsqmc {

  /** an enum denoting index of physical properties */
  enum {WEIGHT=0,     /*!< weight */
	LOCALENERGY,  /*!< local energy, the sum of all the components */
	LOCALPOTENTIAL, /*!< local potential energy = local energy - kinetic energy */
	MULTIPLICITY, /*!< multiplicity, used by DMC for branching */
	LOGPSI,        /*!< log(fabs(psi)) instead of square of the many-body wavefunction \f$|\Psi|^2\f$ */
        SUMRATIO,      /*!< sum of wavefunction ratios for multiple H and Psi */
	SIGN,          /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
	AGE,          /*!< the age of the walker. set to zero when the walker is updated */
	WOSVAR,       /*!< Variance of WOS potential */
	NUMPROPERTIES /*!< the number of properties */
	//SCALED,       /*!< scaling factor for the drift */
	//CAPACITY=15
       };
  
  /**
   *\brief A container class to represent a walker.
   *
   * A walker stores the particle configurations 
   * and a property container.
   * The template (P)articleSet(A)ttribute is a generic container
   * of position types.
   * The template (G)radient(A)ttribute is a generic container
   * of gradients types.
   */
  template<class T, class PA, class GA=PA>
  struct Walker {
    
    ///typedef for the property container, fixed size
    typedef TinyVector<T,NUMPROPERTIES> PropertyContainer_t;

    ///id reserved for forward walking
    int ID;

    ///scalar properties of a walker
    PropertyContainer_t  Properties;

    /**the configuration vector (3N-dimensional vector to store
       the positions of all the particles for a single walker)*/
    PA R;
    
    ///drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
    GA Drift;

    ///dynamic properties in addition to Properties
    Matrix<T> DynProperty;

    ///create a walker for n-particles
    inline explicit Walker(int n) : Properties(0.0) {  
      Properties[WEIGHT] = 1.0;
      Properties[MULTIPLICITY] = 1.0;
      Properties[SIGN] = 1.0;
      resize(n);
    }
    
    ///constructor
    inline Walker() : Properties(0.0) {
      Properties[WEIGHT] = 1.0;
      Properties[MULTIPLICITY] = 1.0;
      Properties[SIGN] = 1.0;
    }

    ///copy constructor
    inline Walker(const Walker& a){
      makeCopy(a);
    }
    
    ///assignment operator
    inline Walker& operator=(const Walker& a) {
      makeCopy(a);
      return *this;
    }

    ///return the number of particles per walker
    inline int size() const { return R.size(); }

    ///resize for n-particles
    inline void resize(int nptcl) {
      R.resize(nptcl); Drift.resize(nptcl); 
    }

    inline void makeCopy(const Walker& a) {    
      resize(a.R.size());
      R = a.R;
      Drift = a.Drift;
      Properties = a.Properties;
      DynProperty.copy(a.DynProperty);
    }

    //return the address of the values of Hamiltonian terms
    inline T* restrict getEnergyBase() {
      return DynProperty.data();
    }

    //return the address of the values of Hamiltonian terms
    inline const T* restrict getEnergyBase() const {
      return DynProperty.data();
    }

    inline T* restrict getEnergyBase(int i) {
      return DynProperty[i];
    }

    //return the address of the values of Hamiltonian terms
    inline const T* restrict getEnergyBase(int i) const {
      return DynProperty[i];
    }

    /** reset the property of a walker
     *@param iw the id of this walker
     *@param psi the wavefunction value
     *@param ene the local energy
     *
     *Assign the values and reset the weight and multiplicity to one to start a run
     */
    inline void resetProperty(T logpsi, T sigN, T ene) {
      Properties(WEIGHT) = 1.0;
      Properties(MULTIPLICITY) = 1.0;
      Properties(LOCALENERGY) = ene;
      Properties(LOGPSI)=logpsi;
      Properties(SIGN)=sigN;
    }

    /** reset the walker weight, multiplicity and age */
    inline void reset() {
      Properties(WEIGHT)=1.0;
      Properties(MULTIPLICITY)=1.0;
      Properties(AGE)=0.0;
    }

    inline void resizeDynProperty(int n, int m) {
      DynProperty.resize(n,m);
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
