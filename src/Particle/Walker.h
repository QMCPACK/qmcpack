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
#include "Utilities/PooledData.h"

namespace ohmmsqmc {

  /** an enum denoting index of physical properties */
  enum {WEIGHT=0,     /*!< weight */
	LOCALENERGY,  /*!< local energy, the sum of all the components */
	MULTIPLICITY, /*!< multiplicity, used by DMC for branching */
	PSISQ,        /*!< square of the many-body wavefunction \f$|\Psi|^2\f$ */
	PSI,          /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
	AGE,          /*!< the age of the walker. set to zero when the walker is updated */
	SCALED,       /*!< scaling factor for the drift */
	NUMPROPERTIES /*!< the number of properties */
       };
  
  /**
   *\brief A container class to represent a walker.
   *
   * A walker stores the particle configurations 
   * and a property container.
   * The template (P)articleSet(A)ttribute is a generic container
   * of position types.
   */
  template<class T, class PA>
  struct Walker {

    typedef TinyVector<T,NUMPROPERTIES> PropertyContainer_t;

    ///scalar properties of a walker
    PropertyContainer_t  Properties;

    /**the configuration vector (3N-dimensional vector to store
       the positions of all the particles for a single walker)*/
    PA R;
    
    ///drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
    PA Drift;

    ///vector to store the constituents of the local energy
    std::vector<T> E;

    ///container for any data that are accessed by FIFO get/put functions
    PooledData<T> Data;

    ///create a walker for n-particles
    inline explicit Walker(int n) {  
      Properties(WEIGHT) = 1.0;
      Properties(MULTIPLICITY) = 1.0;
      Properties(PSI) = 1.0;
      resize(n);
    }
    
    ///constructor
    inline Walker(): Properties(0.0) {
      Properties(WEIGHT) = 1.0;
      Properties(MULTIPLICITY) = 1.0;
      Properties(PSI) = 1.0;
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

    inline void makeCopy(const Walker& a) {    
      resize(a.R.size());
      R = a.R;
      Drift = a.Drift;
      Properties = a.Properties;
      Data = a.Data;
      E = a.E;
    }

    /** reset the property of a walker
     *@param iw the id of this walker
     *@param psi the wavefunction value
     *@param ene the local energy
     *
     *Assign the values and reset the weight and multiplicity to one to start a run
     */
    inline void resetProperty(T psi, T ene) {
      Properties(WEIGHT) = 1.0;
      Properties(MULTIPLICITY) = 1.0;
      Properties(LOCALENERGY) = ene;
      Properties(PSISQ)=psi*psi;
      Properties(PSI)=psi;
    }

    ///resize for n-particles
    inline void resize(int n) {
      R.resize(n);
      Drift.resize(n);
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
