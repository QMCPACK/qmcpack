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

namespace ohmmsqmc {

  /** an enum denoting index of physical properties */
  enum {WEIGHT=0,     /*!< weight */
	LOCALENERGY,  /*!< local energy, the sum of all the components */
	LOCALPOTENTIAL, /*!< local potential energy = local energy - kinetic energy */
	MULTIPLICITY, /*!< multiplicity, used by DMC for branching */
	PSISQ,        /*!< square of the many-body wavefunction \f$|\Psi|^2\f$ */
	PSI,          /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
	AGE,          /*!< the age of the walker. set to zero when the walker is updated */
	SCALED,       /*!< scaling factor for the drift */
	NUMPROPERTIES, /*!< the number of properties */
	CAPACITY=15
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

    typedef TinyVector<T,CAPACITY> PropertyContainer_t;
    //typedef std::vector<T> PropertyContainer_t;

    ///scalar properties of a walker
    PropertyContainer_t  Properties;

    /**the configuration vector (3N-dimensional vector to store
       the positions of all the particles for a single walker)*/
    PA R;
    
    ///drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
    PA Drift;

    ///create a walker for n-particles
    inline explicit Walker(int n) : Properties(0.0) {  
      Properties[WEIGHT] = 1.0;
      Properties[MULTIPLICITY] = 1.0;
      Properties[PSI] = 1.0;
      resize(n);
    }
    
    ///constructor
    inline Walker() : Properties(0.0) {
      Properties[WEIGHT] = 1.0;
      Properties[MULTIPLICITY] = 1.0;
      Properties[PSI] = 1.0;
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

    inline void makeCopy(const Walker& a) {    
      resize(a.R.size());
      R = a.R;
      Drift = a.Drift;
      Properties = a.Properties;
    }

    //return the address of the values of Hamiltonian terms
    inline T* restrict getEnergyBase() {
      return Properties.begin()+NUMPROPERTIES;
    }

    //return the address of the values of Hamiltonian terms
    inline const T* restrict getEnergyBase() const {
      return Properties.begin()+NUMPROPERTIES;
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


    /** reset the walker weight, multiplicity and age */
    inline void reset() {
      Properties(WEIGHT)=1.0;
      Properties(MULTIPLICITY)=1.0;
      Properties(AGE)=0.0;
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
