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
  enum {Weight=0,
	Multiplicity,
	LocalEnergy, 
	LocalPotential, 
	PsiSq, 
	Sign, 
	Age,
	ScaleD, 
	NumProperties};
  
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

    typedef TinyVector<T,NumProperties> PropertyContainer_t;

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
      resize(n);
    }
    
    ///constructor
    inline Walker(): Properties(0.0) {
      Properties(Weight) = 1.0;
      Properties(Multiplicity) = 1.0;
      Properties(Sign) = 1.0;
    }

    ///copy constructor
    inline Walker(const Walker& a): R(a.R), Drift(a.Drift),
                                    Properties(a.Properties),
				    E(a.E),
				    Data(a.Data)
    { }
    
    ///assignment operator
    inline Walker& operator=(const Walker& a) {
      R = a.R;
      Drift = a.Drift;
      Properties = a.Properties;
      Data = a.Data;
      E = a.E;
      return *this;
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
