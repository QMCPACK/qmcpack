//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/***************************************************************************
 *
 * The POOMA Framework
 * 
 * This program was prepared by the Regents of the University of
 * California at Los Alamos National Laboratory (the University) under
 * Contract No.  W-7405-ENG-36 with the U.S. Department of Energy (DOE).
 * The University has certain rights in the program pursuant to the
 * contract and the program should not be copied or distributed outside
 * your organization.  All rights in the program are reserved by the DOE
 * and the University.  Neither the U.S.  Government nor the University
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMA for more details
 *
 ***************************************************************************/

#ifndef PARTICLE_BCONDS_H
#define PARTICLE_BCONDS_H

#include <math.h>

/*! \note Modification of ParticleBConds.h of Pooma r1 to work with
 *  changed TinyVector classes.
 */
//////////////////////////////////////////////////////////////////////
// particle boundary condition functions ...
// ParticleNoBCond = no action is taken
// ParticlePeriodicBCond = shift the vector withtin [-len/2, len/2)
// ParticlePeriodicWrapBCond = wrap around at endpoints of the interval
// ParticleReflectiveBCond = values bounce back from endpoints
//////////////////////////////////////////////////////////////////////

// null BC; value is not changed


//! \param t A position subject to the boundary condition.
//! \param minval mininum range (ignored)
//! \param maxval maximum range (ignored)
//! \return t
template<class T>
inline T 
ParticleNoBCond(const T t, const T /* minval */, const T /* maxval */) {
  return t;
}

//! \param t A position subject to the boundary condition.
//! \param minval mininum range (ignored)
//! \param maxval maximum range (ignored)
//! \return t if |t| < (maxval-minval)/2
//! \note This funtion returns the shortest distance of 
//!  t under periodic boundary condition. minval is ignored for periodic BC.
//!  [-maxval/2, maxval/2) is the actual range
template<class T>
inline T 
ParticlePeriodicBCond(const T t, const T minval, const T maxval) {
  //T t0=fmod(t,maxval);
  //return t0 - static_cast<T>(static_cast<int>(t0*2.0*maxval));
  if(t < -0.5) return t+1.0;
  else if(t>= 0.5) return t-1.0;
  else return t;
}

//! \param t A position subject to the boundary condition.
//! \param minval mininum range (ignored)
//! \param maxval maximum range (ignored)
//! \return t
//! \note Poomar1::ParticlePeriodicBCond is ParticlePeriodicWrapBCond
template<class T>
inline T
ParticlePeriodicWrapBCond(const T t, const T minval, const T maxval) {
  if (t < minval)
    return (maxval - (minval - t));
  else if (t >= maxval)
    return (minval + (t - maxval));
  else
    return t;
}


//! \note reflective BC; values bounce back from endpoints
//! Implemented by Poomar1. 
template<class T>
inline T 
ParticleReflectiveBCond(const T t, const T minval, const T maxval) {
  if (t < minval)
    return (minval + (minval - t));
  else if (t >= maxval)
    return (maxval - (t - maxval));
  else
    return t;
}

//! \note sink BC; particles stick to the selected face
//! Implemented by Poomar1
template<class T>
T ParticleSinkBCond(const T t, const T minval, const T maxval) {
  if (t < minval)
    return minval;
  else if (t >= maxval)
    return maxval;
  else
    return t;
}

////////////////////////////////////////////////////////////////
//
// Fuction classes to apply BC to a TinyVector
//
////////////////////////////////////////////////////////////////
// forward declaration
template<class T, unsigned D> class ParticleBConds;

/*! \class template<class T> struct applyBConds
 *  \brief Empty class. Specialization of applyBConds with specific
 *  ParticleBConds is intended.
 */
template<class T>
struct applyBConds {
};

/*!\class template<class T, unsigned D> struct applyBConds<ParticleBConds<T,D> >
 * \brief Specialization for ParticleBConds<T,D> defined below.
 */
template<class T, unsigned D>
struct applyBConds<ParticleBConds<T,D> > {

  //!< Apply boundary conditions for each direction
  static inline TinyVector<T,D> 
  apply(const ParticleBConds<T,D>& bc, const TinyVector<T,D>& a){
    TinyVector<T,D> res;
    for(int i=0; i<D; i++) res[i] = bc.apply(a[i],i,1.0);
    return res;
  }
};

/*!\class template<class T> struct applyBConds<ParticleBConds<T,1> >
 * \brief Specialization for ParticleBConds<T,D> defined below.
 */
template<class T>
struct applyBConds<ParticleBConds<T,1> > {

  static inline TinyVector<T,1> 
  apply(const ParticleBConds<T,1>& bc, const TinyVector<T,1>& a){
    return 
    TinyVector<T,1>(bc.apply(a[0],0,1));
  }
};

/*!\class template<class T> struct applyBConds<ParticleBConds<T,2> >
 * \brief Specialization for ParticleBConds<T,D> defined below.
 */
template<class T>
struct applyBConds<ParticleBConds<T,2> > {

  static inline TinyVector<T,2> 
  apply(const ParticleBConds<T,2>& bc, const TinyVector<T,2>& a){
    return 
    TinyVector<T,2>(bc.apply(a[0],0,1),bc.apply(a[1],1,1));
  }
};

/*!\class template<class T> struct applyBConds<ParticleBConds<T,3> >
 * \brief Specialization for ParticleBConds<T,D> defined below.
 */
template<class T>
struct applyBConds<ParticleBConds<T,3> > {

  static inline TinyVector<T,3> 
    apply(const ParticleBConds<T,3>& bc, const TinyVector<T,3>& a){
    return 
      TinyVector<T,3>(bc.apply(a[0],0,1),bc.apply(a[1],1,1),bc.apply(a[2],2,1));
  }

  static inline TinyVector<T,3>
    wrap(const ParticleBConds<T,3>& bc, const TinyVector<T,3>& a){
    return 
      TinyVector<T,3>(bc.wrap(a[0],0,1),bc.wrap(a[1],1,1),
		      bc.wrap(a[2],2,1));
  }
};

/*!\class ParticleBConds
 * \brief ParticleBConds is a container for a set of particle boundary condition
 * functions.  Boundary conditions for particles are not objects, but just
 * functions which map a position X -> X', given the minimum and maximum
 * values of the spatial domain.
 ***************************************************************************/
//////////////////////////////////////////////////////////////////////
// general container for a set of particle boundary conditions
template<class T, unsigned Dim>
class ParticleBConds {

public:
  // typedef for a pointer to boundary condition function
  typedef T (*ParticleBCond)(const T, const T, const T);

public:
  // constructor: initialize all BC's to periodic ones
  ParticleBConds() {
    for (int d=(Dim - 1); d >= 0; --d) {
      BCList[d] = ParticlePeriodicBCond; // for distance
      BCTransList[d] = ParticlePeriodicWrapBCond; // for translation, wrapping
    }
  }

  // operator= to copy values from another container
  ParticleBConds<T,Dim>& operator=(const ParticleBConds<T,Dim>& pbc) {
    for (int d=(Dim - 1); d >= 0; --d) {
      BCList[d] = pbc.BCList[d];
      BCTransList[d] = ParticlePeriodicWrapBCond;
    }
    return *this;
  }

  // operator[] to get value of Nth boundary condition
  inline ParticleBCond& operator[](unsigned d) { return BCList[d]; }

  // return a wrpper function for d-direction
  inline ParticleBCond& wrapper(unsigned d) { return BCTransList[d]; }

  // for the given value in the given dimension over the given NDRegion,
  // apply the proper BC and return back the new value
  inline T apply(const T t, const unsigned d, const T len) const {
    return BCList[d](t,0,len);
  }

  // for the given value in the given dimension over the given NDRegion,
  // apply the proper BC and return back the new value
  inline T wrap(const T t, const unsigned d, const T len) const {
    return BCTransList[d](t,0,len);
  }

  // to work with difffent types of data
  template<class T1>
  inline T1 apply(const T1 t, const unsigned d, const T1 len) const {
    return static_cast<T1>(
      BCList[d](static_cast<T>(t), 0, static_cast<T>(len)));
  }

  inline TinyVector<T,Dim> apply(const TinyVector<T,Dim>& x) const {
    return applyBConds<ParticleBConds<T,Dim> >::apply(*this,x);
  }

  inline TinyVector<T,Dim> wrap(const TinyVector<T,Dim>& x) const {

    return applyBConds<ParticleBConds<T,Dim> >::wrap(*this,x);    
  }

private:
  // array storing the function pointers
  // Note that the size is Dim not 2*Dim
  ParticleBCond BCList[Dim];
  ParticleBCond BCTransList[Dim];
};

#endif // PARTICLE_BCONDS_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
