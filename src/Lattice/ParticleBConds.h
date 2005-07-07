//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_PARTICLE_BCONDS_H
#define OHMMS_PARTICLE_BCONDS_H

#include <limits>
#include <math.h>

/**generic ParticleBConds
 */
template<class T, unsigned D> struct ParticleBConds {};

/**class ParticleBConds<T,3> 
 *brief specialization for 3-dimensional supercells
 */
template<class T>
class ParticleBConds<T,3> {

public:

  // typedef for a pointer to boundary condition function
  typedef T (*ParticleBCond)(const T, const T, const T);

public:
  // constructor: initialize all BC's to periodic ones
  ParticleBConds() { }

  ///apply -> displacement
  inline TinyVector<T,3> apply(const TinyVector<T,3>& x) const {
    TinyVector<T,3> dr=x;
    if(dr[0]<-0.5) 
      dr[0]+=1.0;
    else if(dr[0]>=0.5) 
      dr[0]-=1.0;
    if(dr[1]<-0.5) 
      dr[1]+=1.0;
    else if(dr[1]>=0.5) 
      dr[1]-=1.0;
    if(dr[2]<-0.5) 
      dr[2]+=1.0;
    else if(dr[2]>=0.5) 
      dr[2]-=1.0;
    return dr;
 }

  ///wrap -> apply
  inline TinyVector<T,3> wrap(const TinyVector<T,3>& rin) const {
    register T x(rin[0]),y(rin[1]),z(rin[2]);
    const T epsilon = -numeric_limits<T>::epsilon();
    const T plus_one = 1.0;
    if(x<epsilon) x+=plus_one;
    else if(x>=plus_one) x-=plus_one;
    if(y<epsilon) y +=plus_one;
    else if(y>=plus_one) y-= plus_one;
    if(z<epsilon) z +=plus_one;
    else if(z >=plus_one) z -= plus_one;
    return TinyVector<T,3>(x,y,z);
  }

  ///applyBC
  inline void applyBC(const TinyVector<T,3>& rin, TinyVector<T,3>& rout) const {
    const T epsilon = -numeric_limits<T>::epsilon();
    const T plus_one = 1.0;
    rout=rin;
    if(rout[0]<epsilon)        rout[0] += plus_one;
    else if(rout[0]>=plus_one) rout[0] -= plus_one;
    if(rout[1]<epsilon)        rout[1] += plus_one;
    else if(rout[1]>=plus_one) rout[1] -= plus_one;
    if(rout[2]<epsilon)        rout[2] += plus_one;
    else if(rout[2]>=plus_one) rout[2] -= plus_one;
  }
};

#endif // OHMMS_PARTICLE_BCONDS_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
