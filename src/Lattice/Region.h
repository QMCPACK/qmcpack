//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
#ifndef OHMMS_REGINON_H
#define OHMMS_REGION_H

/* \class Region defines a spatial region bound by [Ri,Rf)
   \brief Defined in unit vectors 0 <= Ri, Rf < 1
*/
template<class T, unsigned D>
struct Region {

  typedef T Scalar_t;
  enum {DIM = D};
  T Ri[D], Rf[D];
  
  Region(){ }

  Region(const T* r0, const T* dr) {
    set(r0,dr);
  }

  Region(const Region<T,D>& rg) {
    for(int i=0; i<D; i++) Ri[i] = rg.Ri[i];
    for(int i=0; i<D; i++) Rf[i] = rg.Rf[i];
  }

  Region<T,D>& operator=(const Region<T,D>& rg) {
    for(int i=0; i<D; i++) Ri[i] = rg.Ri[i];
    for(int i=0; i<D; i++) Rf[i] = rg.Rf[i];
    return *this;
  }

  ~Region(){ }

  inline void set(const T* r0, const T*  dr) {
    for(int i=0; i<D; i++) Ri[i] = r0[i];
    for(int i=0; i<D; i++) Rf[i] = r0[i]+dr[i];
  }

  template<class Pos_t>
  inline bool inside(const Pos_t& r) const {
    for(int i=0; i<DIM; i++) { 
      if(r[i] < Ri[i] || r[i] >= Rf[i]) return false;
    }
    return true;
  }
};
#endif
  


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
