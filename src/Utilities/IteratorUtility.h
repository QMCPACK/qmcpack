//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef OHMMS_ITERATOR_UTILITIES_H
#define OHMMS_ITERATOR_UTILITIES_H

namespace qmcplusplus {
  /** delete the pointers in [first,last)
  */
  template<class IT>
    inline void delete_iter(IT first, IT last) {
      while(first != last) { delete *first; ++first;}
    }


  template<class T, unsigned D>
    inline T* get_first_address(ParticleAttrib<TinyVector<T,D> >& a) {
      return &(a[0][0]);
    }

  template<class T, unsigned D>
    inline T* get_last_address(ParticleAttrib<TinyVector<T,D> >& a) {
      return &(a[0][0])+D*a.size();
    }

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
