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
#ifndef OHMMS_BITSET_IBMGCC_H
#define OHMMS_BITSET_IBMGCC_H

template<unsigned N>
struct bitset {

  bool data[N];

  bitset() { reset();}
  inline void reset() { for(int i=0; i<N; i++) data[i] = false;}
  inline void set(int i) {data[i] = true;}
  inline void flip(int i) { 
    data[i] =  (data[i])? false: true;
  }

  inline bool operator[](int i) const { return data[i];}
  inline bool& operator[](int i) { return data[i];}

  inline bool any() const {
    int i=0;
    while(i<N) { if(data[i++]) return true;}
    return false;
  }
};
#endif // OHMMS_BITSET_IBMGCC_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
