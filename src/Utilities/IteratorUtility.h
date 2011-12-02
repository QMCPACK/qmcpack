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
      while(first != last) { if(*first) delete *first; ++first;}
    }


  template<typename IT1, typename IT2>
    inline void accumulate_elements(IT1 first, IT1 last, IT2 res)
    {
      while(first != last) *res++ += *first++;
    }

//  template<typename IT1, typename IT2, typename INT>
//    inline void accumulate_elements(IT1 first, IT2 res, INT n)
//    {
//      for(;n>0; n--) *res++ += *first++;
//    }

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
