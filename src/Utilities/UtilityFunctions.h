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
#ifndef OHMMS_UTILITYFUNCTIONS_H
#define OHMMS_UTILITYFUNCTIONS_H
/**@file UtilityFunctions.h
 *@brief A collection of utility functions. 
 */

/** Partition ntot over npart
 *\param ntot the total size
 *\param npart the number of partitions
 *\param adist an index array
 *
 *Simply divide ntot data among npart partitions so that the size of
 *each group cannot differ by more than 1.  
 *
 *The array adist contains the offset for the i-th group, 
 *i.e., adist[i+1]-adist[i] is the number of elements for the i-th group.
 */
template<class IV>
inline void FairDivide(int ntot, int npart, IV& adist) {
  adist.resize(npart+1);
  int bat=ntot/npart;
  int residue = ntot%npart;
  adist[0] = 0;
  for(int i=1; i<npart; i++) {
    if(i<residue) 
      adist[i] = adist[i-1] + bat+1;
    else 
      adist[i] = adist[i-1] + bat;
  }
  adist[npart]=ntot;
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
