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

/** partition ntot elements among npart
 * @param ntot total number of elements
 * @param npart number of partitions
 * @param adist distribution offset 
 *
 * adist[ip-1]-adist[ip] is the number of elements of ip partition
 * This method makes the zero-th node equal to or less than 1.
 */
template<class IV>
inline void FairDivideLow(int ntot, int npart, IV& adist) {
  int bat=ntot/npart;
  int residue = npart-ntot%npart;
  adist[0] = 0;
  for(int i=0; i<npart; i++) {
    if(i<residue)
      adist[i+1] = adist[i] + bat;
    else
      adist[i+1] = adist[i] + bat+1;
  }
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
