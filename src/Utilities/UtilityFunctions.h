//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include<cstdlib>
#include<tuple>

#ifndef OHMMS_UTILITYFUNCTIONS_H
#define OHMMS_UTILITYFUNCTIONS_H
/**@file UtilityFunctions.h
 *@brief A collection of utility functions.
 */

/** Partition ntot over npart
 *\param me index of boundaries being returned
 *\param ntot the total size
 *\param npart the number of partitions
 *
 *Simply divide ntot data among npart partitions so that the size of
 *each group cannot differ by more than 1.
 *Return the boundaries of the partition associated with partition 'me'
 *
 */
template<typename IType>
inline std::tuple<IType,IType> FairDivideBoundary(IType me, IType ntot, IType npart)
{
  IType bat=ntot/npart;
  IType residue = ntot%npart;
  if(me < residue)
    return std::make_tuple( me*(bat+1), (me+1)*(bat+1) );
  else 
    return std::make_tuple( me*bat + residue, (me+1)*bat + residue );    
}

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
inline void FairDivide(int ntot, int npart, IV& adist)
{
  if(adist.size() != npart+1)
    adist.resize(npart+1);
  int bat=ntot/npart;
  int residue = ntot%npart;
  adist[0] = 0;
  for(int i=1; i<npart; i++)
  {
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
inline void FairDivideLow(int ntot, int npart, IV& adist)
{
  if(adist.size() != npart+1)
    adist.resize(npart+1);
  int bat=ntot/npart;
  int residue = npart-ntot%npart;
  adist[0] = 0;
  for(int i=0; i<npart; i++)
  {
    if(i<residue)
      adist[i+1] = adist[i] + bat;
    else
      adist[i+1] = adist[i] + bat+1;
  }
}

/** partition ntot elements among npart
 * @param me rank [0,ntot)
 * @param ntot total number of elements
 * @param npart number of partitions
 * @param adist distribution offset
 * @return partition id to which me belongs
 *
 * mypart satisfies adist[mypart] <= me < adist[mypart+1]
 */
template<class IV>
inline int FairDivideHigh(int me, int ntot, int npart, IV& adist)
{
  if(adist.size() != npart+1)
    adist.resize(npart+1);
  int bat=ntot/npart;
  int residue = ntot%npart;
  int mypart=0;
  adist[0] = 0;
  for(int i=1; i<npart; i++)
  {
    if(i<residue)
      adist[i] = adist[i-1] + bat+1;
    else
      adist[i] = adist[i-1] + bat;
    if(me>= adist[i] && me<adist[i+1])
      mypart=i;
  }
  adist[npart]=ntot;
  return mypart;
}

/** partition ntot elements among npart
 * @param me rank [0,ntot)
 * @param ntot total number of elements
 * @param npart number of partitions
 * @param adist distribution offset
 * @return partition id to which me belongs
 *
 * mypart satisfies adist[mypart] <= me < adist[mypart+1]
 */
template<class IV>
inline int FairDivideLow(int me, int ntot, int npart, IV& adist)
{
  if(adist.size() != npart+1)
    adist.resize(npart+1);
  int bat=ntot/npart;
  int residue = npart-ntot%npart;
  int mypart=0;
  adist[0] = 0;
  for(int i=0; i<npart; i++)
  {
    if(i<residue)
      adist[i+1] = adist[i] + bat;
    else
      adist[i+1] = adist[i] + bat+1;
    if(me>= adist[i] && me<adist[i+1])
      mypart=i;
  }
  return mypart;
}

#endif
