//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file CommUtilities.h
 * @brief define convenience functions for mpi operations.
 */
#ifndef OHMMS_COMMUNICATION_UTILITY_FUNCTIONS_H
#define OHMMS_COMMUNICATION_UTILITY_FUNCTIONS_H
#include "Message/Communicate.h"
#if defined(HAVE_MPI)
namespace qmcplusplus
{
template<typename IT> inline void wait_all(IT first, IT last)
{
  vector<Communicate::request> r(first,last);
  vector<Communicate::status> st(r.size());
  MPI_Waitall(r.size(),&(r[0]),&(st[0]));
}

template<typename CT> inline void wait_all(CT& requests)
{
  vector<Communicate::status> st(requests.size());
  MPI_Waitall(requests.size(),&(requests[0]),&(st[0]));
}


inline void wait_all(int n, Communicate::request* pending)
{
  vector<Communicate::status> st(n);
  MPI_Waitall(n,pending,&(st[0]));
}

template<typename CT>
inline void cancel(CT& r)
{
  for(int i=0; i<r.size(); i++)
    MPI_Cancel(&r[i]);
}

template<typename IT>
inline void cancel(IT first, IT last)
{
  while(first != last)
  {
    MPI_Cancel(&(*first));
    ++first;
  }
}

template<typename T>
inline void bcast(T& a, Communicate* comm)
{
  comm->bcast(a);
}

}
#else
namespace qmcplusplus
{
template<typename CT>
inline void cancel(CT& r) { }

template<typename T>
inline void bcast(T& a, Communicate* comm) { }
}
#endif
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2298 $   $Date: 2007-11-15 15:09:53 -0600 (Thu, 15 Nov 2007) $
 * $Id: CommOperators.h 2298 2007-11-15 21:09:53Z jnkim $
 ***************************************************************************/
