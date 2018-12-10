//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
  std::vector<Communicate::request> r(first,last);
  std::vector<Communicate::status> st(r.size());
  MPI_Waitall(r.size(),&(r[0]),&(st[0]));
}

template<typename CT> inline void wait_all(CT& requests)
{
  std::vector<Communicate::status> st(requests.size());
  MPI_Waitall(requests.size(),&(requests[0]),&(st[0]));
}


inline void wait_all(int n, Communicate::request* pending)
{
  std::vector<Communicate::status> st(n);
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

