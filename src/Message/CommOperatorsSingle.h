//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
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
#ifndef OHMMS_COMMUNICATION_OPERATORS_SINGLE_H
#define OHMMS_COMMUNICATION_OPERATORS_SINGLE_H

///dummy declarations to be specialized
template<typename T> inline void gsum(T&, int) { }

template<typename T> inline void Communicate::allreduce(T& ) { }

template<typename T> inline void Communicate::reduce(T& ) { }

template<typename T>
inline void Communicate::reduce(T* restrict , T* restrict, int n) { }

template<typename T> inline void Communicate::bcast(T& ) {  }

template<typename T> inline void Communicate::bcast(T* restrict ,int n) { }

template<typename T> inline Communicate::request
Communicate::irecv(int source, int tag, T& )
{
  return 1;
}

template<typename T> inline void
Communicate::send(int dest, int tag, T&) { }

template<typename T> inline void Communicate::gather(T& sb, T& rb, int dest) { }

template<typename T> inline void Communicate::allgather(T& sb, T& rb, int count) { }

template<typename T> inline void Communicate::scatter(T& sb, T& rb, int dest) { }

template<typename T> inline Communicate::request
Communicate::isend(int dest, int tag, T&)
{
  return 1;
}

template<typename T> inline Communicate::request
Communicate::irecv(int source, int tag, T* , int n)
{
  return 1;
}

template<typename T> inline Communicate::request
Communicate::isend(int dest, int tag, T*, int n)
{
  return 1;
}

template<typename T, typename IT>
inline void Communicate::gatherv(T& sb, T& rb, IT&, IT&, int dest)
{
}

template<typename T, typename IT>
inline void Communicate::scatterv(T& sb, T& rb, IT&, IT&, int source)
{
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 2635 $   $Date: 2008-04-25 16:46:48 -0500 (Fri, 25 Apr 2008) $
 * $Id: CommOperators.h 2635 2008-04-25 21:46:48Z kesler $
 ***************************************************************************/
