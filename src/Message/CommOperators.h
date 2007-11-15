//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
#ifndef OHMMS_COMMUNICATION_OPERATORS_H
#define OHMMS_COMMUNICATION_OPERATORS_H
#include "Message/Communicate.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"

template<class T> 
inline void gsum(T&, int) 
{ 
}

template<class T> 
inline void 
Communicate::allreduce(T& ) 
{
}

template<class T> 
inline void 
Communicate::reduce(T& ) 
{
}

/** dummy declaration to be specialized */
template<class T> 
inline void 
Communicate::reduce(T* restrict , T* restrict, int n) 
{
}

/** dummy declaration to be specialized */
template<class T> 
inline void 
Communicate::bcast(T* restrict ,int n) 
{
}

#ifdef HAVE_MPI
template<>
inline void gsum(int& g, int gid) {
  int gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

template<unsigned N>
inline void gsum(APPNAMESPACE::TinyVector<double,N>& g, int gid) 
{
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  APPNAMESPACE::TinyVector<double,N> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(std::vector<int>& g, int gid) {
  std::vector<int> gt(g.size(), 0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(double& g, int gid) {
  double gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

template<unsigned N>
inline void gsum(APPNAMESPACE::TinyVector<int,N>& g, int gid) 
{
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  APPNAMESPACE::TinyVector<int,N> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), N, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(std::vector<double>& g, int gid) {
  std::vector<double> gt(g.size(), 0.0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(APPNAMESPACE::Matrix<double>& g, int gid) 
{
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  vector<double> gt(g.size());
  std::copy(g.begin(),g.end(),gt.begin());
  MPI_Allreduce(g.data(), &gt[0], g.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::copy(gt.begin(),gt.end(),g.data());
}

////////////////////////////////
//template Communicate functions
////////////////////////////////
template<>
inline void 
Communicate::allreduce(int& g) 
{
  int gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_INT, MPI_SUM, myMPI);
}

template<>
inline void 
Communicate::allreduce(double& g) 
{
  double gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_DOUBLE, MPI_SUM, myMPI);
}

template<>
inline void 
Communicate::allreduce(APPNAMESPACE::TinyVector<double,OHMMS_DIM>& g) 
{
  APPNAMESPACE::TinyVector<double,OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(APPNAMESPACE::TinyVector<int,OHMMS_DIM>& g) 
{
  APPNAMESPACE::TinyVector<int,OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_INT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(std::vector<int>& g) 
{
  std::vector<int> gt(g.size(), 0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_INT,MPI_SUM,myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(std::vector<double>& g) 
{
  std::vector<double> gt(g.size(), 0.0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_DOUBLE,MPI_SUM,
      myMPI);
  g = gt;
}

template<>
inline void 
Communicate::reduce(std::vector<double>& g) 
{
  std::vector<double> gt(g.size(), 0);
  MPI_Reduce(&(g[0]),&(gt[0]),g.size(),MPI_DOUBLE,MPI_SUM,0,myMPI);
  if(master()) g = gt;
}

template<>
inline void 
Communicate::allreduce(APPNAMESPACE::Matrix<double>& g) 
{
  vector<double> gt(g.size());
  std::copy(g.begin(),g.end(),gt.begin());
  MPI_Allreduce(g.data(), &gt[0], g.size(), MPI_DOUBLE, MPI_SUM, 
      myMPI);
  std::copy(gt.begin(),gt.end(),g.data());
}

template<>
inline void 
Communicate::reduce(int* restrict g, int* restrict res, int n) 
{
  MPI_Reduce(g, res, n, MPI_INT, MPI_SUM, 0, myMPI);
}

template<>
inline void 
Communicate::reduce(double* restrict g, double* restrict res, int n) 
{
  MPI_Reduce(g, res, n, MPI_DOUBLE, MPI_SUM, 0, myMPI);
}

template<>
inline void
Communicate::bcast(double* restrict x, int n) 
{
  MPI_Bcast(x,n,MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(int* restrict x, int n) 
{
  MPI_Bcast(x,n,MPI_INT,0,myMPI);
}

#endif

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
