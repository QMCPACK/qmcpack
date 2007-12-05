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

///dummy declarations to be specialized
template<typename T> inline void gsum(T&, int) { }

template<typename T> inline void Communicate::allreduce(T& ) { }

template<typename T> inline void Communicate::reduce(T& ) { }

template<typename T> 
inline void 
Communicate::reduce(T* restrict , T* restrict, int n) { }

template<typename T> inline void 
Communicate::bcast(T& ) { }

template<typename T> inline void 
Communicate::bcast(T* restrict ,int n) { }

template<typename T> inline Communicate::request
Communicate::irecv(int source, int tag, T& ) 
{ 
  return 1;
}

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
  if(!d_mycontext) g = gt;
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
Communicate::bcast(int& g) 
{
  MPI_Bcast(&g,1,MPI_INT,0,myMPI);
}

template<>
inline void 
Communicate::bcast(double& g) 
{
  MPI_Bcast(&g,1,MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<double,2>& g) 
{
  MPI_Bcast(g.begin(),2,MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<double,3>& g) 
{
  MPI_Bcast(g.begin(),3,MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<double,4>& g) 
{
  MPI_Bcast(g.begin(),4,MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(std::vector<double>& g) 
{
  MPI_Bcast(&(g[0]),g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(std::vector<int>& g) 
{
  MPI_Bcast(&(g[0]),g.size(),MPI_INT,0,myMPI);
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

template<> inline Communicate::request
Communicate::isend(int dest, int tag, vector<double>& g)
{
  request r;
  MPI_Isend(&(g[0]),g.size(),MPI_DOUBLE,dest,tag, myMPI,&r);
  return r;
}

template<> inline Communicate::request
Communicate::irecv(int source, int tag, vector<double>& g)
{
  request r;
  MPI_Irecv(&(g[0]),g.size(),MPI_DOUBLE,source,tag, myMPI,&r);
  return r;
}
#endif
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
