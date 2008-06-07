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
#ifndef OHMMS_COMMUNICATION_OPERATORS_MPI_H
#define OHMMS_COMMUNICATION_OPERATORS_MPI_H
#include "Utilities/PooledData.h"
///dummy declarations to be specialized
template<typename T> inline void gsum(T&, int) 
{ 
  APP_ABORT("Need specialization for gsum(T&, int)");
}

template<typename T> inline void Communicate::allreduce(T& ) 
{ 
  APP_ABORT("Need specialization for allreduce(T&)");
}

template<typename T> inline void Communicate::reduce(T& ) 
{ 
  APP_ABORT("Need specialization for reduce(T&)");
}

template<typename T> 
inline void 
Communicate::reduce(T* restrict , T* restrict, int n) 
{ 
  APP_ABORT("Need specialization for reduce(T* restrict , T* restrict, int n)");
}

template<typename T> inline void 
Communicate::bcast(T& ) 
{  
  APP_ABORT("Need specialization for bcast(T&)");
}

template<typename T> inline void 
Communicate::bcast(T* restrict ,int n) 
{ 
  APP_ABORT("Need specialization for bcast(T* restrict ,int n)");
}

template<typename T> inline void
Communicate::send(int dest, int tag, T&) 
{ 
  APP_ABORT("Need specialization for send(int, int, T& )");
}

template<typename T> inline void Communicate::gather(T& sb, T& rb, int dest)
{ 
  APP_ABORT("Need specialization for gather(T&, T&, int)");
}

template<typename T> inline void Communicate::scatter(T& sb, T& rb, int dest)
{ 
  APP_ABORT("Need specialization for scatter(T&, T&, int)");
}

template<typename T> inline Communicate::request
Communicate::irecv(int source, int tag, T& ) 
{ 
  APP_ABORT("Need specialization for irecv(int source, int tag, T& )");
  return 1;
}

template<typename T> inline Communicate::request
Communicate::isend(int dest, int tag, T&)
{
  APP_ABORT("Need specialization for isend(int source, int tag, T& )");
  return 1;
}

template<typename T> inline Communicate::request
Communicate::irecv(int source, int tag, T* , int n) 
{ 
  APP_ABORT("Need specialization for irecv(int source, int tag, T*, int )");
  return 1;
}

template<typename T> inline Communicate::request
Communicate::isend(int dest, int tag, T*, int n)
{
  APP_ABORT("Need specialization for isend(int source, int tag, T*, int )");
  return 1;
}

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
inline void gsum(std::vector<double>& g, int gid) 
{
  std::vector<double> gt(g.size(), 0.0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(APPNAMESPACE::Matrix<double>& g, int gid) 
{
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::vector<double> gt(g.size());
  std::copy(g.begin(),g.end(),gt.begin());
  MPI_Allreduce(g.data(), &gt[0], g.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  std::copy(gt.begin(),gt.end(),g.data());
}

template<>
inline void 
Communicate::allreduce(int& g) 
{
  if(d_ncontexts==1) return;
  int gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_INT, MPI_SUM, myMPI);
}

template<>
inline void 
Communicate::allreduce(long& g) 
{
  if(d_ncontexts==1) return;
  long gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_LONG, MPI_SUM, myMPI);
}


template<>
inline void 
Communicate::allreduce(double& g) 
{
  if(d_ncontexts==1) return;
  double gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_DOUBLE, MPI_SUM, myMPI);
}

template<>
inline void 
Communicate::allreduce(APPNAMESPACE::TinyVector<double,OHMMS_DIM>& g) 
{
  if(d_ncontexts==1) return;
  APPNAMESPACE::TinyVector<double,OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_DOUBLE, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(APPNAMESPACE::TinyVector<int,OHMMS_DIM>& g) 
{
  if(d_ncontexts==1) return;
  APPNAMESPACE::TinyVector<int,OHMMS_DIM> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), OHMMS_DIM, MPI_INT, MPI_SUM, myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(std::vector<int>& g) 
{
  if(d_ncontexts==1) return;
  std::vector<int> gt(g.size(), 0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_INT,MPI_SUM,myMPI);
  g = gt;
}

template<>
inline void 
Communicate::allreduce(std::vector<long>& g) 
{
  if(d_ncontexts==1) return;
  std::vector<long> gt(g.size(), 0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_LONG,MPI_SUM,myMPI);
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
  std::vector<double> gt(g.size());
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
Communicate::bcast(bool &g)
{
  int val = g ? 1 : 0;
  MPI_Bcast(&val,1,MPI_INT,0,myMPI);
  g = val != 0;
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<double,2>& g) 
{
  MPI_Bcast(g.begin(),2,MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<int,2>& g) 
{
  MPI_Bcast(g.begin(),2,MPI_INT,0,myMPI);
}

template<>
inline void 
Communicate::bcast(APPNAMESPACE::TinyVector<int,3>& g) 
{
  MPI_Bcast(g.begin(),3,MPI_INT,0,myMPI);
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
Communicate::bcast(APPNAMESPACE::Tensor<double,3>& g) 
{
  MPI_Bcast(&(g[0]),9,MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(APPNAMESPACE::Vector<double>& g)
{
  MPI_Bcast(&(g[0]),g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(APPNAMESPACE::Vector<int>& g)
{
  MPI_Bcast(&(g[0]),g.size(),MPI_INT,0,myMPI);
}

template<>
inline void
Communicate::bcast(APPNAMESPACE::Vector<APPNAMESPACE::TinyVector<double,2> >& g)
{
  MPI_Bcast(&(g[0]),2*g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(APPNAMESPACE::Vector<APPNAMESPACE::TinyVector<double,3> >& g)
{
  MPI_Bcast(&(g[0]),3*g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(Array<double,3> &g)
{
  MPI_Bcast(g.data(), g.size(), MPI_DOUBLE, 0, myMPI);
}

template<>
inline void
Communicate::bcast(Array<complex<double>,3> &g)
{
  MPI_Bcast(g.data(), 2*g.size(), MPI_DOUBLE, 0, myMPI);
}


template<>
inline void 
Communicate::bcast(std::vector<double>& g) 
{
  MPI_Bcast(&(g[0]),g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void 
Communicate::bcast(PooledData<double>& g) 
{
  MPI_Bcast(&(g[0]),g.size(),MPI_DOUBLE,0,myMPI);
}

template<>
inline void
Communicate::bcast(std::vector<APPNAMESPACE::TinyVector<double,2> > &g)
{
  MPI_Bcast(&(g[0][0]), 2*g.size(), MPI_DOUBLE, 0, myMPI);
}

template<>
inline void
Communicate::bcast(std::vector<APPNAMESPACE::TinyVector<double,3> > &g)
{
  MPI_Bcast(&(g[0][0]), 3*g.size(), MPI_DOUBLE, 0, myMPI);
}


template<>
inline void 
Communicate::bcast(std::vector<int>& g) 
{
  MPI_Bcast(&(g[0]),g.size(),MPI_INT,0,myMPI);
}

template<>
inline void 
Communicate::bcast(std::vector<bool>& g) 
{
  std::vector<int> intVec(g.size());
  for (int i=0; i<g.size(); i++)
    intVec[i] = g[i] ? 1 : 0;
  MPI_Bcast(&(intVec[0]),g.size(),MPI_INT,0,myMPI);
  for (int i=0; i<g.size(); i++)
    g[i] = intVec[i] != 0;
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

template<> inline void
Communicate::send(int dest, int tag, std::vector<double>& g)
{
  MPI_Send(&(g[0]),g.size(),MPI_DOUBLE,dest,tag, myMPI);
}

template<> inline Communicate::request
Communicate::isend(int dest, int tag, std::vector<double>& g)
{
  request r;
  MPI_Isend(&(g[0]),g.size(),MPI_DOUBLE,dest,tag, myMPI,&r);
  return r;
}

template<> inline Communicate::request
Communicate::irecv(int source, int tag, std::vector<double>& g)
{
  request r;
  MPI_Irecv(&(g[0]),g.size(),MPI_DOUBLE,source,tag, myMPI,&r);
  return r;
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 2635 $   $Date: 2008-04-25 16:46:48 -0500 (Fri, 25 Apr 2008) $
 * $Id: CommOperators.h 2635 2008-04-25 21:46:48Z kesler $ 
 ***************************************************************************/
