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

template<class T> inline void gsum(T&, int) { }

#ifdef HAVE_MPI
template<>
inline void gsum(int& g, int gid) {
  int gt = g;
  MPI_Allreduce(&(gt), &(g), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

template<unsigned N>
inline void gsum(TinyVector<double,N>& g, int gid) {
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  TinyVector<double,N> gt(g);
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
inline void gsum(TinyVector<int,N>& g, int gid) {
  //TinyVector<double,N> gt = g;
  //MPI_Allreduce(gt.begin(), g.begin(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  TinyVector<int,N> gt(g);
  MPI_Allreduce(g.begin(), gt.begin(), N, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  g = gt;
}

template<>
inline void gsum(std::vector<double>& g, int gid) {
  std::vector<double> gt(g.size(), 0.0);
  MPI_Allreduce(&(g[0]),&(gt[0]),g.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  g = gt;
}


#endif

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
