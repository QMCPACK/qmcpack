//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef OHMMS_COMMUNICATION_OPERATORS_SINGLE_H
#define OHMMS_COMMUNICATION_OPERATORS_SINGLE_H

///dummy declarations to be specialized
template<typename T> inline void gsum(T&, int) { }

template<typename T> inline void Communicate::allreduce(T& ) { }

template<typename T> inline void Communicate::reduce(T& ) { }

template<typename T>
inline void Communicate::reduce(T* restrict , T* restrict, int n) { }

template<typename T>
inline void Communicate::reduce_in_place(T* restrict, int n) { }

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

template<typename T> inline void Communicate::allgather(T& sb, T& rb, int count)
{
  for(size_t i=0; i<count; i++) rb[i]=sb[i];
}

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

template<typename T> 
void allreduce(T&,Communicate::mpi_comm_type comm)
{ }

template<typename T> 
void bcast(T&,Communicate::mpi_comm_type)
{ }

template<typename T> 
void bcast(T* restrict, int n,Communicate::mpi_comm_type comm)
{ }

template<typename T> 
void bcast(T* restrict, int n, int orig, Communicate::mpi_comm_type comm)
{ }

template<typename T> 
void send(T* restrict, int n, int dest, int tag, Communicate::mpi_comm_type comm)
{ }

#ifdef HAVE_MPI
template<typename T> 
void recv(T* restrict, int n, int dest, int tag, Communicate::mpi_comm_type comm, MPI_Status*)
{ }
#endif

template<typename T, typename IT> 
void gatherv(T* sb, T* rb, int n, IT& counts, IT& displ, int dest)
{ }

#ifdef HAVE_MPI
template<typename T, typename IT> 
void gatherv(T* sb, T* rb, int n,IT& counts, IT& displ, int dest, MPI_Comm comm)
{ }
#endif

template<typename T> 
void allgather(T& sb, T& rb, int count, Communicate::mpi_comm_type comm)
{ }

template<typename T> 
void allgather(T* sb, T* rb, int count)
{ }

#ifdef HAVE_MPI
template<typename T, typename IT> 
void scatterv(T* sb, T* rb, int n, IT& counts, IT& displ, int source, MPI_Comm)
{ }
#endif

template<typename T> 
void gsum(T&)
{ }

#ifdef HAVE_MPI
template<typename T> 
void gsum(T&,Communicate::mpi_comm_type comm)
{ }
#endif

template<typename T> 
void gmax(T&,Communicate::mpi_comm_type comm)
{ }

#endif

