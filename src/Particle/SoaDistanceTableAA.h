//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_DTDIMPL_AA_H
#define QMCPLUSPLUS_DTDIMPL_AA_H
#include "simd/algorithm.hpp"

namespace qmcplusplus
{

/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAA: public DTD_BConds<T,D,SC>, public DistanceTableData
{

  int Ntargets;
  int Ntargets_padded;

  SoaDistanceTableAA(ParticleSet& target)
    : DTD_BConds<T,D,SC>(target.Lattice), DistanceTableData(target,target)
  {
    resize(target.getTotalNum());
  }

#if (__cplusplus >= 201103L)
  SoaDistanceTableAA()=delete;
  SoaDistanceTableAA(const SoaDistanceTableAA&)=delete;
#endif
  ~SoaDistanceTableAA() {}

  size_t compute_size(int N)
  {
    const size_t N_padded = getAlignedSize<T>(N);
    const size_t Alignment = getAlignment<T>();
    return (N_padded*(2*N-N_padded+1)+(Alignment-1)*N_padded)/2;
  }

  void resize(int n)
  {
    N[SourceIndex]=N[VisitorIndex]=Ntargets=n;
    Ntargets_padded=getAlignedSize<T>(n);
    Distances.resize(Ntargets,Ntargets_padded);
    const size_t total_size = compute_size(Ntargets);
    memoryPool.resize(total_size*D);
    Displacements.resize(Ntargets); 
    for(int i=0; i<Ntargets; ++i)
      Displacements[i].attachReference(i,total_size,memoryPool.data()+compute_size(i));

    Temp_r.resize(Ntargets);
    Temp_dr.resize(Ntargets);
  }

  inline void evaluate(ParticleSet& P)
  {
    CONSTEXPR T BigR= std::numeric_limits<T>::max();
    //P.RSoA.copyIn(P.R); 
    for(int iat=0; iat<Ntargets; ++iat)
    {
      DTD_BConds<T,D,SC>::computeDistances(P.R[iat], P.RSoA, Distances[iat], Displacements[iat], 0, Ntargets, iat);
      Distances[iat][iat]=BigR; //assign big distance
    }
  }

  inline void evaluate(ParticleSet& P, IndexType jat)
  {
    DTD_BConds<T,D,SC>::computeDistances(P.R[jat], P.RSoA, Distances[jat], Displacements[jat], 0, Ntargets, jat);
    Distances[jat][jat]=std::numeric_limits<T>::max(); //assign a big number
  }

  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew)
  {
    DTD_BConds<T,D,SC>::computeDistances(rnew, P.RSoA, Temp_r.data(),Temp_dr, 0, Ntargets, P.activePtcl);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew)
  {
    //#pragma omp master
    moveOnSphere(P,rnew);
  }

  ///update the iat-th row for iat=[0,iat-1)
  inline void update(IndexType iat)
  {
    if(iat==0) return;
    //update by a cache line
    const int nupdate=getAlignedSize<T>(iat);
    simd::copy_n(Temp_r.data(),nupdate,Distances[iat]);
    for(int idim=0;idim<D; ++idim)
      simd::copy_n(Temp_dr.data(idim),nupdate,Displacements[iat].data(idim));
  }

};
}
#endif
