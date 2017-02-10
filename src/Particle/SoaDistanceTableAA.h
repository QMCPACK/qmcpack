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
  int BlockSize;

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

  void resize(int n)
  {
    Ntargets=n;
    Ntargets_padded=getAlignedSize<T>(n);
    BlockSize=Ntargets_padded*D;
    Distances.resize(Ntargets,Ntargets_padded);

    memoryPool.resize(Ntargets*BlockSize);
    Displacements.resize(Ntargets); 
    for(int i=0; i<Ntargets; ++i)
      Displacements[i].resetByRef(Ntargets,Ntargets_padded,memoryPool.data()+i*BlockSize);

    Temp_r.resize(Ntargets);
    Temp_dr.resize(Ntargets);
  }

  inline void evaluate(ParticleSet& P)
  {
    CONSTEXPR T BigR= std::numeric_limits<T>::max();
    //P.RSoA.copyIn(P.R); 
    for(int iat=0; iat<Ntargets; ++iat)
    {
      DTD_BConds<T,D,SC>::computeDistances(P.R[iat], P.RSoA, Distances[iat], Displacements[iat], 0, Ntargets);
      Distances[iat][iat]=BigR; //assign big distance
    }
  }

  inline void evaluate(ParticleSet& P, IndexType jat)
  {
    activePtcl=jat;
    DTD_BConds<T,D,SC>::computeDistances(P.R[jat], P.RSoA, Distances[jat],Displacements[jat], 0, Ntargets);
    Distances[jat][jat]=std::numeric_limits<T>::max(); //assign a big number
  }

  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew, IndexType jat) 
  {
    DTD_BConds<T,D,SC>::computeDistances(rnew, P.RSoA, Temp_r.data(),Temp_dr, 0, Ntargets);
    Temp_r[jat]=std::numeric_limits<T>::max(); //assign a big number
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat)
  {
    //#pragma omp master
    activePtcl=jat;
    moveOnSphere(P,rnew,jat);
  }

  ///update the iat-th row for iat=[0,iat-1)
  inline void update(IndexType iat)
  {
    //if(iat==0 || iat!=activePtcl) return;
    //update by a cache line
    //const int nupdate=getAlignedSize<T>(iat);
    const int nupdate=NumTargets;
    simd::copy_n(Temp_r.data(),nupdate,Distances[iat]);
    for(int idim=0;idim<D; ++idim)
      simd::copy_n(Temp_dr.data(idim),nupdate,Displacements[iat].data(idim));
  }

};
}
#endif
