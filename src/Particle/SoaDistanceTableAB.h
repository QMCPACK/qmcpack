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

#ifndef QMCPLUSPLUS_SOADTDIMPL_AB_H
#define QMCPLUSPLUS_SOADTDIMPL_AB_H

namespace qmcplusplus
{

  template<bool MINUS>
    struct ColUpdate 
    {
      template<typename T>
        static inline void apply(int n, const T* restrict in, T* restrict out, const int nstride)
        {
          for(int i=0; i<n; ++i) out[i*nstride]=in[i];
        }
    };

  template<>
    struct ColUpdate<true>
    {
      template<typename T>
        static inline void apply(int n, const T* restrict in, T* restrict out, const int nstride)
        {
          for(int i=0; i<n; ++i) out[i*nstride]=-in[i];
        }
    };


/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense case
 */
template<typename T, unsigned D, int SC>
struct SoaDistanceTableAB: public DTD_BConds<T,D,SC>, public DistanceTableData
{

  int Nsources;
  int Ntargets;
  int Ntargets_padded;
  int BlockSize;

  ///constructor using source and target arrays
  SoaDistanceTableAB(const ParticleSet& source, ParticleSet& target)
    : DTD_BConds<T,D,SC>(source.Lattice), DistanceTableData(source,target)
  {
    N[SourceIndex]=Nsources=source.getTotalNum();
    N[VisitorIndex]=Ntargets=target.getTotalNum();
    Ntargets_padded=getAlignedSize<T>(n);
    BlockSize=Ntargets_padded*D;
    Distances.resize(Nsources,Ntargets_padded);

    memoryPool.resize(Nsources*BlockSize);
    Displacements.resize(Nsources); 
    for(int i=0; i<Nsources; ++i)
      Displacements[i].attachReference(Ntargets,Ntargets_padded,memoryPool.data()+i*BlockSize);

    Temp_r.resize(Nsources);
    Temp_dr.resize(Nsources);
  }

  SoaDistanceTableAB()=delete;
  SoaDistanceTableAB(const SoaDistanceTableAB&)=delete;
  ~SoaDistanceTableAB() {}

  inline void evaluate(ParticleSet& P)
  {
    for(int iat=0; iat<Nsources; ++iat)
      DTD_BConds<T,D,SC>::computeDistances(Origin->R[iat], P.RSoA, Distances[iat], Displacements[iat], 0, Ntargets);
  }

  inline void evaluate(ParticleSet& P, IndexType jat)
  {
    moveOnSphere(P,P.R[jat],jat);
    update(jat);
  }

  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew)
  {
    DTD_BConds<T,D,SC>::computeDistances(rnew, Origin->RSoA, Temp_r.data(),Temp_dr, 0, Nsources);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew)
  {
    moveOnSphere(P,rnew,jat);
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType iat)
  {
    CONSTEXPR bool MINUS=true;
    CONSTEXPR bool PLUS=false;

    ColUpdate<PLUS>::apply(Nsources,Temp_r.data(),Distances.data()+iat,Ntargets_padded);
    for(int idim=0;idim<D; ++idim)
      ColUpdate<MINUS>::apply(Nsources,Temp_dr.data(idim),dispAddress(idim,0,iat),BlockSize);
  }

};
}
#endif
