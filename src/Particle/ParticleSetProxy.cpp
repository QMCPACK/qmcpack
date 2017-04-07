//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file ParticleSetProxy.cpp
 *
 * Helper classes for fine-grained or nested threading
 */
#include <Particle/ParticleSetProxy.h>
#include <Lattice/ParticleBConds.h>
#include <simd/simd.hpp>

namespace qmcplusplus
{

//  template <int SC>
//    struct D2OnSphere: public DistTableOnSphere
//  {
//
//    D2OnSphere(const DistanceTableData& dt):DistTableOnSphere(dt)
//    {}
//
//    inline void makeMoves(int iel, const DistanceTableData& dt, const std::vector<PosType>& displs)
//    {
//      const int nS=dt.N[SourceIndex];
//      const int nK=r_m.cols();
//      const PosType* dr_iel = dt.dr_full[iel];
//      for(int jat=0; jat<nS; ++jat)
//      {
//        PosType* restrict dr_ptr=dr_m[jat];
//        RealType* restrict r_ptr=r_m[jat];
//        for(int k=0; k<nK; ++k)
//        {
//          dr_ptr[k]=displs[k]-dr_iel[jat];
//          r_ptr[k]=std::sqrt(dot(dr_ptr[k],dr_ptr[k]));
//        }
//      }
//    }
//  };
//
//  template <int SC>
//    struct D1OnSphere: public DistTableOnSphere
//  {
//    D1OnSphere(const DistanceTableData& dt):DistTableOnSphere(dt)
//    {}
//
//    inline void makeMoves(int iel, const DistanceTableData& dt, const std::vector<PosType>& displs)
//    {
//      const int nC=dt.N[SourceIndex];
//      const int nV=dt.N[VisitorIndex];
//      const int nK=r_m.cols();
//      for(int iat=0, loc=iel; iat<nC; ++iat,loc+=nV)
//      {
//        const PosType& dr0(dt.dr(loc));
//        PosType* restrict dr_ptr=dr_m[iat];
//        RealType* restrict r_ptr=r_m[iat];
//        for(int k=0; k<nK; ++k)
//        {
//          dr_ptr[k]=dr0+displs[k];
//          r_ptr[k]=std::sqrt(dot(dr_ptr[k],dr_ptr[k]));
//        }
//      }
//    }
//  };

template <typename T, unsigned D, int SC>
struct D2OnSphere: public DistTableOnSphere
{
  DTD_BConds<T,D,SC> myBconds;

  D2OnSphere(const ParticleSet& P, const DistanceTableData& dt)
    :DistTableOnSphere(dt), myBconds(P.Lattice)
  {}

  inline void makeMoves(int iel, const DistanceTableData& dt, const std::vector<PosType>& displs)
  {
    const int nS=dt.N[SourceIndex];
    const int nK=r_m.cols();
    const PosType* dr_iel = dt.dr_full[iel];
    for(int jat=0; jat<nS; ++jat)
    {
      PosType* restrict dr_ptr=dr_m[jat];
      for(int k=0; k<nK; ++k)
        dr_ptr[k]=displs[k]-dr_iel[jat];
    }
    myBconds.evaluate_rsquared(dr_m.data(),r_m.data(), r_m.size());
    simd::sqrt(r_m.data(), r_m.size());
  }
};

template <typename T, unsigned D, int SC>
struct D1OnSphere: public DistTableOnSphere
{
  DTD_BConds<T,D,SC> myBconds;

  D1OnSphere(const ParticleSet& P, const DistanceTableData& dt)
    :DistTableOnSphere(dt), myBconds(P.Lattice)
  {}

  inline void makeMoves(int iel, const DistanceTableData& dt, const std::vector<PosType>& displs)
  {
    const int nC=dt.N[SourceIndex];
    const int nV=dt.N[VisitorIndex];
    const int nK=r_m.cols();
    for(int iat=0, loc=iel; iat<nC; ++iat,loc+=nV)
    {
      const PosType& dr0(dt.dr(loc));
      PosType* restrict dr_ptr=dr_m[iat];
      for(int k=0; k<nK; ++k)
        dr_ptr[k]=displs[k]+dr0;
    }
    myBconds.evaluate_rsquared(dr_m.data(),r_m.data(), r_m.size());
    simd::sqrt(r_m.data(), r_m.size());
  }
};

/** constructor
 * @param pt MC ParticleSet
 *
 * Create DistTableOnSphere using the existing distance tables of pt
 */
MultiMoveHandler::MultiMoveHandler(const ParticleSet& pt)
  : ActivePtcl(-1), GroupID(0)
{
  myTables.reserve(pt.DistTables.size());
  myTables.push_back(new D2OnSphere<double,3,SUPERCELL_BULK>(pt, *(pt.DistTables[0])));
  for(int t=1; t<pt.DistTables.size(); ++t)
    myTables.push_back(new D1OnSphere<double,3,SUPERCELL_BULK>(pt, *(pt.DistTables[t])));
}

MultiMoveHandler::~MultiMoveHandler()
{
  delete_iter(myTables.begin(), myTables.end());
}

}
