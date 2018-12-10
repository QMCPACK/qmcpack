//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_SYMMETRICDISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_SYMMETRICDISTANCETABLEDATAIMPL_H

namespace qmcplusplus
{

/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense-symmetric case
 *
 * SymmetricDTD stands for Symmetric Distance Table Data.
 * The source and target sets are identical and the sum is over
 * distict pairs, i.e., \f$\sum_{i}\sum_{j<i}\f$.
 * The template parameter BC provides BC::apply member function
 * to evaluate Cartesian distances.
 */
template<typename T, unsigned D, int SC>
struct SymmetricDTD
    : public DTD_BConds<T,D,SC>, public DistanceTableData
{


  //blitz::Array<IndexType,2> IJ;
  ///constructor using source and target arrays
  SymmetricDTD(const ParticleSet& source, const ParticleSet& target)
    : DTD_BConds<T,D,SC>(source.Lattice), DistanceTableData(source,target)
  {
    reset(Origin->getTotalNum(),1);
  }

  inline void reset(int m, int nactive)
  {
    if(m != N[SourceIndex] || nactive != N[WalkerIndex])
    {
      N[SourceIndex]=m;
      N[VisitorIndex]=m;
      int nn = m*(m-1)/2;
      M.resize(m+1);
      J.resize(nn);
      IJ.resize(m*m);
      PairID.resize(m*m);
      resize(nn,nactive);
      M[0] = 0;
      int nsp = Origin->groups();
      int ij = 0;
      for(int i=0; i<m; i++)
      {
        for(int j=i+1; j<m; j++, ij++)
        {
          J[ij] = j;
          PairID[ij] = Origin->GroupID[j]+nsp*Origin->GroupID[i];
          IJ[i*m+j] = ij;
          IJ[j*m+i] = ij;
          //@warning: using a simple pair scheme
          // Upper packed-storage scheme
          //PairID[ij] = Origin.GroupID[j]+(2*nsp-Origin.GroupID[i]-1)*(Origin.GroupID[i])/2;
        }
        M[i+1] = ij;
      }
      npairs_m = ij;
    }
  }


  inline virtual void nearest_neighbors(int n,int neighbors,std::vector<ripair>& ri,bool transposed=false)
  {
    int m = N[VisitorIndex];
    int shift = n*m;
    for(int i=0; i<n; ++i)
    {
      ri[i].first  = r_m[IJ[shift+i]];
      ri[i].second = i;
    }
    ri[n].first  = std::numeric_limits<RealType>::max();
    ri[n].second = n;
    shift = M[n];
    for(int i=n+1; i<m; ++i)
    {
      ri[i].first  = r_m[shift+i];
      ri[i].second = i;
    }
    partial_sort(ri.begin(),ri.begin()+neighbors,ri.end());
  }

  ///evaluate the Distance Table using a set of Particle Positions
  // inline void evaluate(const WalkerSetRef& W) {
  //   int copies = W.walkers();
  //   int visitors = W.particles();
  //   for(int iw=0; iw<copies; iw++) {
  //     int nn=0;
  //     for(int i=0; i<visitors-1; i++) {
  //       PosType a0 = W.R(iw,i);
  //       for(int j=i+1; j<visitors; j++, nn++) {
  //         PosType drij = W.R(iw,j)-a0;
  //         RealType sep = std::sqrt(BC::apply(Origin.Lattice,drij));
  //#ifdef USE_FASTWALKER
  //         r2_m(nn,iw) = sep;
  //         rinv2_m(nn,iw) = 1.0/sep;
  //         dr2_m(nn,iw) = drij;
  //#else
  //         r2_m(iw,nn) = sep;
  //         rinv2_m(iw,nn) = 1.0/sep;
  //         dr2_m(iw,nn) = drij;
  //#endif
  //       }
  //     }
  //   }
  // }

  inline void evaluate(ParticleSet& P)
  {
    const int n = N[SourceIndex];
    for(int i=0,ij=0; i<n; i++)
      for(int j=i+1; j<n; j++, ij++)
        dr_m[ij]=P.R[j]-P.R[i];
    //old with static type
    //BC::apply(Origin.Lattice,dr_m,r_m,rinv_m);
    DTD_BConds<T,D,SC>::apply_bc(dr_m,r_m,rinv_m);
  }

  inline void evaluate(ParticleSet& P, int jat)
  {
    APP_ABORT("  No need to call SymmetricDTD::evaluate(ParticleSet& P, int jat)");
    //based on full evaluation. Only compute it if jat==0
    if(jat==0) evaluate(P);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew)
  {
    for(int iat=0; iat<N[SourceIndex]; ++iat)
    {
      PosType drij(rnew - P.R[iat]);
      RealType sep=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      Temp[iat].r1=sep;
      Temp[iat].rinv1=1.0/sep;
      Temp[iat].dr1=drij;
    }
  }

  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew)
  {
    for(int iat=0; iat<N[SourceIndex]; ++iat)
    {
      PosType drij(rnew - P.R[iat]);
      Temp[iat].r1=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      Temp[iat].dr1=drij;
    }
  }

  ///update the stripe for jat-th particle
  inline void update(IndexType jat)
  {
    int nn=jat;
    for(int iat=0; iat<jat; iat++,nn+=N[SourceIndex])
    {
      int loc =IJ[nn];
      r_m[loc] = Temp[iat].r1;
      rinv_m[loc]= 1.0/Temp[iat].r1;
      //rinv_m[loc]= Temp[iat].rinv1;
      dr_m[loc]= Temp[iat].dr1;
    }
    for(int nn=M[jat]; nn<M[jat+1]; nn++)
    {
      int iat = J[nn];
      r_m[nn] = Temp[iat].r1;
      rinv_m[nn]= 1.0/Temp[iat].r1;
      //rinv_m[nn]= Temp[iat].rinv1;
      dr_m[nn]= -1.0*Temp[iat].dr1;
    }
  }
};

}
#endif
