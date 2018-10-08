//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_ASYMMETRICDISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_ASYMMETRICDISTANCETABLEDATAIMPL_H

namespace qmcplusplus
{

/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense-asymmetric case
 *
 * AsymmetricDTD stands for Asymmetric Distance Table Data with
 * distinct source and target sets.
 * The template parameter BC provides BC::apply member function
 * to evaluate Cartesian distances.
 */
template<typename T, unsigned D, int SC>
struct AsymmetricDTD
    : public DTD_BConds<T,D,SC>, public DistanceTableData
{
  const ParticleSet& Target;
  AsymmetricDTD(const ParticleSet& source,
                const ParticleSet& target)
    : DTD_BConds<T,D,SC>(source.Lattice), DistanceTableData(source,target)
    , Target(target)
  {
    reset(Origin->getTotalNum(),Target.getTotalNum(),1);
  }

  /*!\fn  void reset(int n1, int n2, int nactive){
   *\param n1 the number of sources
   *\param n2 the number of targets
   *\param nactive the number of copies of the targets
   *\brief Resize the internal data and assign indices
   */
  inline void reset(int n1, int n2, int nactive)
  {
    if( n1!=N[SourceIndex] || n2 != N[VisitorIndex]
        || nactive != N[WalkerIndex])
    {
      N[SourceIndex] = n1;
      N[VisitorIndex] = n2;
      int m = n1*n2;
      if(m)
      {
        M.resize(n1+1);
        J.resize(m);
        PairID.resize(m);
        resize(m,nactive);
        M[0] = 0;
        int ij = 0;
        for(int i=0; i<n1; i++)
        {
          for(int j=0; j<n2; j++, ij++)
          {
            J[ij] = j;
            PairID[ij] = Origin->GroupID[i];
          }
          M[i+1] = M[i]+n2;
        }
        npairs_m = n1*n2;
      }
    }
  }
       
  inline virtual void nearest_neighbor(std::vector<ripair>& ri,bool transposed=false) const
  {
    if(transposed)
    {
      for(int n=0; n<ri.size(); ++n)
        ri[n].first = std::numeric_limits<RealType>::max();
      const int m = N[SourceIndex];
      const int nv = N[VisitorIndex];
      int shift = 0;
      for(int i=0; i<m; ++i,shift+=nv)
        for(int n=0; n<ri.size(); ++n)
        {
          ripair& rin = ri[n];
          RealType rp = r_m[shift+n];
          if(rp<rin.first)
          {
            rin.first  = rp;
            rin.second = i;
          }
        }
    }
    else
    {
      const int m = N[VisitorIndex];
      for(int n=0; n<ri.size(); ++n)
      {
        const int shift = M[n];
        ripair& rin = ri[n];
        rin.first = std::numeric_limits<RealType>::max();
        for(int i=0; i<m; ++i)
        {
          RealType rp = r_m[shift+i];
          if(rp<rin.first)
          {
            rin.first  = rp;
            rin.second = i;
          }
        }
      }
    }
  }

  inline virtual void nearest_neighbors(int n,int neighbors,std::vector<ripair>& ri,bool transposed=false)
  {
    if(transposed)
    {
      const int m = N[SourceIndex];
      const int nv = N[VisitorIndex];
      int shift = 0;
      for(int i=0; i<m; ++i,shift+=nv)
      {
        ri[i].first  = r_m[shift+n];
        ri[i].second = i;
      }
    }
    else
    {
      const int m = N[VisitorIndex];
      const int shift = M[n];
      for(int i=0; i<m; ++i)
      {
        ri[i].first  = r_m[shift+i];
        ri[i].second = i;
      }
    }
    partial_sort(ri.begin(),ri.begin()+neighbors,ri.end());
  }

  virtual void nearest_neighbors_by_spec(int n,int neighbors,int spec_start,std::vector<ripair>& ri,bool transposed=false)
  {
    if(transposed)
    {
      const int nv = N[VisitorIndex];
      int shift = spec_start*nv;
      for(int i=0; i<ri.size(); ++i,shift+=nv)
      {
        ri[i].first  = r_m[shift+n];
        ri[i].second = i;
      }
    }
    else
    {
      const int shift = M[n]+spec_start;
      for(int i=0; i<ri.size(); ++i)
      {
        ri[i].first  = r_m[shift+i];
        ri[i].second = i;
      }
    }
    partial_sort(ri.begin(),ri.begin()+neighbors,ri.end());
  }

  ///not so useful inline but who knows
  inline void evaluate(ParticleSet& P)
  {
    const int ns=N[SourceIndex];
    const int nt=N[VisitorIndex];
    for(int i=0; i<ns; i++)
      for(int j=0; j<nt; j++)
        dr_m[i*nt+j]=P.R[j]-Origin->R[i];
    DTD_BConds<T,D,SC>::apply_bc(dr_m,r_m,rinv_m);
  }

  inline void evaluate(ParticleSet& P, int jat)
  {
    APP_ABORT("  No need to call AsymmetricDTD::evaluate(ParticleSet& P, int jat)");
    //based on full evaluation. Only compute it if jat==0
    if(jat==0) evaluate(P);
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew)
  {
    for(int iat=0; iat<N[SourceIndex]; iat++)
    {
      PosType drij(rnew-Origin->R[iat]);
      //RealType sep2(BC::apply(Origin.Lattice,drij));
      RealType sep(std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij)));
      Temp[iat].r1=sep;
      Temp[iat].rinv1=1.0/sep;
      Temp[iat].dr1=drij;
    }
  }

  ///evaluate the temporary pair relations
  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew)
  {
    for(int iat=0; iat<N[SourceIndex]; iat++)
    {
      PosType drij(rnew-Origin->R[iat]);
      Temp[iat].r1=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      Temp[iat].dr1=drij;
    }
  }

  inline void update(IndexType jat)
  {
    for(int iat=0,loc=jat; iat<N[SourceIndex]; iat++, loc += N[VisitorIndex])
    {
      r_m[loc]=Temp[iat].r1;
      rinv_m[loc]=1/Temp[iat].r1;
      //rinv_m[loc]=Temp[iat].rinv1;
      dr_m[loc]=Temp[iat].dr1;
    }
  }

  size_t get_neighbors(int iat, RealType rcut, int* restrict jid, RealType* restrict dist, PosType* restrict displ) const
  {
    size_t nn=0;
    const int nt=N[VisitorIndex];
    for(int jat=0,loc=iat*nt; jat<nt; ++jat,++loc)
    {
      RealType rij=r_m[loc];
      if(rij<rcut) 
      {//make the compact list
        jid[nn]=jat;
        dist[nn]=rij;
        displ[nn]=dr_m[loc];
        nn++;
      }
    }
    return nn;
  }
};

}
#endif
