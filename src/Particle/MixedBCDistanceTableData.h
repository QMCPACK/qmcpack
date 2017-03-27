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
    
    



#ifndef QMCPLUSPLUS_MIXEDBC_DISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_MIXEDBC_DISTANCETABLEDATAIMPL_H

namespace qmcplusplus
{

/** class to handle box condition
 * @tparam T type of values
 * @tparam D dimension
 * @tparam ORTHO true, if orthorhombic
 *
 * Generic declaration to be specialized
 */
template<typename T, unsigned D, bool ORTHO> struct BoxBCHandler {};

/** specialization for 3D orthorhombic cell */
template<typename T>
struct BoxBCHandler<T,3,true>
{
  T L0, L0Inv;
  T L1, L1Inv;
  T L2, L2Inv;
  BoxBCHandler(const CrystalLattice<T,3>& lat)
    : L0(lat.R(0,0)),L1(lat.R(1,1)),L2(lat.R(2,2))
  {
    L0Inv=1.0/L0;
    L1Inv=1.0/L1;
    L2Inv=1.0/L2;
  }

  inline void put2box(const TinyVector<T,3>& pin, TinyVector<T,3>& pout)
  {
    pout[0]=pin[0]*L0Inv;
    pout[0]=L0*(pout[0]-std::floor(pout[0]));
    pout[1]=pin[1]*L1Inv;
    pout[1]=L1*(pout[1]-std::floor(pout[1]));
    pout[2]=pin[2]*L2Inv;
    pout[2]=L2*(pout[2]-std::floor(pout[2]));
  }

  template<typename PAT1, typename PAT2>
  inline void put2box(const PAT1& pin, PAT2& pout, int n)
  {
    for(int i=0; i<n; ++i)
      put2box(pin[i],pout[i]);
  }
};

/** specialization for 3D general*/
template<typename T>
struct BoxBCHandler<T,3,false>
{
  CrystalLattice myLattice;
  BoxBCHandler(const CrystalLattice<T,3>& lat)
    : myLattice(lat)
  { }

  inline void put2box(const TinyVector<T,3>& pin, TinyVector<T,3>& pout)
  {
    TinyVector<T,3> ru(myLattice.toUnit(pin));
    ru[0]-=std::floor(pout[0]);
    ru[1]-=std::floor(pout[1]);
    ru[2]-=std::floor(pout[2]);
    pout=myLattice.toCart(ru);
  }

  template<typename PAT1, typename PAT2>
  inline void put2box(const PAT1& pin, PAT2& pout, int n)
  {
    for(int i=0; i<n; ++i)
      put2box(pin[i],pout[i]);
  }
};

template<typename T, unsigned D, bool ORTHO>
struct MixedBCDTD
    : public BoxBCHandler<T,D,ORTHO>, public DistanceTableData
{
  const ParticleSet& Target;
  std::vector<PosType> RinBox;

  MixedBCDTD(const ParticleSet& source,
             const ParticleSet& target)
    : BoxBCHandler<T,D,ORTHO>(target.Lattice)
    , DistanceTableData(source,target) , Target(target)
  {
    create(1);
  }

  void create(int walkers)
  {
    int nw = (walkers>0)? walkers:1;
    reset(Origin.getTotalNum(),Target.getTotalNum(),nw);
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
            PairID[ij] = Origin.GroupID[i];
          }
          M[i+1] = M[i]+n2;
        }
        npairs_m = n1*n2;
      }
      RinBox.resize(n2);
    }
  }

  ///not so useful inline but who knows
  inline void evaluate(const ParticleSet& P)
  {
    BoxBCHandler<T,D,ORTHO>::put2box(P.R,RinBox,N[VistiorIndex]);
    for(int i=0,ij=0; i<N[SourceIndex]; i++)
      for(int j=0; j<N[VisitorIndex]; j++,ij++)
      {
        dr_m[ij]=RinBox[j]-Origin.R[i];
        r_m[ij]=std::sqrt(dr_m[ij],dr_m[ij]);
        rinv_m[ij]=1.0/r_m[ij];
      }
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat)
  {
    activePtcl=jat;
    PosType tpos;
    BoxBCHandler<T,D,ORTHO>::put2box(rnew,tpos);
    for(int iat=0, loc=jat; iat<N[SourceIndex]; iat++,loc+=N[VisitorIndex])
    {
      PosType drij(tpos-Origin.R[iat]);
      Temp[iat].r1=std::sqrt(dot(drij,drij));
      Temp[iat].rinv1=1.0/Temp[iat].r1;
      Temp[iat].dr1=drij;
    }
  }

  ///evaluate the temporary pair relations
  inline void moveby(const ParticleSet& P, const PosType& displ, IndexType jat)
  {
  }


  ///evaluate the temporary pair relations
  inline void moveOnSphere(const ParticleSet& P, const PosType& displ, IndexType jat)
  {
  }

  inline void update(IndexType jat)
  {
    for(int iat=0,loc=jat; iat<N[SourceIndex]; iat++, loc += N[VisitorIndex])
    {
      r_m[loc]=Temp[iat].r1;
      rinv_m[loc]=1/Temp[iat].r1;
      dr_m[loc]=Temp[iat].dr1;
    }
  }
};

}
#endif
