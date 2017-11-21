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
    
    
#ifndef QMCPLUSPLUS_GENERIC_THREEBODYJASTROWORBITAL_H
#define QMCPLUSPLUS_GENERIC_THREEBODYJASTROWORBITAL_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include <fstream>
#include <iostream>

namespace qmcplusplus
{

template<class FT>
class ThreeBodyJastrowOrbital: public OrbitalBase
{
  const ParticleSet& CenterRef;
  ///table index
  IndexType myTableIndex;
  ///size of centers
  IndexType NumCenters;
  ///size of quantum particle
  IndexType NumPtcls;
  ///temporary array to store phi(k,i)=Func[k](r_ik)
  Matrix<RealType> phi2;
  Matrix<PosType> dphi2;
  Matrix<RealType> d2phi2;

  RealType curVal, curLap;
  PosType curGrad;
  ValueVectorType U,d2U;
  GradVectorType dU;
  ValueType *FirstAddressOfdU, *LastAddressOfdU;
  std::vector<FT*> Funique;
  std::vector<FT*> Fs;
  Array<RealType,3> C;
public:

  ///constructor
  ThreeBodyJastrowOrbital(ParticleSet& ions, ParticleSet& els):
    CenterRef(ions),
    FirstAddressOfdU(0), LastAddressOfdU(0)
  {
    if(ions.tag() == els.tag())
      myTableIndex=0;
    else
      myTableIndex=els.addTable(ions);
    NumCenters=ions.getTotalNum();
    NumPtcls=els.getTotalNum();
    C.resize(NumCenters,NumPtcls,NumPtcls);
    C=0.0;
    for(int k=0; k<NumCenters; ++k)
      for(int i=0; i<NumPtcls; ++i)
        for(int j=0; j<NumPtcls; ++j)
          C(k,i,j)=(i==j)?0.0:1.0;
    phi2.resize(NumCenters,NumPtcls);
    dphi2.resize(NumCenters,NumPtcls);
    d2phi2.resize(NumCenters,NumPtcls);
    Funique.resize(ions.getSpeciesSet().getTotalNum(),0);
    Fs.resize(NumCenters,0);
  }

  ~ThreeBodyJastrowOrbital() { }

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  void checkInVariables(opt_variables_type& active)
  {
  }

  void checkOutVariables(const opt_variables_type& active)
  {
  }

  void resetParameters(const opt_variables_type& active)
  {
  }


  ValueType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    LogValue=0.0;
    const DistanceTableData& d2(*P.DistTables[0]);
    const DistanceTableData& d1(*P.DistTables[myTableIndex]);
    phi2=0.0;
    dphi2=0.0;
    d2phi2=0.0;
    RealType dudr,d2udr2;
    for(int k=0; k<NumCenters; ++k)
    {
      FT* func=FT[k];
      for(int nn=d1.M[k]; nn<d1.M[k+1]; ++nn)
      {
        int j=d1.j[nn];
        phi2(k,j)=func->evaluate(d1->r(nn),dudr,d2udr2);
        dudr *= d_table->rinv(nn);
        dphi2(k,j)=dudr*d1->dr(nn);
        d2phi2(k,j)=d2udr2+2.0*dudr;
      }
    }
    if(myTableIndex ==0)
    {
      //symmetrize phi2, dphi2, d2phi2
    }
    const RealType* restrict cptr=C.data();
    LogValue=0.0;
    for(int k=0; k<NumCenters; ++k)
      for(int i=0; i<NumPtcls; ++i)
      {
        RealType v_ki=phi2(k,j);
        PosType dv_ki=dphi2(k,j);
        RealType d2v_ki=d2phi2(k,j);
        for(int j=0; j<NumPtcls; ++j)
        {
          RealType c=*cptr++;
          RealType v_kj=phi2(k,j);
          LogValue+=c*v_ki*v_kj;
          G[i]+=c*dv_ki*v_kj;
          L[i]+=c*d2v_ki*v_kj;
          G[j]+=c*v_ki*dphi2(k,j);
          L[j]+=c*v_ki*d2phi2(k,j);
        }
      }
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return exp(evaluateLog(P,G,L));
  }

  /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
   * @param P active particle set
   * @param iat particle that has been moved.
   */
  inline ValueType ratio(ParticleSet& P, int iat)
  {
    return 0;
  }

  inline void restore(int iat) {}

  void addFunc(int source_type, FT* afunc)
  {
    if(Funique[source_type])
    {
      APP_ABORT("  ThreeBodyJastrowOrbital duplicate functor ");
    }
    else
    {
      Funique[source_type]=afunc;
      for(int i=0; i<Fs.size(); i++)
        if(CenterRef.GroupID[i] == source_type)
          Fs[i]=afunc;
    }
  }


  /** equivalent to evalaute with additional data management */
  void registerData(ParticleSet& P, WFBufferType& buf)
  {
  }

  ValueType updateBuffer(ParticleSet& P, WFBufferType& buf)
  {
    return LogValue;
  }

  /** copy the current data from a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   *
   *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
  }

  /** return the current value and copy the current data to a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   */
  inline ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf)
  {
  }

  void acceptMove(ParticleSet& P, int iat)
  {
  }

};
}
#endif
