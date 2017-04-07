//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GENERIC_TWOBODYJASTROW_SHARED_H
#define QMCPLUSPLUS_GENERIC_TWOBODYJASTROW_SHARED_H
#include "OhmmsPETE/OhmmsMatrix.h"
#include <numeric>
namespace qmcplusplus
{
/** @ingroup OrbitalComponent
 * @brief Specialization for two-body Jastrow function using one functor
 *
 *Identical function \f$ u(r_{ij})\f$ for all the pair types.
 */
template<class FT>
class TwoBodyJastrow<FT,true>: public OrbitalBase
{

  const DistanceTableData* d_table;

  int N,NN;
  ValueType DiffVal, DiffValSum;
  ValueVectorType U,d2U,curLap,curVal;
  GradVectorType dU,curGrad;
  ValueType *FirstAddressOfdU, *LastAddressOfdU;

public:

  typedef FT FuncType;
  FT F;

  ///constructor
  TwoBodyJastrow(ParticleSet& p, DistanceTableData* dtable): d_table(dtable)
  {
    N=p.getTotalNum();
    NN=N*N;
    U.resize(NN+1);
  }

  ~TwoBodyJastrow()
  {
    DEBUGMSG("TwoBodyJastrow::~TwoBodyJastrow")
  }


  ///reset the value of the Two-Body Jastrow functions
  void reset()
  {
    F.reset();
  }

  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P)
  {
    LOGMSG("TwoBodyJastrow::targetParticleSet")
    d_table = DistanceTable::getTable(DistanceTable::add(P));
  }

  /** implements the virtual functions of OrbitalBase
  @param P the particle set
  @param G returns the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
  @param L returns the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$
  @return \f$exp(-J({\bf R}))\f$
  @note The DistanceTableData contains only distinct pairs of the
  particles belonging to one set, e.g., SymmetricDTD.
  */
  inline ValueType evaluateLog(ParticleSet& P,
                               ParticleSet::ParticleGradient_t& G,
                               ParticleSet::ParticleLaplacian_t& L)
  {
    LogValue=0.0;
    ValueType dudr, d2udr2;
    PosType gr;
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //LogValue -= F.evaluate(d_table->r(nn), dudr, d2udr2);
        ValueType uij = F.evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue -= uij;
        U[i*N+j]=uij;
        U[j*N+i]=uij;
        //multiply 1/r
        dudr *= d_table->rinv(nn);
        gr = dudr*d_table->dr(nn);
        //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
        ValueType lap = d2udr2+2.0*dudr;
        //multiply -1
        G[i] += gr;
        G[j] -= gr;
        L[i] -= lap;
        L[j] -= lap;
      }
    }
    return LogValue;
  }

  inline ValueType evaluate(ParticleSet& P,
                            ParticleSet::ParticleGradient_t& G,
                            ParticleSet::ParticleLaplacian_t& L)
  {
    return exp(evaluateLog(P,G,L));
  }

  /** evalaute ratio only
   *@param P the active particle set
   *@param iat the index of the particle that is moved
   */
  ValueType ratio(ParticleSet& P, int iat)
  {
    ValueType d(0.0);
    for(int jat=0, ij=iat*N; jat<N; jat++,ij++)
    {
      if(iat != jat)
      {
        d += U[ij]-F.evaluate(d_table->Temp[jat].r1);
      }
    }
    return exp(d);
    //for(int jat=0; jat<N; jat++) {
    //  if(iat != jat) {
    //    d += F.evaluate(d_table->Temp[jat].r0) -F.evaluate(d_table->Temp[jat].r1);
    //  }
    //}
    //return exp(d);
  }

  /** later merge the loop */
  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
    register ValueType dudr, d2udr2,u;
    register PosType gr;
    DiffVal = 0.0;
    for(int jat=0, ij=iat*N; jat<N; jat++,ij++)
    {
      if(jat==iat)
      {
        curVal[jat] = 0.0;
        curGrad[jat]=0.0;
        curLap[jat]=0.0;
      }
      else
      {
        curVal[jat] = F.evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
        dudr *= d_table->Temp[jat].rinv1;
        curGrad[jat] = -dudr*d_table->Temp[jat].dr1;
        curLap[jat] = -(d2udr2+2.0*dudr);
        DiffVal += (U[ij]-curVal[jat]);
      }
    }
    GradType sumg,dg;
    ValueType suml=0.0,dl;
    for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N)
    {
      sumg += (dg=curGrad[jat]-dU[ij]);
      suml += (dl=curLap[jat]-d2U[ij]);
      dG[jat] -= dg;
      dL[jat] += dl;
    }
    dG[iat] += sumg;
    dL[iat] += suml;
    return exp(DiffVal);
  }

  /** later merge the loop */
  ValueType logRatio(ParticleSet& P, int iat,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL)
  {
    register ValueType dudr, d2udr2,u;
    register PosType gr;
    DiffVal = 0.0;
    for(int jat=0, ij=iat*N; jat<N; jat++,ij++)
    {
      if(jat==iat)
      {
        curVal[jat] = 0.0;
        curGrad[jat]=0.0;
        curLap[jat]=0.0;
      }
      else
      {
        curVal[jat] = F.evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
        dudr *= d_table->Temp[jat].rinv1;
        curGrad[jat] = -dudr*d_table->Temp[jat].dr1;
        curLap[jat] = -(d2udr2+2.0*dudr);
        DiffVal += (U[ij]-curVal[jat]);
      }
    }
    GradType sumg,dg;
    ValueType suml=0.0,dl;
    for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N)
    {
      sumg += (dg=curGrad[jat]-dU[ij]);
      suml += (dl=curLap[jat]-d2U[ij]);
      dG[jat] -= dg;
      dL[jat] += dl;
    }
    dG[iat] += sumg;
    dL[iat] += suml;
    return DiffVal;
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
    DiffValSum += DiffVal;
    for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N)
    {
      dU[ij]=curGrad[jat];
      dU[ji]=-1.0*curGrad[jat];
      d2U[ij]=d2U[ji] = curLap[jat];
      U[ij] =  U[ji] = curVal[jat];
    }
  }

  /* evalaute ratio
   *@param P the active particle set
   *@param iat the index of the particle that is moved
   *
   *@note
   *DiffVal=\f$\sum_{j}(\phi(r_{ij}^0) - \phi(r_{ij}))\f$
   *The negative sign is taken into account in the expression.
   *JK tried to remove negations as much as possible (probably waste of time).
   *This leads to a rather messy +/- convetions.
  ValueType ratio(ParticleSet& P, int iat) {
    register ValueType dudr, d2udr2,u;
    DiffVal = 0.0;
    for(int jat=0, ij=iat*N; jat<N; jat++,ij++) {
  if(jat!=iat) {
  DiffVal += U[ij]-F.evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
  }
    }
    return exp(DiffVal);
  }
   */

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat)
  {
    DiffValSum += DiffVal;
    GradType sumg,dg;
    ValueType suml=0.0,dl;
    for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N)
    {
      sumg += (dg=curGrad[jat]-dU[ij]);
      suml += (dl=curLap[jat]-d2U[ij]);
      dU[ij]=curGrad[jat];
      dU[ji]=-1.0*curGrad[jat];
      d2U[ij]=d2U[ji] = curLap[jat];
      U[ij] =  U[ji] = curVal[jat];
      dG[jat] -= dg;
      dL[jat] += dl;
    }
    dG[iat] += sumg;
    dL[iat] += suml;
  }

  /** equivalent to evalaute but the main function is to store data */
  ValueType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    N=d_table->size(VisitorIndex);
    NN=N*N;
    U.resize(NN+1);
    d2U.resize(NN);
    dU.resize(NN);
    curGrad.resize(N);
    curLap.resize(N);
    curVal.resize(N);
    ValueType dudr, d2udr2,u;
    LogValue=0.0;
    PosType gr;
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //ValueType sumu = F.evaluate(d_table->r(nn));
        u = F.evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue -= u;
        dudr *= d_table->rinv(nn);
        gr = dudr*d_table->dr(nn);
        //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
        ValueType lap = d2udr2+2.0*dudr;
        int ij = i*N+j, ji=j*N+i;
        U[ij]=u;
        U[ji]=u;
        dU[ij] = gr;
        dU[ji] = -1.0*gr;
        d2U[ij] = -lap;
        d2U[ji] = -lap;
        //add gradient and laplacian contribution
        P.G[i] += gr;
        P.G[j] -= gr;
        P.L[i] -= lap;
        P.L[j] -= lap;
      }
    }
    U[NN]= LogValue;
    buf.add(U.begin(), U.end());
    buf.add(d2U.begin(), d2U.end());
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    buf.add(FirstAddressOfdU,LastAddressOfdU);
    return LogValue;
  }

  /** equivalent to evalaute but the main function is to store data */
  ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    ValueType dudr, d2udr2,u;
    LogValue=0.0;
    PosType gr;
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //ValueType sumu = F.evaluate(d_table->r(nn));
        u = F.evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue -= u;
        dudr *= d_table->rinv(nn);
        gr = dudr*d_table->dr(nn);
        //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
        ValueType lap = d2udr2+2.0*dudr;
        int ij = i*N+j, ji=j*N+i;
        U[ij]=u;
        U[ji]=u;
        dU[ij] = gr;
        dU[ji] = -1.0*gr;
        d2U[ij] = -lap;
        d2U[ji] = -lap;
        //add gradient and laplacian contribution
        P.G[i] += gr;
        P.G[j] -= gr;
        P.L[i] -= lap;
        P.L[j] -= lap;
      }
    }
    U[NN]= LogValue;
    buf.put(U.begin(), U.end());
    buf.put(d2U.begin(), d2U.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    return LogValue;
  }

  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(U.begin(), U.end());
    buf.get(d2U.begin(), d2U.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    //ready to accumulate the differences when a move is acepted
    DiffValSum=0.0;
  }

  inline ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf)
  {
    ValueType x = (U[NN] += DiffValSum);
    buf.put(U.begin(), U.end());
    buf.put(d2U.begin(), d2U.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    return exp(x);
  }

#ifdef USE_FASTWALKER
  inline void evaluate(WalkerSetRef& W,
                       ValueVectorType& psi,
                       WalkerSetRef::WalkerGradient_t& G,
                       WalkerSetRef::WalkerLaplacian_t& L)
  {
    ValueType dudr, d2udr2;
    int nw = W.walkers();
    const DistanceTableData::IndexVectorType& M = d_table->M;
    const DistanceTableData::IndexVectorType& J = d_table->J;
    const DistanceTableData::IndexVectorType& PairID = d_table->PairID;
    for(int iw=0; iw<nw; iw++)
    {
      ValueType sumu = 0.0;
      for(int i=0; i<d_table->size(SourceIndex); i++)
      {
        for(int nn=M[i]; nn<M[i+1]; nn++)
        {
          int j = J[nn];
          sumu += F.evaluate(d_table->r(iw,nn), dudr, d2udr2);
          dudr *= d_table->rinv(iw,nn);
          PosType gr = dudr*d_table->dr(iw,nn);
          ValueType lap = d2udr2+2.0*dudr;
          G(iw,i) += gr;
          G(iw,j) -= gr;
          L(iw,i) -= lap;
          L(iw,j) -= lap;
        }
      }
      psi[iw]*= exp(-sumu);
    }
  }
#else
  inline void evaluate(WalkerSetRef& W,
                       ValueVectorType& psi,
                       WalkerSetRef::WalkerGradient_t& G,
                       WalkerSetRef::WalkerLaplacian_t& L)
  {
    ValueType dudr, d2udr2;
    int nw = W.walkers();
    const DistanceTableData::IndexVectorType& M = d_table->M;
    const DistanceTableData::IndexVectorType& J = d_table->J;
    const DistanceTableData::IndexVectorType& PairID = d_table->PairID;
    std::vector<ValueType> sumu(nw,0.0);
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=M[i]; nn<M[i+1]; nn++)
      {
        int j = J[nn];
        for(int iw=0; iw<nw; iw++)
        {
          sumu[iw] += F.evaluate(d_table->r(iw,nn), dudr, d2udr2);
          dudr *= d_table->rinv(iw,nn);
          PosType gr = dudr*d_table->dr(iw,nn);
          ValueType lap = d2udr2+2.0*dudr;
          G(iw,i) += gr;
          G(iw,j) -= gr;
          L(iw,i) -= lap;
          L(iw,j) -= lap;
        }
      }
    }
    for(int iw=0; iw<nw; iw++)
      psi[iw]*= exp(-sumu[iw]);
  }
#endif

};
}
#endif

