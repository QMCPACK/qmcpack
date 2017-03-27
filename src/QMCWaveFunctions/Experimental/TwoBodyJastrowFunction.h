//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GENERIC_TWOBODYJASTROW_H
#define QMCPLUSPLUS_GENERIC_TWOBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief A generic TwoBodyJastrow functions
 *
 *  This is a dummy class to enable specializations for ShareFunction=true|false.
 *The first template parameter FT is a functional \f$ u(r_{ij}) \f$, e.g.
 *PadeJastrow<T>
 *Requriement of the template function is
 *ValueType evaluate(ValueType r, ValueType& dudr, ValueType& d2udr2).
 *The second template parameter SharedFunction is a boolean.
 *SharedFunction=false means that each pair type (ij) has a unique
 *functions, while SharedFunction=true means that all pair types
 *share the same functional

  * The Two-Body Jastrow has the form
  * \f$ \psi = \exp{\left(-J_2\right)}\f$ where
   \f[ J_{2}({\bf R}) = \sum_{i<j}^N u(r_{ij}) \f]
   where \f[ r_{ij} = |{\bf r}_i - {\bf r}_j| \f]
   and the summnation is over all distinct pairs,
   as indicated by the sum \f$\sum_{i<j}\f$.
   *TwoBodyJastrow<FT,bool> is specialized for TwoBodyJastrow<FT,true>
   *and TwoBodyJastrow<FT,false>.
   *
   To calculate the Gradient use the identity
   \f[
   {\bf \nabla}_k(r_{ik}) = -\frac{{\bf r_{ik}}}{r_{ik}}.
   \f]

   \f[
   {\bf \nabla}_k(J_{ee}({\bf R})) =
   {\bf \nabla}_k \left[ \sum_{i=1}^{k-1} u(r_{ik})
   + \sum_{i=k+1}^N u(r_{ki})
   \right]
   \f]
   by application of the chain rule
   \f[
   {\bf \nabla}_k(J_{ee}({\bf R})) =
   \sum_{i=1}^{k-1} \frac{du}{dr_{ik}}({\bf \hat{r_{ik}}})+
   +\sum_{i=k+1}^N \frac{du}{dr_{ki}}({\bf \hat{r_{ki}}})
   \f]
   which finally leads to the result
   \f[
   {\bf \nabla}_k(J_{ee}({\bf R}))=
   -\sum_{i \neq k}^N \frac{du}{dr_{ik}}{\bf \hat{r_{ik}}}.
   \f]
   To calculate the Laplacian, use the identity
   \f[
   \nabla^2_k(r_{ik})=\frac{2}{r_{ik}},
   \f]
   and the vector product rule
   \f[
   \nabla^2 \cdot (f{\bf A}) =
   f({\bf \nabla} \cdot {\bf A})+{\bf A}\cdot({\bf \nabla} f)
   \f]
   \f[
   \nabla^2_k (J_{ee}({\bf R}))  =
   -\sum_{i \neq k}^N \left(\frac{du}{dr_{ik}}\right) {\bf \nabla}_k
   \cdot {\bf \hat{r_{ki}}} - {\bf \hat{r_{ki}}} \cdot
   \left(\frac{d^2u}{dr_{ik}^2}\right){\bf \hat{r_{ki}}}
   \f]
   which can be simplified to
   \f[\nabla^2_k (J_{ee}({\bf R})) = \sum_{i \neq k}^N \left(\frac{2}{r_{ik}}\frac{du}{dr_{ik}}
   + \frac{d^2u}{dr_{ik}^2}\right) \f]
 *
 *
 *@warning Note that -1 is multiplied globally after summing the pair terms.
 */
template<class FT, bool SharedFunction>
class TwoBodyJastrow: public OrbitalBase {};

/** @ingroup OrbitalComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 *
 *Each pair-type can have distinct function \f$u(r_{ij})\f$.
 *For electrons, distinct pair correlation functions are used
 *for spins up-up/down-down and up-down/down-up.
 */
template<class FT>
class TwoBodyJastrow<FT,false>: public OrbitalBase
{

  const DistanceTableData* d_table;

  int N,NN;
  ValueType DiffVal, DiffValSum;
  ValueVectorType U,d2U,curLap,curVal;
  GradVectorType dU,curGrad;
  ValueType *FirstAddressOfdU, *LastAddressOfdU;
  Matrix<int> PairID;

public:

  typedef FT FuncType;
  ///container for the Jastrow functions
  std::vector<FT*> F;

  ///constructor
  TwoBodyJastrow(ParticleSet& p, DistanceTableData* dtable): d_table(dtable)
  {
    N=p.getTotalNum();
    NN=N*N;
    U.resize(NN+1);
    PairID.resize(N,N);
    int nsp=p.groups();
    for(int i=0; i<N; i++)
      for(int j=0; j<N; j++)
        PairID(i,j) = p.GroupID[i]*nsp+p.GroupID[j];
  }

  ~TwoBodyJastrow()
  {
    DEBUGMSG("TwoBodyJastrow::~TwoBodyJastrow")
    //for(int i=0; i<F.size(); i++) delete F[i];
  }

  ///reset the value of all the Two-Body Jastrow functions
  void reset()
  {
    for(int i=0; i<F.size(); i++)
      F[i]->reset();
  }

  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::getTable(DistanceTable::add(P));
  }

  /**
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@param G returns the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
   *@param L returns the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$
   *@return \f$exp(-J({\bf R}))\f$
   *@brief While evaluating the value of the Jastrow for a set of
   *particles add the gradient and laplacian contribution of the
   *Jastrow to G(radient) and L(aplacian) for local energy calculations
   *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f]
   *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
   *@note The DistanceTableData contains only distinct pairs of the
   *particles belonging to one set, e.g., SymmetricDTD.
   */
  ValueType evaluateLog(ParticleSet& P,
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
        //LogValue -= F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        ValueType uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue -= uij;
        U[i*N+j]=uij;
        U[j*N+i]=uij; //save for the ratio
        //multiply 1/r
        dudr *= d_table->rinv(nn);
        gr = dudr*d_table->dr(nn);
        //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)
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

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return exp(evaluateLog(P,G,L));
  }

  ValueType ratio(ParticleSet& P, int iat)
  {
    ValueType d(0.0);
    const int* pairid(PairID[iat]);
    for(int jat=0, ij=iat*N; jat<N; jat++,ij++)
    {
      if(iat != jat)
      {
        d += U[ij]-F[pairid[jat]]->evaluate(d_table->Temp[jat].r1);
      }
    }
    return exp(d);
    //for(int jat=0; jat<N; jat++) {
    //  if(iat != jat) {
    //    FT *func(F[pairid[jat]]);
    //    d += func->evaluate(d_table->Temp[jat].r0) -func->evaluate(d_table->Temp[jat].r1);
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
    const int* pairid = PairID[iat];
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
        curVal[jat] = F[pairid[jat]]->evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
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
    const int* pairid = PairID[iat];
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
        curVal[jat] = F[pairid[jat]]->evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
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


  inline void update(ParticleSet& P,
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

  inline ValueType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    std::cerr <<"REGISTERING 2 BODY JASTROW "<< std::endl;
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
    PairID.resize(d_table->size(SourceIndex),d_table->size(SourceIndex));
    int nsp=P.groups();
    for(int i=0; i<d_table->size(SourceIndex); i++)
      for(int j=0; j<d_table->size(SourceIndex); j++)
        PairID(i,j) = P.GroupID[i]*nsp+P.GroupID[j];
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //ValueType sumu = F.evaluate(d_table->r(nn));
        u = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
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

  inline ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    ValueType dudr, d2udr2,u;
    LogValue=0.0;
    PosType gr;
    PairID.resize(d_table->size(SourceIndex),d_table->size(SourceIndex));
    int nsp=P.groups();
    for(int i=0; i<d_table->size(SourceIndex); i++)
      for(int j=0; j<d_table->size(SourceIndex); j++)
        PairID(i,j) = P.GroupID[i]*nsp+P.GroupID[j];
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //ValueType sumu = F.evaluate(d_table->r(nn));
        u = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
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

  inline void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(U.begin(), U.end());
    buf.get(d2U.begin(), d2U.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
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
  void evaluate(WalkerSetRef& W,
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
          sumu += F[PairID[nn]]->evaluate(d_table->r(iw,nn), dudr, d2udr2);
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
  void evaluate(WalkerSetRef& W,
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
          sumu[iw] += F[PairID[nn]]->evaluate(d_table->r(iw,nn), dudr, d2udr2);
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
#include "QMCWaveFunctions/TwoBodyJastrowFunction.Shared.h"
#endif

