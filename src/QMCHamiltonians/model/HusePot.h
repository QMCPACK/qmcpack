//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HUSEPOT_H
#define QMCPLUSPLUS_HUSEPOT_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *@brief HusePot for the indentical source and target particle sets.
 *
 */
struct HusePot: public QMCHamiltonianBase
{

  RealType m,V;
  RealType pf,K,L;
  RealType root3;

  HusePot(ParticleSet& P):m(0),V(0),pf(0),K(0),L(0)
  {
    root3=std::sqrt(3);
    const DistanceTableData* d_table = DistanceTable::add(P,DT_AOS);
  }

  ~HusePot() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    //d_table = DistanceTable::add(P);
    //PtclRef=&P;
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    const DistanceTableData* d_table=P.DistTables[0];
    Value = 0.0;
    for(int nn=0; nn<d_table->getTotNadj(); ++nn)
    {
      Value += f(d_table->r(nn));
    }
    return Value*=pf;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t
  registerData(ParticleSet& P, BufferType& buffer)
  {
    NewValue=evaluate(P);
    buffer.add(Value);
    return Value;
  }

  inline Return_t
  updateBuffer(ParticleSet& P, BufferType& buffer)
  {
    NewValue=evaluate(P);
    buffer.put(Value);
    return Value;
  }

  inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.get(Value);
    NewValue=Value;
  }

  inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
  {
    buffer.put(Value);
  }

  inline Return_t
  evaluatePbyP(ParticleSet& P, int active)
  {
    APP_ABORT("HusePot::evaluatePbyP");
    return 0.0;
    //const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[0]->Temp);
    //Return_t del=0.0;
    //for(int iat=0; iat<Centers; ++iat) del+=(temp[iat].rinv1-temp[iat].rinv0);
    //return NewValue=Value+Q*del;
  }

  inline Return_t f(Return_t r)
  {
    Return_t x=r-root3;
    if (x>0)
      return 0;
    Return_t x3=x*x*x;
    Return_t x4=x3*x;
    Return_t x5=x4*x;
    return L*x5+K*x4-x3;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet Tattrib;
    Tattrib.add(m,"m");
    Tattrib.add(V,"V");
    Tattrib.put(cur);
    app_log()<<"HusePot parameters"<< std::endl;
    app_log()<<"  m: "<<m<< std::endl;
    app_log()<<"  V: "<<V<< std::endl;
    if(m > 0.1*std::pow(root3-1,3))
      APP_ABORT("m max is 0.1*std::pow(root3-1,3) ~ 0.0392304845 ");
    L=4*m*std::pow((root3-1),-5)-std::pow(root3-1,-2);
    K=5*m*std::pow((root3-1),-4)-2.0/(root3-1);
    Return_t rc= root3-0.6*(root3-1)/(1-m*4*std::pow(root3-1,-3));
    Return_t f_rc =1.0/f(rc);
    app_log()<<"  Huse H: "<<f(1)*f_rc<< std::endl;
    pf=V*f_rc;
    return true;
  }

  bool get(std::ostream& os) const
  {
    //os << "HusePot: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    HusePot* cl = new HusePot(qp);
    cl->pf=pf;
    cl->L=L;
    cl->K=K;
    return cl;
  }

  //#ifdef USE_FASTWALKER
  //    inline void
  //    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
  //      std::vector<ValueType> e(W.walkers(),0.0);
  //      for(int nn = 0; nn< d_table->getTotNadj(); nn++) {
  //	for(int iw=0; iw<W.walkers(); iw++) {
  //	  e[iw] += d_table->rinv(iw,nn);
  //	}
  //      }
  //      for(int iw=0; iw<W.walkers(); iw++) { LE[iw] += C*e[iw];}
  //    }
  //#else
  //    inline void
  //    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
  //      for(int iw=0; iw<W.walkers(); iw++) {
  //	RealType e =0.0;
  //	for(int nn = 0; nn< d_table->getTotNadj(); nn++)
  //      e += d_table->rinv(iw,nn);
  //	LE[iw] += C*e;
  //      }
  //    }
  //#endif
};
}
#endif


