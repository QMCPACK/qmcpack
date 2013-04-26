//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_GAUSSIANPOT_H
#define QMCPLUSPLUS_GAUSSIANPOT_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *@brief GaussianPot for the indentical source and target particle sets.
 *
 */
struct GaussianPot: public QMCHamiltonianBase
{

  ///number of particle
  int Centers;
//     core radius
  RealType d;
//     core strength
  RealType Q;

  GaussianPot(ParticleSet& P)
  {
    Centers=P.getTotalNum();
    d = 1.0;
    Q=1.0;
    const DistanceTableData* d_table = DistanceTable::add(P);
  }

  ~GaussianPot() { }

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
      Value +=  std::exp(-1.0*(d_table->r(nn)*d)*(d_table->r(nn)*d));
    }
    return Value*=Q;
    //return C*std::accumulate(d_table->rinv.data(),
    //	  	       d_table->rinv.data()+d_table->getTotNadj(),
    //		       0.0);
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
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
    APP_ABORT("GaussianPot::evaluatePbyP");
    return 0.0;
    //const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[0]->Temp);
    //Return_t del=0.0;
    //for(int iat=0; iat<Centers; ++iat) del+=(temp[iat].rinv1-temp[iat].rinv0);
    //return NewValue=Value+Q*del;
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet Tattrib;
    Tattrib.add(d,"length");
    Tattrib.add(Q,"mag");
    Tattrib.put(cur);
    app_log()<<"GaussianPot parameters"<<endl;
    app_log()<<"  length: "<<d<<endl;
    d=1.0/d;
    app_log()<<"  mag: "<<Q<<endl;
    return true;
  }

  bool get(std::ostream& os) const
  {
    //os << "GaussianPot: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    GaussianPot* cl = new GaussianPot(qp);
    cl->d=d;
    cl->Q=Q;
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

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 4236 $   $Date: 2009-09-29 12:43:33 -0500 (Tue, 29 Sep 2009) $
 * $Id: CoulombPotential.h 4236 2009-09-29 17:43:33Z jeongnim.kim $
 ***************************************************************************/

