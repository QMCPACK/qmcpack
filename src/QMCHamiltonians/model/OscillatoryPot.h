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
    
    
#ifndef QMCPLUSPLUS_OSCILLPOT_H
#define QMCPLUSPLUS_OSCILLPOT_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 *@brief OscillatoryPotential for the indentical source and target particle sets.
 *
 */
struct OscillatoryPotential: public QMCHamiltonianBase
{
  RealType v0, k0, r0, r1, rm0, rm1;


  OscillatoryPotential(ParticleSet& P)
  {
    v0=-1.0;
    k0=1.0;
    r0=1.0;
    r1=1.0;
    const DistanceTableData* d_table = DistanceTable::add(P,DT_AOS);
  }

  ~OscillatoryPotential() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    //d_table = DistanceTable::add(P,DT_AOS);
    //PtclRef=&P;
  }

  inline Return_t
  evaluate(ParticleSet& P)
  {
    const DistanceTableData* d_table=P.DistTables[0];
    Value = 0.0;
    for(int nn=0; nn<d_table->getTotNadj(); ++nn)
    {
      Return_t x(d_table->r(nn));
      Return_t x2(x*x);
      Value += std::cos(k0*x)*std::exp(-x2*rm0)/std::sqrt(rm1*x2+1);
    }
    Value*=v0;
    return Value;
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
    APP_ABORT("OscillatoryPotential::evaluatePbyP");
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
    Tattrib.add(v0,"v0");
    Tattrib.add(r0,"r0");
    Tattrib.add(r1,"r1");
    Tattrib.add(k0,"k0");
    Tattrib.put(cur);
    rm0=1.0/(r0*r0);
    rm1=1.0/(r1*r1);
    return true;
  }

  bool get(std::ostream& os) const
  {
    //os << "OscillatoryPotential: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    OscillatoryPotential* cl = new OscillatoryPotential(*this);
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


