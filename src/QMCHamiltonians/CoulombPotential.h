//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_COULOMBPOTENTIAL_H
#define QMCPLUSPLUS_COULOMBPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include <numeric>

namespace qmcplusplus
{

/** CoulombPotential
 * @tparam T type of the elementary data
 *
 * Hamiltonian operator for the Coulomb interaction for both AA and AB type for open systems.
 */
template<typename T>
struct CoulombPotential: public QMCHamiltonianBase
{
  ///true, if CoulombAA for quantum particleset
  bool is_active;
  ///distance table index, 0 indicate AA type
  int myTableIndex;
  ///number of centers
  int nCenters;
  ///source particle set
  ParticleSet* Pa;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal,1>* Va_sample;
  Array<TraceReal,1>* Vb_sample;
#endif
  ParticleSet* Pb;

  /** constructor
   * @param s source particleset
   * @param t target particleset
   * @param active if true, new Value is computed whenver evaluate is used.
   *
   * if t==0, t=s and AA interaction is used.
   */
  inline CoulombPotential(ParticleSet* s, ParticleSet* t, bool active, bool copy=false)
    : Pa(s),Pb(t),is_active(active)
  {
    set_energy_domain(potential);
    if(t)
      two_body_quantum_domain(*s,*t);
    else
      two_body_quantum_domain(*s,*s);
    nCenters=s->getTotalNum();
    if(t) // add source particle to target distance table
      myTableIndex=t->addTable(*s);//add source to the target distance table list
    else // a-a
      myTableIndex=s->addTable(*s);
    if(!is_active) //precompute the value
    {
      if(!copy) s->DistTables[0]->evaluate(*s);
      Value=evaluateAA(s->DistTables[0],s->Z.first_address());
    }
  }

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities()
  {
    request.contribute_array(myName);
  }

  virtual void checkout_particle_quantities(TraceManager& tm)
  {
    streaming_particles = request.streaming_array(myName);
    if(streaming_particles)
    {
      Va_sample = tm.checkout_real<1>(myName,*Pa);
      if(Pb)
        Vb_sample = tm.checkout_real<1>(myName,*Pb);
      else
        if(!is_active)
          evaluate_spAA(Pa->DistTables[0],Pa->Z.first_address());
    }
  }

  virtual void delete_particle_quantities()
  {
    if(streaming_particles)
    {
      delete Va_sample;
      if(Pb)
        delete Vb_sample;
    }
  }
#endif

  /** evaluate AA-type interactions */
  inline T evaluateAA(const DistanceTableData* d, const T* restrict Z)
  {
    T res=0.0;
#if !defined(REMOVE_TRACEMANAGER)
    if(streaming_particles)
      res = evaluate_spAA(d,Z);
    else
#endif
    {
      const int* restrict M=d->M.data();
      const int* restrict J=d->J.data();
      for(int iat=0; iat<nCenters; ++iat)
      {
        T q=Z[iat];
        for(int nn=M[iat]; nn<M[iat+1]; ++nn)
          res+=static_cast<RealType>(q*Z[J[nn]]*d->rinv(nn));
      }
    }
    return res;
  }


  inline T evaluateAB(const DistanceTableData* d, const T* restrict Za, const T* restrict Zb)
  {
    T res=0.0;
#if !defined(REMOVE_TRACEMANAGER)
    if(streaming_particles)
      res = evaluate_spAB(d,Za,Zb);
    else
#endif
    {
      const int* restrict M=d->M.data();
      const int* restrict J=d->J.data();
      for(int iat=0; iat<nCenters; ++iat)
      {
        T q=Za[iat];
        for(int nn=M[iat]; nn<M[iat+1]; ++nn)
          res+=static_cast<RealType>(q*Zb[J[nn]]*d->rinv(nn));
      }
    }
    return res;
  }


#if !defined(REMOVE_TRACEMANAGER)
  /** evaluate AA-type interactions */
  inline T evaluate_spAA(const DistanceTableData* d, const T* restrict Z)
  {
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    T res=0.0;
    RealType pairpot;
    Array<RealType,1>& Va_samp = *Va_sample;
    Va_samp = 0.0;
    for(int iat=0; iat<nCenters; ++iat)
    {
      T q=Z[iat];
      for(int nn=M[iat],it=0; nn<M[iat+1]; ++nn,it++)
      {
        pairpot = static_cast<RealType>(.5*q*Z[J[nn]]*d->rinv(nn));
        Va_samp(iat)+=pairpot;
        Va_samp(it) +=pairpot;
        res += 2.0*pairpot;
      }
    }
#if defined(TRACE_CHECK)
    T Vnow = res;
    T Vsum = Va_samp.sum();
    T Vorig = evaluateAA_orig(d,Z);
    if(abs(Vsum-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"accumtest: CoulombPotential::evaluateAA()"<<endl;
      app_log()<<"accumtest:   tot:"<< Vnow <<endl;
      app_log()<<"accumtest:   sum:"<< Vsum <<endl;
      APP_ABORT("Trace check failed");
    }
    if(abs(Vorig-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"versiontest: CoulombPotential::evaluateAA()"<<endl;
      app_log()<<"versiontest:   orig:"<< Vorig <<endl;
      app_log()<<"versiontest:    mod:"<< Vnow <<endl;
      APP_ABORT("Trace check failed");
    }
#endif
    return res;
  }


  inline T evaluate_spAB(const DistanceTableData* d, const T* restrict Za, const T* restrict Zb)
  {
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    T res=0.0;
    RealType pairpot;
    Array<RealType,1>& Va_samp = *Va_sample;
    Array<RealType,1>& Vb_samp = *Vb_sample;
    Va_samp = 0.0;
    Vb_samp = 0.0;
    for(int iat=0; iat<nCenters; ++iat)
    {
      T q=Za[iat];
      for(int nn=M[iat],it=0; nn<M[iat+1]; ++nn,it++)
      {
        pairpot = static_cast<RealType>(.5*q*Zb[J[nn]]*d->rinv(nn));
        Va_samp(iat)+=pairpot;
        Vb_samp(it) +=pairpot;
        res += 2.0*pairpot;
      }
    }
#if defined(TRACE_CHECK)
    RealType Vnow  = res;
    RealType Vasum = Va_samp.sum();
    RealType Vbsum = Vb_samp.sum();
    RealType Vsum  = Vasum+Vbsum;
    RealType Vorig = evaluateAB_orig(d,Za,Zb);
    if(abs(Vsum-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"accumtest: CoulombPotential::evaluateAB()"<<endl;
      app_log()<<"accumtest:   tot:"<< Vnow <<endl;
      app_log()<<"accumtest:   sum:"<< Vsum <<endl;
      APP_ABORT("Trace check failed");
    }
    if(abs(Vasum-Vbsum)>TraceManager::trace_tol)
    {
      app_log()<<"sharetest: CoulombPotential::evaluateAB()"<<endl;
      app_log()<<"sharetest:   a share:"<< Vasum <<endl;
      app_log()<<"sharetest:   b share:"<< Vbsum <<endl;
      APP_ABORT("Trace check failed");
    }
    if(abs(Vorig-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"versiontest: CoulombPotential::evaluateAB()"<<endl;
      app_log()<<"versiontest:   orig:"<< Vorig <<endl;
      app_log()<<"versiontest:    mod:"<< Vnow <<endl;
      APP_ABORT("Trace check failed");
    }
#endif
    return res;
  }
#endif



  /** evaluate AA-type interactions */
  inline T evaluateAA_orig(const DistanceTableData* d, const T* restrict Z)
  {
    T res=0.0;
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    for(int iat=0; iat<nCenters; ++iat)
    {
      T q=Z[iat];
      for(int nn=M[iat]; nn<M[iat+1]; ++nn)
        res+=static_cast<RealType>(q*Z[J[nn]]*d->rinv(nn));
    }
    return res;
  }


  inline T evaluateAB_orig(const DistanceTableData* d, const T* restrict Za, const T* restrict Zb)
  {
    T res=0.0;
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    for(int iat=0; iat<nCenters; ++iat)
    {
      T q=Za[iat];
      for(int nn=M[iat]; nn<M[iat+1]; ++nn)
        res+=static_cast<RealType>(q*Zb[J[nn]]*d->rinv(nn));
    }
    return res;
  }


  void resetTargetParticleSet(ParticleSet& P)
  {
    //myTableIndex is the same
  }

  ~CoulombPotential() { }

  void update_source(ParticleSet& s)
  {
    if(myTableIndex == 0 && (s.tag() == Pa->tag() || s.parent() == Pa->tag()))
      {
        Value=evaluateAA(s.DistTables[myTableIndex],s.Z.first_address());
      }
  }
   
  inline Return_t evaluate(ParticleSet& P)
  {
    if(is_active)
    {
      if(myTableIndex)
        Value=evaluateAB(P.DistTables[myTableIndex],Pa->Z.first_address(),P.Z.first_address());
      else
        Value=evaluateAA(P.DistTables[myTableIndex],P.Z.first_address());
    }
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    if(myTableIndex)
      os << "CoulombAB source=" << Pa->getName() << endl;
    else
      os << "CoulombAA source/target " << Pa->getName() << endl;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    if(myTableIndex)
      return new CoulombPotential(Pa,&qp,true);
    else
    {
      if(is_active)
        return new CoulombPotential(&qp,0,true);
      else
        // Ye Luo April 16th, 2015
        // avoid recomputing ion-ion DistanceTable when reusing ParticleSet 
        return new CoulombPotential(Pa,0,false,true);
    }
  }
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

