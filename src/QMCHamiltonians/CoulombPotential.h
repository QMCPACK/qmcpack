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
#ifndef QMCPLUSPLUS_COULOMBPOTENTIAL_H
#define QMCPLUSPLUS_COULOMBPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include <numeric>

namespace qmcplusplus {

  /** @ingroup hamiltonian
   *@brief CoulombPotential for the different source and target particle sets.
   *
   * \f[ H = \sum_i \frac{Z(i)q}{r} \f] 
   * where \f$ Z(i) \f$ is the effective charge of the Ith 
   * ion and \f$ q \f$ is the charge of the set of quantum particles.
   * For instance, \f$ q = -1 \f$ for electrons and 
   * \f$ q = 1 \f$ for positrons.
   *
   * @warning need to be generalized by checking visitor.Species.
   */
  struct CoulombPotentialAB: public QMCHamiltonianBase {

    ///number of ions
    int Centers;
    ParticleSet& sourcePtcl;
    DistanceTableData* d_table;
    ///container for the ion charges
    vector<RealType> Z;

    CoulombPotentialAB(ParticleSet& ions, ParticleSet& els): 
      sourcePtcl(ions), d_table(0) { 
      d_table = DistanceTable::add(ions,els);
      //index for attribute charge
      SpeciesSet& tspecies(ions.getSpeciesSet());
      int iz = tspecies.addAttribute("charge");
      Centers = ions.getTotalNum();
      Z.resize(Centers);
      RealType C = -1.0; 
      for(int iat=0; iat<Centers;iat++) {
        Z[iat] = tspecies(iz,ions.GroupID[iat])*C;
      }
    }
    
    void resetTargetParticleSet(ParticleSet& P)  
    {
      d_table = DistanceTable::add(sourcePtcl,P);
    }

    ~CoulombPotentialAB() { }

    inline Return_t evaluate(ParticleSet& P) {
      Value=0.0;
      for(int iat=0; iat<Centers; ++iat) {
        for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; ++nn)
          Value+=Z[iat]*d_table->rinv(nn);
      }
      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
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
    }

    inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
    {
      buffer.put(Value);
    }

    inline Return_t 
    evaluatePbyP(ParticleSet& P, int active)
    {
      APP_ABORT("CoulombPotential::evaluatePbyP");
      return 0.0;
      //const std::vector<DistanceTableData::TempDistType> &temp(d_table->Temp);
      //Return_t del=0.0;
      //for(int iat=0; iat<Centers; ++iat) del+=Z[iat]*(temp[iat].rinv1-temp[iat].rinv0);
      //return NewValue=Value+del;
    }

    bool put(xmlNodePtr cur) 
    {
      return true;
    }

    bool get(std::ostream& os) const 
    {
      os << "CoulombPotentialAB potential: " << sourcePtcl.getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new CoulombPotentialAB(sourcePtcl, qp); 
    }

    //#ifdef USE_FASTWALKER
    //    inline void 
    //    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    //      register int nw = W.walkers();
    //      ValueVectorType e(nw);
    //      for(int iat=0; iat<Centers; iat++) {
    //        e=0.0;
    //        for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++) {
    //          for(int iw=0; iw<nw; iw++) {
    //	    e[iw] += d_table->rinv(iw,nn);
    //	  }
    //	}
    //        ValueType z=Z[iat];
    //        for(int iw=0; iw<W.walkers(); iw++) { LE[iw] += z*e[iw];} 
    //      }
    //    }
    //#else
    //    inline void 
    //    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    //      for(int iat=0; iat<Centers; iat++) {
    //        RealType z=Z[iat];
    //        for(int iw=0; iw<W.walkers(); iw++) {
    //	  RealType esub = 0.0;
    //	  for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++) {
    //	    esub += d_table->rinv(iw,nn);
    //	  }
    //	  LE[iw] += z*esub; ///multiply z
    //	}
    //      }
    //    }
    //#endif
  };

  /** @ingroup hamiltonian
   *@brief CoulombPotential for the indentical source and target particle sets. 
   *
   * \f[ H = \sum_i \frac{q^2}{r} \f] 
   * where \f$ q \f$ is the charge of the set of quantum 
   * particles.  For instance, \f$ q = -1 \f$ for electrons 
   * and \f$ q = 1 \f$ for positrons.
   */
  struct CoulombPotentialAA: public QMCHamiltonianBase {

    ///number of particle
    int Centers;
    ///Charge factor=q*q
    RealType Q;
    //DistanceTableData* d_table;
    //ParticleSet* PtclRef;
    //ElecElecPotential(RealType c=1.0): C(c){}
    //CoulombPotentialAA(ParticleSet& P):d_table(NULL),PtclRef(&P) {
    CoulombPotentialAA(ParticleSet& P) 
    {
      Centers=P.getTotalNum();
      Q = 1.0;//cheating, need to fix this
      const DistanceTableData* d_table = DistanceTable::add(P);
    }

    ~CoulombPotentialAA() { }

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
      for(int nn=0; nn<d_table->getTotNadj(); ++nn) { Value += d_table->rinv(nn); }
      return Value*=Q;
      //return C*std::accumulate(d_table->rinv.data(), 
      //	  	       d_table->rinv.data()+d_table->getTotNadj(),
      //		       0.0);
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
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
      APP_ABORT("CoulombAA::evaluatePbyP");
      return 0.0;
      //const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[0]->Temp);
      //Return_t del=0.0;
      //for(int iat=0; iat<Centers; ++iat) del+=(temp[iat].rinv1-temp[iat].rinv0);
      //return NewValue=Value+Q*del;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      //os << "CoulombPotentialAA: " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new CoulombPotentialAA(qp);
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
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

