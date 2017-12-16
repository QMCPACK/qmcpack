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
    
    
#ifndef QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_H
#define QMCPLUSPLUS_ATOMIC_HARTREEFOCK_STO_H
#include "Particle/DistanceTable.h"
#include "Numerics/SlaterTypeOrbital.h"
#include "Numerics/SphericalTensor.h"
#include "QMCWaveFunctions/DummyBasisSet.h"

namespace qmcplusplus
{

template<class T, class POS>
struct ComboSTO
{

  std::string Name;

  typedef T value_type;
  typedef SphericalTensor<T,POS> SphericalHarmonics_t;
  typedef GenericSTO<T> RadialOrbital_t;

  int LM;
  SphericalHarmonics_t& Ylm;
  std::vector<RadialOrbital_t*> Rnl;
  std::vector<T> C;

  ComboSTO(const ComboSTO& aSTO):
    Name(aSTO.Name), LM(aSTO.LM), Ylm(aSTO.Ylm), C(aSTO.C)
  {
    Rnl.resize(C.size(),NULL);
    for(int i=0; i<Rnl.size(); i++)
      Rnl[i] = aSTO.Rnl[i];
  }

  ComboSTO(int lm, SphericalHarmonics_t& ylm,
           const std::vector<RadialOrbital_t*>& rnl,
           T* c):
    LM(lm), Ylm(ylm)
  {
    C.insert(C.begin(),c, c+rnl.size());
    Rnl.resize(rnl.size(),NULL);
    for(int i=0; i<rnl.size(); i++)
      Rnl[i] = rnl[i];
  }

  /*return the value of the radial function: used by Transform2GridFunctor **/
  inline T f(T r)
  {
    T rnl = 0.0;
    for(int nl=0; nl<C.size(); nl++)
      rnl += C[nl]*Rnl[nl]->f(r);
    return rnl;
  }

  /*return the derivate of the radial function: used by Transform2GridFunctor **/
  inline T df(T r)
  {
    T drnl = 0.0;
    for(int nl=0; nl<C.size(); nl++)
      drnl += C[nl]*Rnl[nl]->df(r);
    return drnl;
  }

  inline T evaluate(T r, T rinv, T& drnl, T& d2rnl)
  {
    T rnl = 0.0;
    d2rnl = 0.0;
    drnl=0.0;
    for(int nl=0; nl<C.size(); nl++)
    {
      Rnl[nl]->evaluateAll(r,rinv);
      rnl += C[nl]*Rnl[nl]->Y;
      drnl += C[nl]*Rnl[nl]->dY;
      d2rnl += C[nl]*Rnl[nl]->d2Y;
    }
    return rnl;
  }

  /** return the value only */
  inline T evaluate()
  {
    T rnl = 0.0;
    for(int nl=0; nl<C.size(); nl++)
    {
      rnl += C[nl]*Rnl[nl]->Y;
    }
    return Ylm.getYlm(LM)*rnl;
  }

  inline T operator()(T r, T rinv, const POS& dr, POS& dy, T& d2y)
  {
    T rnl = 0.0, d2rnl = 0.0, drnl=0.0;
    for(int nl=0; nl<C.size(); nl++)
    {
      rnl += C[nl]*Rnl[nl]->Y;
      drnl += C[nl]*Rnl[nl]->dY;
      d2rnl += C[nl]*Rnl[nl]->d2Y;
    }
    T drnloverr = rinv*drnl;
    T ang = Ylm.getYlm(LM);
    POS gr_rad(drnloverr*dr);
    POS gr_ang(Ylm.getGradYlm(LM));
    dy= ang*gr_rad+rnl*gr_ang;
    d2y = ang*(2.0*drnloverr+d2rnl)+2.0*dot(gr_rad,gr_ang);
    return ang*rnl;
  }

};

struct HFAtomicSTOSet: public QMCTraits
{

  typedef SphericalTensor<RealType,PosType> SphericalHarmonics_t;
  typedef GenericSTO<RealType>              RadialOrbital_t;
  typedef ComboSTO<RealType,PosType>        SPO_t;

  typedef DummyBasisSet                     BasisSet_t;

  ///unique set of orbitals
  SphericalHarmonics_t Ylm;
  std::vector<RadialOrbital_t*> RnlPool;
  std::vector<SPO_t*> Orbital;

  ///default constructor for He: test only
  HFAtomicSTOSet();

  explicit HFAtomicSTOSet(int lmax): Ylm(lmax) { }

  ///reference to a DistanceTableData
  const DistanceTableData* d_table;

  inline void reset() { }

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(d_table->origin(),P);
  }

  inline int size() const
  {
    return Orbital.size();
  }

  template<class VV>
  inline void evaluate(const ParticleSet& P, int iat, VV& phi)
  {
    RealType r(d_table->Temp[0].r1);
    RealType rinv(d_table->Temp[0].rinv1);
    PosType dr(d_table->Temp[0].dr1);
    Ylm.evaluate(dr);
    for(int nl=0; nl<RnlPool.size(); nl++)
    {
      RnlPool[nl]->evaluate(r,rinv);
    }
    for(int j=0; j<Orbital.size(); j++)
    {
      phi[j] = Orbital[j]->evaluate();
    }
  }

  // evaluate the single-particle orbital value of iat-th el
  template<class VV, class GV>
  inline void evaluate(const ParticleSet& P, int iat, VV& phi, GV& dphi, VV& d2phi )
  {
    RealType r(d_table->Temp[0].r1);
    RealType rinv(d_table->Temp[0].rinv1);
    PosType dr(d_table->Temp[0].dr1);
    Ylm.evaluateAll(dr);
    for(int nl=0; nl<RnlPool.size(); nl++)
    {
      RnlPool[nl]->evaluateAll(r,rinv);
    }
    for(int j=0; j<Orbital.size(); j++)
    {
      phi[j] = (*Orbital[j])(r,rinv,dr,dphi[j],d2phi[j]);
    }
  }

  template<class VM, class GM>
  inline void
  evaluate(const ParticleSet& P, int first, int last,
           VM& logdet, GM& dlogdet, VM& d2logdet)
  {
    int nptcl = last-first;
    int nn = first;///first pair of the particle subset
    for(int i=0; i<nptcl; i++, nn++)
    {
      RealType r(d_table->r(nn));
      RealType rinv(d_table->rinv(nn));
      PosType dr(d_table->dr(nn));
      Ylm.evaluateAll(dr);
      for(int nl=0; nl<RnlPool.size(); nl++)
      {
        RnlPool[nl]->evaluateAll(r,rinv);
      }
      for(int j=0; j<Orbital.size(); j++)
      {
        logdet(j,i) = (*Orbital[j])(r,rinv,dr,dlogdet(i,j),d2logdet(i,j));
      }
    }
  }


//    template<class VM, class GM>
//    inline void
//    evaluate(const WalkerSetRef& W, int first, int last,
//	     std::vector<VM>& logdet, std::vector<GM>& dlogdet, std::vector<VM>& d2logdet) {
//
//#ifdef USE_FASTWALKER
//      int nptcl = last-first;
//      for(int i=0,nn=first; i<nptcl; i++, nn++) {
//        for(int iw=0; iw<W.walkers(); iw++) {
//          RealType r = d_table->r(iw,nn);
//          RealType rinv = d_table->rinv(iw,nn);
//          PosType dr = d_table->dr(iw,nn);
//          Ylm.evaluateAll(dr);
//          for(int nl=0; nl<RnlPool.size(); nl++) RnlPool[nl]->evaluateAll(r,rinv);
//          for(int j=0; j<Orbital.size(); j++)
//            logdet[iw](j,i) = (*Orbital[j])(r,rinv,dr,dlogdet[iw](i,j),d2logdet[iw](i,j));
//        }
//      }
//#else
//      int nptcl = last-first;
//      for(int iw=0; iw<W.walkers(); iw++) {
//	int nn = first;///first pair of the particle subset
//	for(int i=0; i<nptcl; i++, nn++) {
//	  RealType r = d_table->r(iw,nn);
//	  RealType rinv = d_table->rinv(iw,nn);
//	  PosType dr = d_table->dr(iw,nn);
//	  Ylm.evaluateAll(dr);
//	  for(int nl=0; nl<RnlPool.size(); nl++)  {
//	    RnlPool[nl]->evaluateAll(r,rinv);
//	  }
//	  for(int j=0; j<Orbital.size(); j++)  {
//	  logdet[iw](j,i)
//	    = (*Orbital[j])(r,rinv,dr,dlogdet[iw](i,j),d2logdet[iw](i,j));
//	  }
//	}
//      }
//#endif
//    }

  void setTable(DistanceTableData* dtable)
  {
    d_table =dtable;
  }

};
}
#endif
