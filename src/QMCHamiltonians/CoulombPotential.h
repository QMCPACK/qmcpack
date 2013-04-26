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
  ParticleSet* PtclA;

  /** constructor
   * @param s source particleset
   * @param t target particleset
   * @param quantum if true, new Value is computed whenver evaluate is used.
   *
   * if t==0, t=s and AA interaction is used.
   */
  inline CoulombPotential(ParticleSet* s, ParticleSet* t, bool quantum)
    : PtclA(s),is_active(quantum)
  {
    nCenters=s->getTotalNum();
    if(t) // add source particle to target distance table
      myTableIndex=t->addTable(*s);//add source to the target distance table list
    else // a-a
      myTableIndex=s->addTable(*s);
    if(!is_active) //precompute the value
    {
      s->DistTables[0]->evaluate(*s);
      Value=evaluateAA(s->DistTables[0],s->Z.first_address());
    }
  }

  /** evaluate AA-type interactions */
  inline T evaluateAA(const DistanceTableData* d, const T* restrict Z)
  {
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    T res=0.0;
    for(int iat=0; iat<nCenters; ++iat)
    {
      T q=Z[iat];
      for(int nn=M[iat]; nn<M[iat+1]; ++nn)
        res+=static_cast<RealType>(q*Z[J[nn]]*d->rinv(nn));
    }
    return res;
  }

  inline T evaluateAB(const DistanceTableData* d, const T* restrict Za, const T* restrict Zb)
  {
    const int* restrict M=d->M.data();
    const int* restrict J=d->J.data();
    T res=0.0;
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

  inline Return_t evaluate(ParticleSet& P)
  {
    if(is_active)
    {
      if(myTableIndex)
        Value=evaluateAB(P.DistTables[myTableIndex],PtclA->Z.first_address(),P.Z.first_address());
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
      os << "CoulombAB source=" << PtclA->getName() << endl;
    else
      os << "CoulombAA source/target " << PtclA->getName() << endl;
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    if(myTableIndex)
      return new CoulombPotential(PtclA,&qp,true);
    else
    {
      if(is_active)
        return new CoulombPotential(&qp,0,true);
      else
        return new CoulombPotential(PtclA,0,false);
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

