//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ParticleSetProxy.h
 *
 * Helper classes to fine-grained or nested threading
 */
#ifndef QMCPLUSPLUS_PARTICLESETPROXY_H
#define QMCPLUSPLUS_PARTICLESETPROXY_H

#include <Particle/ParticleSet.h>
#include <Particle/DistanceTableData.h>
#include <Utilities/IteratorUtility.h>

namespace qmcplusplus
{

/** Distance table to handle multiple displacements of one active particle
 *
 * j is the index of the target particle that are replaced by \f$\{\delta\}\f$
 * - r_m(i,k) distance \f$|{\bf r}_j+\delta_k-{\bf r}_i|\f$
 * - dr_m(i,k) displacement vector \f${\bf r}_j+\delta_k-{\bf r}_i \f$
 */
struct DistTableOnSphere
{
  enum
  {
    SourceIndex=DistanceTableData::SourceIndex,
    VisitorIndex=DistanceTableData::VisitorIndex
  };

  typedef DistanceTableData::RealType RealType;
  typedef DistanceTableData::PosType PosType;

  Matrix<RealType> r_m;
  Matrix<PosType>  dr_m;

  inline DistTableOnSphere(const DistanceTableData& dt)
  {
    r_m.resize(dt.centers(),12);
    dr_m.resize(dt.centers(),12);
  }

  inline void resize(int nknot)
  {
    this->resize(r_m.rows(), nknot);
  }

  inline void resize(int n, int nknot)
  {
    if(n)
    {
      r_m.resize(n,nknot);
      dr_m.resize(n,nknot);
    }
    else
      APP_ABORT("DistTableOnSphere initialized with an empty particleset");
  }
  virtual void makeMoves(int iel, const DistanceTableData& dt, const vector<PosType>& displs)=0;
};

/** proxy class to handle multiple moves
 *
 * Used by Non-local PP
 * This object cannot modify the target ParticleSet or its DistTables.
 */
struct MultiMoveHandler
{

  typedef ParticleSet::ParticlePos_t ParticlePos_t;
  typedef ParticleSet::SingleParticlePos_t PosType;

  ///ID of the moved particle
  int ActivePtcl;
  ///species index for the active particle
  int GroupID;
  ///new positions of activePtcl
  ParticlePos_t R;
  ///matching tables
  vector<DistTableOnSphere*>  myTables;

  MultiMoveHandler(const ParticleSet& pt);
  ~MultiMoveHandler();

  inline void resize(int nknot)
  {
    R.resize(nknot);
    for(int t=0; t<myTables.size(); ++t)
      myTables[t]->resize(nknot);
  }

  template<typename PA>
  inline void makeMoves(ParticleSet& P, int iat, const PA& displs)
  {
    ActivePtcl=iat;
    GroupID=myPtcl->GroupID[iat];
    for(int k=0; k<displs.size(); ++k)
      R[k]=P.R[iat]+displs[k];
    for(int t=0; t<myTables.size(); ++t)
      myTables[t]->makeMoves(iat, *(P.DistTables[t]), displs);
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5545 $   $Date: 2012-07-10 08:04:38 -0500 (Tue, 10 Jul 2012) $
 * $Id: ParticleSetProxy.h 5545 2012-07-10 13:04:38Z jnkim $
 ***************************************************************************/
