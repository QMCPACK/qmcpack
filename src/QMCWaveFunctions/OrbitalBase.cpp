//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
//#include "QMCWaveFunctions/ProxyOrbital.h"

namespace qmcplusplus
{
OrbitalBase::OrbitalBase():
  IsOptimizing(false),Optimizable(true), UpdateMode(ORB_WALKER), //UseBuffer(true), //Counter(0),
  LogValue(1.0),PhaseValue(0.0),OrbitalName("OrbitalBase"), derivsDone(false), parameterType(0)
#if !defined(ENABLE_SMARTPOINTER)
  ,dPsi(0), ionDerivs(false)
#endif
{ 
  HaveRatiosForVP=false;
}

// OrbitalBase::OrbitalBase(const OrbitalBase& old):
//   Optimizable(old.Optimizable), UseBuffer(old.UseBuffer),
//   dPsi(old.dPsi),dLogPsi(old.dLogPsi),d2LogPsi(old.d2LogPsi),
//   OrbitalName(old.OrbitalName),myVars(old.myVars)
// {
//   //
//   //if(dLogPsi.size()) dLogPsi.resize(dLogPsi.size());
//   //if(d2LogPsi.size()) dLogPsi.resize(d2LogPsi.size());
//   //if(dPsi) dPsi=old.dPsi->makeClone();
// }


void OrbitalBase::setDiffOrbital(DiffOrbitalBasePtr d)
{
#if defined(ENABLE_SMARTPOINTER)
  dPsi=DiffOrbitalBasePtr(d);
#else
  dPsi=d;
#endif
}

void OrbitalBase::evaluateDerivatives(ParticleSet& P,
                                      const opt_variables_type& active,
                                      vector<RealType>& dlogpsi, vector<RealType>& dhpsioverpsi)
{
#if defined(ENABLE_SMARTPOINTER)
  if (dPsi.get())
#else
  if (dPsi)
#endif
    dPsi->evaluateDerivatives(P, active, dlogpsi, dhpsioverpsi);
}

///** makeClone uses optVars  to determine if it will make a clone (deepcopy)
// * or make a ProxyOrbital.
// */
//OrbitalBasePtr OrbitalBase::makeClone(ParticleSet& tpq,  int i)
//{
//  int loc=-1;
//  int iii=0;
//  while(loc<0 && iii<myVars.size())
//  {
//    if(myVars.Index[iii]==i) loc=iii;
//    ++iii;
//  }
//  if(loc<0)
//    return makeProxy(tpq,this);
//  else
//    return makeClone(tpq,true);
//}

/*@todo makeClone should be a pure virtual function
 */
OrbitalBasePtr OrbitalBase::makeClone(ParticleSet& tpq) const
{
  APP_ABORT("Implement OrbitalBase::makeClone "+OrbitalName+ " class.");
  return 0;
}

//void OrbitalBase::copyFrom(const OrbitalBase& old)
//{
//  APP_ABORT("OrbitalBase::copyFrom needs to be implemented by a derived class.");
//}

OrbitalBasePtr OrbitalBase::makeProxy(ParticleSet& tpq)
{
  return 0;
//    OrbitalBase* proxy= new ProxyOrbital(tpq,org);
//#if defined(ENABLE_SMARTPOINTER)
//    return OrbitalBasePtr(proxy);
//#else
//    return proxy;
//#endif
}

OrbitalBase::RealType OrbitalBase::KECorrection()
{
  return 0.0;
}

void OrbitalBase::get_ratios(ParticleSet& P, vector<ValueType>& ratios)
{
  ostringstream o;
  o << "OrbitalBase::get_ratios is not implemented by " << OrbitalName;
  APP_ABORT(o);
}

void OrbitalBase::evaluateRatios(VirtualParticleSet& P, vector<ValueType>& ratios)
{
  ostringstream o;
  o << "OrbitalBase::evaluateRatios is not implemented by " << OrbitalName;
  APP_ABORT(o);
}

void OrbitalBase::evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
    vector<ValueType>& ratios, Matrix<ValueType>& dratios)
{
  //default is only ratios and zero derivatives
  evaluateRatios(VP,ratios);
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

