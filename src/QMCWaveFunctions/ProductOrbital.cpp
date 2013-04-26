//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
/** @file ProductOrbital.cpp
 * @brief Definitions of ProductOrbital
 */
#include "QMCWaveFunctions/ProductOrbital.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

ProductOrbital::~ProductOrbital()
{
  delete_iter(Psi.begin(), Psi.end());
  delete Constraints;
}

void ProductOrbital::resetParameters(const opt_variables_type& active)
{
  //APP_ABORT("ProductOrbital::resetParameters is incomplete");
  //Constraints->resetParameters(active);
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->resetParameters(active);
}

void ProductOrbital::checkOutVariables(const opt_variables_type& o)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->checkOutVariables(o);
}

void ProductOrbital::checkInVariables(opt_variables_type& o)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->checkInVariables(o);
}

void ProductOrbital::reportStatus(ostream& os)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->reportStatus(os);
}

void ProductOrbital::resetTargetParticleSet(ParticleSet& P)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->resetTargetParticleSet(P);
}

ProductOrbital::RealType ProductOrbital::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->evaluateLog(P,G,L);
  return LogValue;
}

ProductOrbital::ValueType ProductOrbital::ratio(ParticleSet& P, int iat
    , ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL)
{
  ValueType r(1.0);
  for(int i=0; i<Psi.size(); i++)
    r *= Psi[i]->ratio(P,iat,dG,dL);
  return r;
}

ProductOrbital::ValueType ProductOrbital::ratio(ParticleSet& P, int iat)
{
  ValueType r(1.0);
  for(int i=0; i<Psi.size(); i++)
    r *= Psi[i]->ratio(P,iat);
  return r;
}

//ProductOrbital::ValueType
//  ProductOrbital::logRatio(ParticleSet& P, int iat,
//      ParticleSet::ParticleGradient_t& dG,
//      ParticleSet::ParticleLaplacian_t& dL) {
//    ValueType r(0.0);
//    for(int i=0; i<Psi.size(); i++)
//      r += Psi[i]->logRatio(P,iat,dG,dL);
//    return r;
//  }

void ProductOrbital::acceptMove(ParticleSet& P, int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->acceptMove(P,iat);
}

void ProductOrbital::restore(int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->restore(iat);
}

void ProductOrbital::update(ParticleSet& P
                            , ParticleSet::ParticleGradient_t& dG , ParticleSet::ParticleLaplacian_t& dL, int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->update(P,dG,dL,iat);
}

ProductOrbital::RealType ProductOrbital::registerData(ParticleSet& P, BufferType& buf)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->registerData(P,buf);
  return LogValue;
}

ProductOrbital::RealType ProductOrbital::updateBuffer(ParticleSet& P
    , BufferType& buf, bool fromscratch)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->updateBuffer(P,buf,fromscratch);
  return LogValue;
}

void ProductOrbital::copyFromBuffer(ParticleSet& P, BufferType& buf)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->copyFromBuffer(P,buf);
}

ProductOrbital::RealType ProductOrbital::evaluateLog(ParticleSet& P,BufferType& buf)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->evaluateLog(P,buf);
  return LogValue;
}

OrbitalBase* ProductOrbital::makeClone(ParticleSet& tpq) const
{
  app_warning() << "  ProductOrbital::makeClone for long-range breakup stuff won't work." << endl;
  ProductOrbital* myclone=new ProductOrbital(*this);
  for(int i=0; i<Psi.size(); ++i)
  {
    myclone->Psi[i]=Psi[i]->makeClone(tpq);
  }
  return myclone;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
