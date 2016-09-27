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
    
    

#include "QMCWaveFunctions/AFMSPOSet.h"
#include "Numerics/OhmmsBlas.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

//   void
//   AFMSPOSet::addParameter ( std::string id, int iorb, int basis)
//   {
//
//   }

bool
AFMSPOSet::put (xmlNodePtr node, SPOPool_t &spo_pool)
{
  std::string gsName, basisName, opt("true"), varname("theta");
  OhmmsAttributeSet attrib;
  attrib.add (gsName,    "gs_sposet");
  attrib.add (basisName, "basis_sposet");
//     attrib.add (pm, "sign");
  attrib.add (theta, "theta");
  attrib.add (varname, "prefix");
  attrib.add (opt, "optimize");
  attrib.put (node);
  Optimizable = ((opt=="yes")||(opt=="true"));
//     if (N == 0) {
//       app_error() << "You must specify \"size\" attribute for linearopt sposet.\n";
//       abort();
//     }
  /////////////////////////////////////
  // First, find ground-state SPOSet //
  /////////////////////////////////////
  if (gsName == "")
  {
    app_error() << "You must supply \"gs_sposet\".  Aborting.\n";
    abort();
  }
  SPOPool_t::iterator iter = spo_pool.find(gsName);
  if (iter == spo_pool.end())
  {
    app_error() << "No sposet named \"" << gsName << "\" found.  Abort.\n";
    abort();
  }
  else
  {
    app_log() << "  Found ground-state SPOSet \"" << gsName << "\".\n";
    GSOrbitals = iter->second;
  }
  //////////////////////////////////////
  // Now, find basis SPOSet from pool //
  //////////////////////////////////////
  iter = spo_pool.find(basisName);
  if (iter == spo_pool.end())
  {
    app_error() << "No sposet named \"" << basisName << "\" found.  Abort.\n";
    abort();
  }
  else
  {
    BasisOrbitals = iter->second;
    app_log() << "  Found basis SPOSet \"" << basisName << "\".\n";
  }
  N = BasisOrbitals->getOrbitalSetSize();
  OrbitalSetSize = N;
  int M = GSOrbitals->getOrbitalSetSize();
  if (N!=M)
  {
    app_error() << "sposet sizes do not match \"" << N << ",  "<<M<<". ABORT."<< std::endl;
    abort();
  }
  resize(N);
  if(Optimizable)
    myVars.insert(varname,theta,true,optimize::SPO_P);
  resetTheta(theta);
  return SPOSetBase::put(node);
}


void
AFMSPOSet::resetTargetParticleSet(ParticleSet& P)
{
  GSOrbitals->resetTargetParticleSet(P);
  BasisOrbitals->resetTargetParticleSet(P);
}

void
AFMSPOSet::setOrbitalSetSize(int norbs)
{
  OrbitalSetSize = norbs;
}


void
AFMSPOSet::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
}

void
AFMSPOSet::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
}

void
AFMSPOSet::resetParameters(const opt_variables_type& active)
{
  if(Optimizable)
  {
    int loc=myVars.where(0);
    if (loc>=0)
      myVars[0]=active[loc];
    theta=active[loc];
  }
  resetTheta(theta);
}

// Obsolete
void
AFMSPOSet::evaluateDerivatives
(ParticleSet& P, int iat, const opt_variables_type& active,
 ValueMatrix_t& d_phi, ValueMatrix_t& d_lapl_phi)
{
  app_error() << "AFMSPOSet::evaluateDerivatives" << std::endl;
  abort();
}

void
AFMSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  GSOrbitals->evaluate(P,iat,GSVal);
  BasisOrbitals->evaluate(P,iat,BasisVal);
  for (int i=0; i<N; i++)
    psi[i] = costheta*GSVal[i] + pm*sintheta*BasisVal[i];
}

void
AFMSPOSet::evaluate(const ParticleSet& P, const PosType& r,
                    std::vector<RealType> &psi)
{
  app_error() << "AFMSPOSet::evaluate(const ParticleSet& P, const PosType& r, std::vector<RealType> &psi)\n  should not be called.  Abort.\n";
  abort();
}

void
AFMSPOSet::evaluate(const ParticleSet& P, int iat,
                    ValueVector_t& psi, GradVector_t& dpsi,
                    ValueVector_t& d2psi)
{
  GSOrbitals->evaluate(P,iat,GSVal,GSGrad,GSLapl);
  BasisOrbitals->evaluate(P,iat,BasisVal,BasisGrad,BasisLapl);
  for (int iorb=0; iorb<N; iorb++)
  {
    psi  [iorb] = costheta*GSVal[iorb] + pm*sintheta*BasisVal[iorb];
    d2psi[iorb] = costheta*GSLapl[iorb]+ pm*sintheta*BasisLapl[iorb];
    for(int x(0); x<OHMMS_DIM; x++)
      dpsi [iorb][x] = costheta*GSGrad[iorb][x] + pm*sintheta*BasisGrad[iorb][x];
  }
}

void
AFMSPOSet::evaluate(const ParticleSet& P, int iat,
                    ValueVector_t& psi, GradVector_t& dpsi,
                    HessVector_t& gg_psi)
{
  app_error() << "Need specialization for AFMSPOSet::evaluate(HessVector_t)  Abort.\n";
  abort();
}


void AFMSPOSet::evaluateForDeriv (const ParticleSet& P, int first, int last,
                                  ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  GSOrbitals->evaluate_notranspose(P, first, last, GSValMatrix, GSGradMatrix, GSLaplMatrix);
  BasisOrbitals->evaluate_notranspose(P, first, last, BasisValMatrix, BasisGradMatrix, BasisLaplMatrix);
  logdet = dsintheta*BasisValMatrix + dcostheta*GSValMatrix;
  d2logdet = dsintheta*BasisLaplMatrix + dcostheta*GSLaplMatrix;
  // Gradient part.
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        dlogdet(i,j)[dim] = dsintheta*BasisGradMatrix(i,j)[dim] + dcostheta*GSGradMatrix(i,j)[dim];
//     for (int i=0; i<N; i++) for (int j=0; j<N; j++) d2logdet(i,j) -= 2.0*dot(dlogdet(i,j),dlogdet(i,j));
}


void
AFMSPOSet::evaluate_notranspose
(const ParticleSet& P, int first, int last,
 ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  GSOrbitals->evaluate_notranspose(P, first, last, GSValMatrix, GSGradMatrix, GSLaplMatrix);
  BasisOrbitals->evaluate_notranspose(P, first, last, BasisValMatrix, BasisGradMatrix, BasisLaplMatrix);
  logdet = costheta*GSValMatrix +pm*sintheta*BasisValMatrix;
  d2logdet = costheta*GSLaplMatrix +pm*sintheta*BasisLaplMatrix;
  // Gradient part.
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        dlogdet(i,j)[dim] = costheta*GSGradMatrix(i,j)[dim] + pm*sintheta*BasisGradMatrix(i,j)[dim];
}

SPOSetBase*
AFMSPOSet::makeClone() const
{
  SPOSetBase *gs, *basis;
  AFMSPOSet *clone;
  gs = GSOrbitals->makeClone();
  basis = BasisOrbitals->makeClone();
  clone = new AFMSPOSet(N,gs,basis);
  clone->setOrbitalSetSize(N);
  clone->Optimizable=Optimizable;
  clone->myVars=myVars;
  clone->resetTheta(theta);
  clone->pm=pm;
  return clone;
}


}
