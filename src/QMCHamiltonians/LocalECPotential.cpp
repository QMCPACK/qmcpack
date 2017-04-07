//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

LocalECPotential::LocalECPotential(const ParticleSet& ions, ParticleSet& els):
  IonConfig(ions),Peln(els),Pion(ions)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions,els);
  NumIons=ions.getTotalNum();
  myTableIndex=els.addTable(ions,DT_SOA_PREFERRED);
  //allocate null
  PPset.resize(ions.getSpeciesSet().getTotalNum(),0);
  PP.resize(NumIons,nullptr);
  Zeff.resize(NumIons,0.0);
  gZeff.resize(ions.getSpeciesSet().getTotalNum(),0);
}

///destructor
LocalECPotential::~LocalECPotential()
{
  delete_iter(PPset.begin(),PPset.end());
  //map<int,RadialPotentialType*>::iterator pit(PPset.begin()), pit_end(PPset.end());
  //while(pit != pit_end) {
  //  delete (*pit).second; ++pit;
  //}
}

void LocalECPotential::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(IonConfig,DT_SOA_PREFERRED);
  if(tid != myTableIndex)
  {
    APP_ABORT("  LocalECPotential::resetTargetParticleSet found a different distance table index.");
  }
}

void LocalECPotential::add(int groupID, RadialPotentialType* ppot, RealType z)
{
  PPset[groupID]=ppot;
  gZeff[groupID]=z;
  for(int iat=0; iat<PP.size(); iat++)
  {
    if(IonConfig.GroupID[iat]==groupID)
    {
      PP[iat]=ppot;
      Zeff[iat]=z;
    }
  }
}

#if !defined(REMOVE_TRACEMANAGER)
void LocalECPotential::contribute_particle_quantities()
{
  request.contribute_array(myName);
}

void LocalECPotential::checkout_particle_quantities(TraceManager& tm)
{
  streaming_particles = request.streaming_array(myName);
  if( streaming_particles)
  {
    Ve_sample = tm.checkout_real<1>(myName,Peln);
    Vi_sample = tm.checkout_real<1>(myName,Pion);
  }
}

void LocalECPotential::delete_particle_quantities()
{
  if( streaming_particles)
  {
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif


LocalECPotential::Return_t
LocalECPotential::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if( streaming_particles)
    Value = evaluate_sp(P);
  else
#endif
  {
    const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
    Value=0.0;
    if(d_table.DTType==DT_SOA)
    {
      for(int iat=0; iat<NumIons; ++iat)
      {
        RadialPotentialType* ppot(PP[iat]);
        if(ppot==nullptr) continue;
        RealType esum(0);
        const RealType* restrict dist=d_table.r_m2[iat];
        for(int nj=0, jmax=d_table.nadj(iat); nj<jmax; ++nj)
        {
          esum += ppot->splint(dist[nj])/dist[nj];
        }
        Value -= esum*Zeff[iat];
      }
    }
    else
    {
      //loop over all the ions
      for(int iat=0; iat<NumIons; iat++)
      {
        RadialPotentialType* ppot(PP[iat]);
        if(ppot==nullptr) continue;
        RealType esum(0);
        for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
          esum += ppot->splint(d_table.r(nn))*d_table.rinv(nn);
        //count the sign and effective charge
        Value -= esum*Zeff[iat];
      }
    }
  }
  return Value;
}



#if !defined(REMOVE_TRACEMANAGER)
LocalECPotential::Return_t
LocalECPotential::evaluate_sp(ParticleSet& P)
{
  const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
  Value=0.0;
  Array<RealType,1>& Ve_samp = *Ve_sample;
  Array<RealType,1>& Vi_samp = *Vi_sample;
  Ve_samp = 0.0;
  Vi_samp = 0.0;
  //loop over all the ions
  for(int iat=0; iat<NumIons; iat++)
  {
    RadialPotentialType* ppot(PP[iat]);
    if(ppot==0)
      continue;
    Return_t esum(0.0),pairpot;
    //loop over all the electrons
    for(int nn=d_table.M[iat],iel=0; nn<d_table.M[iat+1]; ++nn,iel++)
    {
      pairpot = -.5*Zeff[iat]*ppot->splint(d_table.r(nn))*d_table.rinv(nn);
      Vi_samp(iat) += pairpot;
      Ve_samp(iel) += pairpot;
      esum += pairpot;
    }
    Value += esum;
  }
  Value *= 2.0;
#if defined(TRACE_CHECK)
  RealType Vnow  = Value;
  RealType Visum = Vi_samp.sum();
  RealType Vesum = Ve_samp.sum();
  RealType Vsum  = Vesum+Visum;
  RealType Vorig = evaluate_orig(P);
  if(std::abs(Vsum-Vnow)>TraceManager::trace_tol)
  {
    app_log()<<"accumtest: LocalECPotential::evaluate()"<< std::endl;
    app_log()<<"accumtest:   tot:"<< Vnow << std::endl;
    app_log()<<"accumtest:   sum:"<< Vsum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vesum-Visum)>TraceManager::trace_tol)
  {
    app_log()<<"sharetest: LocalECPotential::evaluate()"<< std::endl;
    app_log()<<"sharetest:   e share:"<< Vesum << std::endl;
    app_log()<<"sharetest:   i share:"<< Visum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if(std::abs(Vorig-Vnow)>TraceManager::trace_tol)
  {
    app_log()<<"versiontest: LocalECPotential::evaluate()"<< std::endl;
    app_log()<<"versiontest:   orig:"<< Vorig << std::endl;
    app_log()<<"versiontest:    mod:"<< Vnow << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return Value;
}
#endif



LocalECPotential::Return_t
LocalECPotential::evaluate_orig(ParticleSet& P)
{
  const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
  Value=0.0;
  //loop over all the ions
  for(int iat=0; iat<NumIons; iat++)
  {
    RadialPotentialType* ppot(PP[iat]);
    if(ppot==0)
      continue;
    Return_t esum(0.0);
    for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
      esum += ppot->splint(d_table.r(nn))*d_table.rinv(nn);
    //count the sign and effective charge
    Value -= esum*Zeff[iat];
  }
  return Value;
}


LocalECPotential::Return_t
LocalECPotential::registerData(ParticleSet& P, BufferType& buffer)
{
  PPart.resize(P.getTotalNum());
  NewValue=Value=evaluateForPbyP(P);
  buffer.add(PPart.begin(),PPart.end());
  buffer.add(Value);
  return Value;
}

LocalECPotential::Return_t
LocalECPotential::updateBuffer(ParticleSet& P, BufferType& buffer)
{
  NewValue=Value=evaluateForPbyP(P);
  buffer.put(PPart.begin(),PPart.end());
  buffer.put(Value);
  return Value;
}

void LocalECPotential::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.get(PPart.begin(),PPart.end());
  buffer.get(Value);
}

void LocalECPotential::copyToBuffer(ParticleSet& P, BufferType& buffer)
{
  buffer.put(PPart.begin(),PPart.end());
  buffer.put(Value);
}

LocalECPotential::Return_t
LocalECPotential::evaluateForPbyP(ParticleSet& P)
{
  const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
  PPart=0.0;
  Return_t res=0.0;
  for(int iat=0; iat<NumIons; ++iat)
  {
    RadialPotentialType* ppot(PP[iat]);
    if(ppot)
    {
      Return_t z=-Zeff[iat];
      for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
      {
        Return_t e= z*ppot->splint(d_table.r(nn))*d_table.rinv(nn);
        PPart[d_table.J[nn]]+=e;
        res+=e;
      }
    }
  }
  return res;
}

LocalECPotential::Return_t
LocalECPotential::evaluatePbyP(ParticleSet& P, int active)
{
  const std::vector<DistanceTableData::TempDistType> &temp(P.DistTables[myTableIndex]->Temp);
  PPtmp=0.0;
  for(int iat=0; iat<NumIons; ++iat)
  {
    if(PP[iat])
      PPtmp -= Zeff[iat]*PP[iat]->splint(temp[iat].r1)*temp[iat].rinv1;
  }
  return NewValue=Value+PPtmp-PPart[active];
}

void LocalECPotential::acceptMove(int active)
{
  PPart[active]=PPtmp;
  Value=NewValue;
}

QMCHamiltonianBase* LocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  LocalECPotential* myclone=new LocalECPotential(IonConfig,qp);
  for(int ig=0; ig<PPset.size(); ++ig)
  {
    if(PPset[ig])
    {
      RadialPotentialType* ppot=PPset[ig]->makeClone();
      myclone->add(ig,ppot,gZeff[ig]);
    }
  }
  return myclone;
}
}

