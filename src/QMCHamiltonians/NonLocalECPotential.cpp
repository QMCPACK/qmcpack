//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

void NonLocalECPotential::resetTargetParticleSet(ParticleSet& P)
{
}

/** constructor
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\param psi trial wavefunction
 */
NonLocalECPotential::NonLocalECPotential(ParticleSet& ions, ParticleSet& els,
    TrialWaveFunction& psi, bool computeForces, bool useVP):
  IonConfig(ions), Psi(psi), UseTMove(TMOVE_OFF), myRNG(&Random),
  nonLocalOps(els.getTotalNum()), ComputeForces(computeForces),
  UseVP(useVP), ForceBase(ions,els), Peln(els)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions,els);
  myTableIndex=els.addTable(ions,DT_SOA_PREFERRED);
  NumIons=ions.getTotalNum();
  //els.resizeSphere(NumIons);
  PP.resize(NumIons,nullptr);
  prefix="FNL";
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum(),0);
  PulayTerm.resize(NumIons);
  
  UpdateMode.set(NONLOCAL,1);
}

///destructor
NonLocalECPotential::~NonLocalECPotential()
{
  delete_iter(PPset.begin(),PPset.end());
  //map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
  //while(pit != pit_end) {
  //   delete (*pit).second; ++pit;
  //}
}


#if !defined(REMOVE_TRACEMANAGER)
void NonLocalECPotential::contribute_particle_quantities()
{
  request.contribute_array(myName);
}

void NonLocalECPotential::checkout_particle_quantities(TraceManager& tm)
{
  streaming_particles = request.streaming_array(myName);
  if( streaming_particles)
  {
    Ve_sample = tm.checkout_real<1>(myName,Peln);
    Vi_sample = tm.checkout_real<1>(myName,IonConfig);
    for(int iat=0; iat<NumIons; iat++)
    {
      if(PP[iat])
      {
      PP[iat]->streaming_particles = streaming_particles;
      PP[iat]->Ve_sample = Ve_sample;
      PP[iat]->Vi_sample = Vi_sample;
      }
    }
  }
}

void NonLocalECPotential::delete_particle_quantities()
{
  if( streaming_particles)
  {
    for(int iat=0; iat<NumIons; iat++)
    {
      if(PP[iat])
      {
        PP[iat]->streaming_particles = false;
        PP[iat]->Ve_sample = NULL;
        PP[iat]->Vi_sample = NULL;
      }
    }
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif

NonLocalECPotential::Return_t
NonLocalECPotential::evaluate(ParticleSet& P)
{
  evaluate(P, false);
  return Value;
}

NonLocalECPotential::Return_t
NonLocalECPotential::evaluateWithToperator(ParticleSet& P)
{
  if( UseTMove==TMOVE_V0 || UseTMove==TMOVE_V3 )
  {
    nonLocalOps.reset();
    evaluate(P, true);
  }
  else
    evaluate(P, false);
  return Value;
}

void
NonLocalECPotential::evaluate(ParticleSet& P, bool Tmove)
{
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);
  Value=0.0;
#if !defined(REMOVE_TRACEMANAGER)
  if( streaming_particles)
  {
    (*Ve_sample) = 0.0;
    (*Vi_sample) = 0.0;
  }
#endif
  for(int ipp=0; ipp<PPset.size(); ipp++)
    if(PPset[ipp]) PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  if (ComputeForces)
  {
    for(int iat=0; iat<NumIons; iat++)
      if(PP[iat])
      {
        if(Tmove)
          Value += PP[iat]->evaluate(P,Psi,iat, Txy, forces[iat]);
        else
          Value += PP[iat]->evaluate(P,IonConfig,iat,Psi, forces[iat], PulayTerm[iat]);
      }
  }
  else
  {
    const auto myTable = P.DistTables[myTableIndex];
    if(myTable->DTType == DT_SOA)
    {
      for(int jel=0; jel<P.getTotalNum(); jel++)
      {
        const auto &dist  = myTable->Distances[jel];
        const auto &displ = myTable->Displacements[jel];
        for(int iat=0; iat<NumIons; iat++)
          if(PP[iat]!=nullptr && dist[iat]<PP[iat]->Rmax)
            Value += PP[iat]->evaluateOne(P,iat,Psi,jel,dist[iat],RealType(-1)*displ[iat],Tmove,Txy);
      }
    }
    else
    {
      for(int iat=0; iat<NumIons; iat++)
      {
        if(PP[iat]==nullptr) continue;
        for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++)
        {
          const RealType r(myTable->r(nn));
          if(r>PP[iat]->Rmax) continue;
          Value += PP[iat]->evaluateOne(P,iat,Psi,iel,r,myTable->dr(nn),Tmove,Txy);
        }
      }
    }
  }
#if defined(TRACE_CHECK)
  if( streaming_particles)
  {
    Return_t Vnow = Value;
    RealType Visum = Vi_sample->sum();
    RealType Vesum = Ve_sample->sum();
    RealType Vsum  = Vesum+Visum;
    if(std::abs(Vsum-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"accumtest: NonLocalECPotential::evaluate()"<< std::endl;
      app_log()<<"accumtest:   tot:"<< Vnow << std::endl;
      app_log()<<"accumtest:   sum:"<< Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    if(std::abs(Vesum-Visum)>TraceManager::trace_tol)
    {
      app_log()<<"sharetest: NonLocalECPotential::evaluate()"<< std::endl;
      app_log()<<"sharetest:   e share:"<< Vesum << std::endl;
      app_log()<<"sharetest:   i share:"<< Visum << std::endl;
      APP_ABORT("Trace check failed");
    }
  }
#endif
}

void
NonLocalECPotential::computeOneElectronTxy(ParticleSet& P, const int ref_elec)
{
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);
  const auto myTable = P.DistTables[myTableIndex];
  if(myTable->DTType == DT_SOA)
  {
    const auto &dist  = myTable->Distances[ref_elec];
    const auto &displ = myTable->Displacements[ref_elec];
    for(int iat=0; iat<NumIons; iat++)
      if(PP[iat]!=nullptr && dist[iat]<PP[iat]->Rmax)
        PP[iat]->evaluateOne(P,iat,Psi,ref_elec,dist[iat],RealType(-1)*displ[iat],true,Txy);
  }
  else
  {
    for(int iat=0; iat<NumIons; iat++)
    {
      if(PP[iat]==nullptr) continue;
      for(int nn=myTable->M[iat],iel=0; nn<myTable->M[iat+1]; nn++,iel++)
      {
        const RealType r(myTable->r(nn));
        if(r<PP[iat]->Rmax && iel==ref_elec)
          PP[iat]->evaluateOne(P,iat,Psi,iel,r,myTable->dr(nn),true,Txy);
      }
    }
  }
}

int
NonLocalECPotential::makeNonLocalMovesPbyP(ParticleSet& P)
{
  int NonLocalMoveAccepted = 0;
  RandomGenerator_t& RandomGen(*myRNG);
  if(UseTMove==TMOVE_V0)
  {
    int ibar = nonLocalOps.selectMove(RandomGen());
    //make a non-local move
    if(ibar)
    {
      int iat=nonLocalOps.id(ibar);
      P.setActive(iat);
      if(P.makeMoveAndCheck(iat,nonLocalOps.delta(ibar)))
      {
        GradType grad_iat;
        Psi.ratioGrad(P,iat,grad_iat);
        Psi.acceptMove(P,iat);
        P.acceptMove(iat);
        NonLocalMoveAccepted++;
      }
    }
  }
  else if(UseTMove==TMOVE_V1)
  {
    GradType grad_iat;
    //make a non-local move per particle
    for(int ig=0; ig<P.groups(); ++ig) //loop over species
    {
      for (int iat=P.first(ig); iat<P.last(ig); ++iat)
      {
        nonLocalOps.reset();
        computeOneElectronTxy(P,iat);
        int ibar = nonLocalOps.selectMove(RandomGen());
        if(ibar)
        {
          P.setActive(iat);
          if(P.makeMoveAndCheck(iat,nonLocalOps.delta(ibar)))
          {
            Psi.ratioGrad(P,iat,grad_iat);
            Psi.acceptMove(P,iat);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
    }
  }
  else if(UseTMove==TMOVE_V3)
  {
    nonLocalOps.group_by_elec();
    GradType grad_iat;
    size_t NonLocalMoveAcceptedTemp = 0;
    //make a non-local move per particle
    for(int ig=0; ig<P.groups(); ++ig) //loop over species
    {
      for (int iat=P.first(ig); iat<P.last(ig); ++iat)
      {
        const NonLocalData *oneTMove = nonLocalOps.selectMove(RandomGen(), iat);
        if(oneTMove)
        {
          P.setActive(iat);
          if(P.makeMoveAndCheck(iat,oneTMove->Delta))
          {
            Psi.ratioGrad(P,iat,grad_iat);
            Psi.acceptMove(P,iat);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
    }
  }
  return NonLocalMoveAccepted;
}

void
NonLocalECPotential::add(int groupID, NonLocalECPComponent* ppot)
{
  //map<int,NonLocalECPComponent*>::iterator pit(PPset.find(groupID));
  //ppot->myTable=d_table;
  //if(pit  == PPset.end()) {
  //  for(int iat=0; iat<PP.size(); iat++) {
  //    if(IonConfig.GroupID[iat]==groupID) PP[iat]=ppot;
  //  }
  //  PPset[groupID]=ppot;
  //}
  ppot->myTableIndex=myTableIndex;
  for(int iat=0; iat<PP.size(); iat++)
    if(IonConfig.GroupID[iat]==groupID)
      PP[iat]=ppot;
  PPset[groupID]=ppot;
  if(UseVP && ppot->VP==0) ppot->initVirtualParticle(Peln);
}

QMCHamiltonianBase* NonLocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  NonLocalECPotential* myclone=new NonLocalECPotential(IonConfig,qp,psi,
      ComputeForces, UseVP);
  for(int ig=0; ig<PPset.size(); ++ig)
  {
    if(PPset[ig])
    {
      NonLocalECPComponent* ppot=PPset[ig]->makeClone(qp);
      myclone->add(ig,ppot);
    }
  }
  return myclone;
}


void NonLocalECPotential::addObservables(PropertySetType& plist
    , BufferType& collectables)
{
  QMCHamiltonianBase::addValue(plist);
  if (ComputeForces)
  {
    if(FirstForceIndex<0)
      FirstForceIndex=plist.size();
    for(int iat=0; iat<Nnuc; iat++)
    {
      for(int x=0; x<OHMMS_DIM; x++)
      {
        std::ostringstream obsName1, obsName2;
        obsName1 << "FNL" << "_" << iat << "_" << x;
        plist.add(obsName1.str());
        obsName2 << "FNL_Pulay" << "_" << iat << "_" << x;
        plist.add(obsName2.str());
      }
    }
  }
}

void
NonLocalECPotential::registerObservables(std::vector<observable_helper*>& h5list,
    hid_t gid) const
{
  QMCHamiltonianBase::registerObservables(h5list, gid);
  if (ComputeForces)
  {
    std::vector<int> ndim(2);
    ndim[0]=Nnuc;
    ndim[1]=OHMMS_DIM;
    observable_helper* h5o1 = new observable_helper("FNL");
    h5o1->set_dimensions(ndim,FirstForceIndex);
    h5o1->open(gid);
    h5list.push_back(h5o1);
    observable_helper* h5o2 = new observable_helper("FNL_Pulay");
    h5o2->set_dimensions(ndim,FirstForceIndex+Nnuc*OHMMS_DIM);
    h5o2->open(gid);
    h5list.push_back(h5o2);
  }
}

void
NonLocalECPotential::setObservables(QMCTraits::PropertySetType& plist)
{
  QMCHamiltonianBase::setObservables(plist);
  if (ComputeForces)
  {
    int index = FirstForceIndex;
    for(int iat=0; iat<Nnuc; iat++)
    {
      for(int x=0; x<OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


void
NonLocalECPotential::setParticlePropertyList(QMCTraits::PropertySetType& plist,
    int offset)
{
  QMCHamiltonianBase::setParticlePropertyList (plist, offset);
  if (ComputeForces)
  {
    int index = FirstForceIndex + offset;
    for(int iat=0; iat<Nnuc; iat++)
    {
      for(int x=0; x<OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


}
