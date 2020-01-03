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
#include "QMCHamiltonians/NonLocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"

#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
void NonLocalECPotential::resetTargetParticleSet(ParticleSet& P) {}

/** constructor
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\param psi trial wavefunction
 */
NonLocalECPotential::NonLocalECPotential(ParticleSet& ions,
                                         ParticleSet& els,
                                         TrialWaveFunction& psi,
                                         bool computeForces,
                                         bool useVP)
    : ForceBase(ions, els),
      myRNG(nullptr),
      IonConfig(ions),
      Psi(psi),
      Peln(els),
      ElecNeighborIons(els),
      IonNeighborElecs(ions),
      UseTMove(TMOVE_OFF),
      nonLocalOps(els.getTotalNum()),
      ComputeForces(computeForces)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions, els);
  myTableIndex = els.addTable(ions, DT_SOA_PREFERRED);
  NumIons      = ions.getTotalNum();
  //els.resizeSphere(NumIons);
  PP.resize(NumIons, nullptr);
  prefix = "FNL";
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum(), 0);
  PulayTerm.resize(NumIons);
  UpdateMode.set(NONLOCAL, 1);
  Ve_samp_tmp.resize(els.getTotalNum());
  Vi_samp_tmp.resize(ions.getTotalNum());
}

///destructor
NonLocalECPotential::~NonLocalECPotential()
{
  delete_iter(PPset.begin(), PPset.end());
  //map<int,NonLocalECPComponent*>::iterator pit(PPset.begin()), pit_end(PPset.end());
  //while(pit != pit_end) {
  //   delete (*pit).second; ++pit;
  //}
}


#if !defined(REMOVE_TRACEMANAGER)
void NonLocalECPotential::contribute_particle_quantities() { request.contribute_array(myName); }

void NonLocalECPotential::checkout_particle_quantities(TraceManager& tm)
{
  streaming_particles = request.streaming_array(myName);
  if (streaming_particles)
  {
    Ve_sample = tm.checkout_real<1>(myName, Peln);
    Vi_sample = tm.checkout_real<1>(myName, IonConfig);
    for (int iat = 0; iat < NumIons; iat++)
    {
      if (PP[iat])
      {
        PP[iat]->streaming_particles = streaming_particles;
        //PP[iat]->Ve_sample           = Ve_sample;
        //PP[iat]->Vi_sample           = Vi_sample;
      }
    }
  }
}

void NonLocalECPotential::delete_particle_quantities()
{
  if (streaming_particles)
  {
    for (int iat = 0; iat < NumIons; iat++)
    {
      if (PP[iat])
      {
        PP[iat]->streaming_particles = false;
        //PP[iat]->Ve_sample           = NULL;
        //PP[iat]->Vi_sample           = NULL;
      }
    }
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif

NonLocalECPotential::Return_t NonLocalECPotential::evaluate(ParticleSet& P)
{
  evaluateImpl(P, false);
  return Value;
}

void NonLocalECPotential::mw_evaluate(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list)
{
  mw_evaluateImpl(O_list, P_list, false);
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateWithToperator(ParticleSet& P)
{
  if (UseTMove == TMOVE_V0 || UseTMove == TMOVE_V3)
    evaluateImpl(P, true);
  else
    evaluateImpl(P, false);
  return Value;
}

void NonLocalECPotential::mw_evaluateWithToperator(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list)
{
  if (UseTMove == TMOVE_V0 || UseTMove == TMOVE_V3)
    mw_evaluateImpl(O_list, P_list, true);
  else
    mw_evaluateImpl(O_list, P_list, false);
}

void NonLocalECPotential::evaluateImpl(ParticleSet& P, bool Tmove)
{
  if (Tmove)
    nonLocalOps.reset();
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);
  Value = 0.0;
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles)
  {
    (*Ve_sample) = 0.0;
    (*Vi_sample) = 0.0;
  }
  auto& Ve_samp = Ve_samp_tmp;
  auto& Vi_samp = Vi_samp_tmp;
  Ve_samp = 0.0;
  Vi_samp = 0.0;
#endif
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  const auto& myTable = P.getDistTable(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  if (ComputeForces)
  {
    forces = 0;
    if (myTable.DTType == DT_SOA)
    {
      for (int jel = 0; jel < P.getTotalNum(); jel++)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < NumIons; iat++)
          if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
          {
            RealType pairpot = PP[iat]->evaluateOneWithForces(P, iat, Psi, jel, dist[iat], -displ[iat], forces[iat]);
            if (Tmove) PP[iat]->contributeTxy(jel, Txy);
            Value += pairpot;
            NeighborIons.push_back(iat);
            IonNeighborElecs.getNeighborList(iat).push_back(jel);
          }
      }
    }
    else
    {
      APP_ABORT("NonLocalECPotential::evaluate():  Forces not imlpemented for AoS build\n");
    }
  }
  else
  {
    if (myTable.DTType == DT_SOA)
    {
      for (int jel = 0; jel < P.getTotalNum(); jel++)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < NumIons; iat++)
          if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
          {
            RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, jel, dist[iat], -displ[iat]);
            if (Tmove) PP[iat]->contributeTxy(jel, Txy);
            Value += pairpot;
            NeighborIons.push_back(iat);
            IonNeighborElecs.getNeighborList(iat).push_back(jel);
            Ve_samp(jel) = 0.5*pairpot;
            Vi_samp(iat) = 0.5*pairpot;
          }
      }
    }
    else
    {
#ifndef ENABLE_SOA
      for (int iat = 0; iat < NumIons; iat++)
      {
        if (PP[iat] == nullptr)
          continue;
        std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
        for (int nn = myTable.M[iat], iel = 0; nn < myTable.M[iat + 1]; nn++, iel++)
        {
          const RealType r(myTable.r(nn));
          if (r > PP[iat]->getRmax())
            continue;
          RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, iel, r, myTable.dr(nn));
          if (Tmove) PP[iat]->contributeTxy(iel, Txy);
          Value += pairpot;
          NeighborElecs.push_back(iel);
          ElecNeighborIons.getNeighborList(iel).push_back(iat);
          Ve_samp(iel) = 0.5*pairpot;
          Vi_samp(iat) = 0.5*pairpot;
        }
      }
#endif
    }
  }

#if defined(TRACE_CHECK)
  if (streaming_particles)
  {
    Return_t Vnow  = Value;
    RealType Visum = Vi_sample->sum();
    RealType Vesum = Ve_sample->sum();
    RealType Vsum  = Vesum + Visum;
    if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "accumtest: NonLocalECPotential::evaluate()" << std::endl;
      app_log() << "accumtest:   tot:" << Vnow << std::endl;
      app_log() << "accumtest:   sum:" << Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    if (std::abs(Vesum - Visum) > TraceManager::trace_tol)
    {
      app_log() << "sharetest: NonLocalECPotential::evaluate()" << std::endl;
      app_log() << "sharetest:   e share:" << Vesum << std::endl;
      app_log() << "sharetest:   i share:" << Visum << std::endl;
      APP_ABORT("Trace check failed");
    }
  }
#endif
}

void NonLocalECPotential::mw_evaluateImpl(const RefVector<OperatorBase>& O_list, const RefVector<ParticleSet>& P_list, bool Tmove)
{
  if (Tmove)
    nonLocalOps.reset();

  auto& P = P_list[0].get();
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);

  Value = 0.0;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->randomize_grid(*myRNG);

  //loop over all the ions
  const auto& myTable = P.getDistTable(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  if (ComputeForces)
  {
    APP_ABORT("NonLocalECPotential::mw_evaluateImpl(): Forces not imlpemented\n");
  }
  else
  {
    if (myTable.DTType == DT_SOA)
    {
      for (int jel = 0; jel < P.getTotalNum(); jel++)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < NumIons; iat++)
          if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
          {
            RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, jel, dist[iat], -displ[iat]);
            if (Tmove) PP[iat]->contributeTxy(jel, Txy);
            Value += pairpot;
            NeighborIons.push_back(iat);
            IonNeighborElecs.getNeighborList(iat).push_back(jel);
          }
      }
    }
    else
    {
      APP_ABORT("NonLocalECPotential::mw_evaluateImpl(): not imlpemented for AoS builds\n");
    }
  }
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateWithIonDerivs(ParticleSet& P,
                                                                         ParticleSet& ions,
                                                                         TrialWaveFunction& psi,
                                                                         ParticleSet::ParticlePos_t& hf_terms,
                                                                         ParticleSet::ParticlePos_t& pulay_terms)
{
  //We're going to ignore psi and use the internal Psi.
  //
  //Dummy vector for now.  Tmoves not implemented
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);
  bool Tmove = false;

  forces    = 0;
  PulayTerm = 0;

  Value = 0.0;

  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  const auto& myTable = P.getDistTable(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  if (myTable.DTType == DT_SOA)
  {
    for (int jel = 0; jel < P.getTotalNum(); jel++)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          Value += PP[iat]->evaluateOneWithForces(P, ions, iat, Psi, jel, dist[iat], -displ[iat],
                                                  forces[iat], PulayTerm);
          if (Tmove) PP[iat]->contributeTxy(jel, Txy);
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }
  else
  {
    APP_ABORT("NonLocalECPotential::evaluate():  Forces not imlpemented for AoS build\n");
  }
  hf_terms    -= forces;
  pulay_terms -= PulayTerm;
  return Value;
}

void NonLocalECPotential::computeOneElectronTxy(ParticleSet& P, const int ref_elec)
{
  nonLocalOps.reset();
  std::vector<NonLocalData>& Txy(nonLocalOps.Txy);
  const auto& myTable                  = P.getDistTable(myTableIndex);
  const std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(ref_elec);

  if (myTable.DTType == DT_SOA)
  {
    const auto& dist  = myTable.getDistRow(ref_elec);
    const auto& displ = myTable.getDisplRow(ref_elec);
    for (int atom_index = 0; atom_index < NeighborIons.size(); atom_index++)
    {
      const int iat = NeighborIons[atom_index];
      PP[iat]->evaluateOne(P, iat, Psi, ref_elec, dist[iat], -displ[iat]);
      PP[iat]->contributeTxy(ref_elec, Txy);
    }
  }
  else
  {
#ifndef ENABLE_SOA
    for (int atom_index = 0; atom_index < NeighborIons.size(); atom_index++)
    {
      const int iat = NeighborIons[atom_index];
      int nn        = myTable.M[iat] + ref_elec;
      PP[iat]->evaluateOne(P, iat, Psi, ref_elec, myTable.r(nn), myTable.dr(nn));
      PP[iat]->contributeTxy(ref_elec, Txy);
    }
#endif
  }
}

int NonLocalECPotential::makeNonLocalMovesPbyP(ParticleSet& P)
{
  int NonLocalMoveAccepted = 0;
  RandomGenerator_t& RandomGen(*myRNG);
  if (UseTMove == TMOVE_V0)
  {
    const NonLocalData* oneTMove = nonLocalOps.selectMove(RandomGen());
    //make a non-local move
    if (oneTMove)
    {
      int iat = oneTMove->PID;
      if (P.makeMoveAndCheck(iat, oneTMove->Delta))
      {
        GradType grad_iat;
        Psi.calcRatioGrad(P, iat, grad_iat);
        Psi.acceptMove(P, iat);
        P.acceptMove(iat);
        NonLocalMoveAccepted++;
      }
    }
  }
  else if (UseTMove == TMOVE_V1)
  {
    GradType grad_iat;
    //make a non-local move per particle
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
      for (int iat = P.first(ig); iat < P.last(ig); ++iat)
      {
        computeOneElectronTxy(P, iat);
        const NonLocalData* oneTMove = nonLocalOps.selectMove(RandomGen());
        if (oneTMove)
        {
          if (P.makeMoveAndCheck(iat, oneTMove->Delta))
          {
            Psi.calcRatioGrad(P, iat, grad_iat);
            Psi.acceptMove(P, iat);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
  }
  else if (UseTMove == TMOVE_V3)
  {
    elecTMAffected.assign(P.getTotalNum(), false);
    nonLocalOps.group_by_elec();
    GradType grad_iat;
    //make a non-local move per particle
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
      for (int iat = P.first(ig); iat < P.last(ig); ++iat)
      {
        const NonLocalData* oneTMove;
        if (elecTMAffected[iat])
        {
          // recompute Txy for the given electron effected by Tmoves
          computeOneElectronTxy(P, iat);
          oneTMove = nonLocalOps.selectMove(RandomGen());
        }
        else
          oneTMove = nonLocalOps.selectMove(RandomGen(), iat);
        if (oneTMove)
        {
          if (P.makeMoveAndCheck(iat, oneTMove->Delta))
          {
            Psi.calcRatioGrad(P, iat, grad_iat);
            Psi.acceptMove(P, iat);
            // mark all affected electrons
            markAffectedElecs(P.getDistTable(myTableIndex), iat);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
  }

  if (NonLocalMoveAccepted > 0)
    Psi.completeUpdates();

  return NonLocalMoveAccepted;
}

void NonLocalECPotential::markAffectedElecs(const DistanceTableData& myTable, int iel)
{
  std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(iel);
  for (int iat = 0; iat < NumIons; iat++)
  {
    if (PP[iat] == nullptr)
      continue;
    RealType old_distance = 0.0;
    RealType new_distance = 0.0;
    if (myTable.DTType == DT_SOA)
    {
      old_distance = myTable.getDistRow(iel)[iat];
      new_distance = myTable.getTempDists()[iat];
    }
    else
    {
#ifndef ENABLE_SOA
      old_distance = myTable.r(myTable.M[iat] + iel);
      new_distance = myTable.Temp[iat].r1;
#endif
    }
    bool moved = false;
    // move out
    if (old_distance < PP[iat]->getRmax() && new_distance >= PP[iat]->getRmax())
    {
      moved                           = true;
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      auto iter_at                    = std::find(NeighborIons.begin(), NeighborIons.end(), iat);
      auto iter_el                    = std::find(NeighborElecs.begin(), NeighborElecs.end(), iel);
      *iter_at                        = NeighborIons.back();
      *iter_el                        = NeighborElecs.back();
      NeighborIons.pop_back();
      NeighborElecs.pop_back();
      elecTMAffected[iel] = true;
    }
    // move in
    if (old_distance >= PP[iat]->getRmax() && new_distance < PP[iat]->getRmax())
    {
      moved                           = true;
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      NeighborElecs.push_back(iel);
      NeighborIons.push_back(iat);
    }
    // move around
    if (moved || (old_distance < PP[iat]->getRmax() && new_distance < PP[iat]->getRmax()))
    {
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      for (int jel = 0; jel < NeighborElecs.size(); ++jel)
        elecTMAffected[NeighborElecs[jel]] = true;
    }
  }
}

void NonLocalECPotential::addComponent(int groupID, NonLocalECPComponent* ppot)
{
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
    {
      PP[iat] = ppot;
      ppot->Ve_sample = &Ve_samp_tmp;
      ppot->Vi_sample = &Vi_samp_tmp;
    }
  PPset[groupID] = ppot;
}

OperatorBase* NonLocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  NonLocalECPotential* myclone = new NonLocalECPotential(IonConfig, qp, psi, ComputeForces);
  for (int ig = 0; ig < PPset.size(); ++ig)
  {
    if (PPset[ig])
    {
      NonLocalECPComponent* ppot = PPset[ig]->makeClone(qp);
      myclone->addComponent(ig, ppot);
    }
  }
  return myclone;
}


void NonLocalECPotential::addObservables(PropertySetType& plist, BufferType& collectables)
{
  OperatorBase::addValue(plist);
  if (ComputeForces)
  {
    if (FirstForceIndex < 0)
      FirstForceIndex = plist.size();
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        std::ostringstream obsName1, obsName2;
        obsName1 << "FNL"
                 << "_" << iat << "_" << x;
        plist.add(obsName1.str());
        //        obsName2 << "FNL_Pulay" << "_" << iat << "_" << x;
        //        plist.add(obsName2.str());
      }
    }
  }
}

void NonLocalECPotential::registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const
{
  OperatorBase::registerObservables(h5list, gid);
  if (ComputeForces)
  {
    std::vector<int> ndim(2);
    ndim[0]                 = Nnuc;
    ndim[1]                 = OHMMS_DIM;
    observable_helper* h5o1 = new observable_helper("FNL");
    h5o1->set_dimensions(ndim, FirstForceIndex);
    h5o1->open(gid);
    h5list.push_back(h5o1);
    //    observable_helper* h5o2 = new observable_helper("FNL_Pulay");
    //    h5o2->set_dimensions(ndim,FirstForceIndex+Nnuc*OHMMS_DIM);
    //    h5o2->open(gid);
    //    h5list.push_back(h5o2);
  }
}

void NonLocalECPotential::setObservables(QMCTraits::PropertySetType& plist)
{
  OperatorBase::setObservables(plist);
  if (ComputeForces)
  {
    int index = FirstForceIndex;
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        //    plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


void NonLocalECPotential::setParticlePropertyList(QMCTraits::PropertySetType& plist, int offset)
{
  OperatorBase::setParticlePropertyList(plist, offset);
  if (ComputeForces)
  {
    int index = FirstForceIndex + offset;
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        //        plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


} // namespace qmcplusplus
