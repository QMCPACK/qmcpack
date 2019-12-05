//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <numeric>
#include <iomanip>
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/createDistanceTable.h"
#include "LongRange/StructFact.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"

//#define PACK_DISTANCETABLES

namespace qmcplusplus
{
//using namespace particle_info;

#ifdef QMC_CUDA
template<>
int ParticleSet::Walker_t::cuda_DataSize = 0;
#endif

const TimerNameList_t<ParticleSet::PSTimers> ParticleSet::PSTimerNames = {{PS_newpos, "ParticleSet::computeNewPosDT"},
                                                                          {PS_donePbyP, "ParticleSet::donePbyP"},
                                                                          {PS_setActive, "ParticleSet::setActive"},
                                                                          {PS_update, "ParticleSet::update"}};

ParticleSet::ParticleSet()
    : quantum_domain(classical),
      IsGrouped(true),
      SameMass(true),
      ThreadID(0),
      activePtcl(-1),
      SK(0),
      myTwist(0.0),
      ParentName("0"),
      TotalNum(0)
{
  initPropertyList();
  setup_timers(myTimers, PSTimerNames, timer_level_fine);
}

ParticleSet::ParticleSet(const ParticleSet& p)
    : IsGrouped(p.IsGrouped),
      SameMass(true),
      ThreadID(0),
      activePtcl(-1),
      mySpecies(p.getSpeciesSet()),
      SK(0),
      myTwist(0.0),
      ParentName(p.parentName())
{
  set_quantum_domain(p.quantum_domain);
  assign(p); //only the base is copied, assumes that other properties are not assignable
  //need explicit copy:
  Mass = p.Mass;
  Z    = p.Z;
  //std::ostringstream o;
  //o<<p.getName()<<ObjectTag;
  //this->setName(o.str());
  //app_log() << "  Copying a particle set " << p.getName() << " to " << this->getName() << " groups=" << groups() << std::endl;
  myName              = p.getName();
  PropertyList.Names  = p.PropertyList.Names;
  PropertyList.Values = p.PropertyList.Values;
  PropertyHistory     = p.PropertyHistory;
  Collectables        = p.Collectables;
  //construct the distance tables with the same order
  for (int i = 0; i < p.DistTables.size(); ++i)
  {
    addTable(p.DistTables[i]->origin(), p.DistTables[i]->DTType);
    DistTables[i]->Need_full_table_loadWalker = p.DistTables[i]->Need_full_table_loadWalker;
  }
  if (p.SK)
  {
    LRBox = p.LRBox;               //copy LRBox
    SK    = new StructFact(*p.SK); //safe to use the copy constructor
    //R.InUnit=p.R.InUnit;
    //createSK();
    //SK->DoUpdate=p.SK->DoUpdate;
  }
  setup_timers(myTimers, PSTimerNames, timer_level_fine);
  myTwist = p.myTwist;

  RSoA = p.RSoA;
  G    = p.G;
  L    = p.L;
}

ParticleSet::~ParticleSet()
{
  DEBUG_MEMORY("ParticleSet::~ParticleSet");
  delete_iter(DistTables.begin(), DistTables.end());
  if (SK)
    delete SK;
}

void ParticleSet::create(int numPtcl) { resize(numPtcl); }

void ParticleSet::create(const std::vector<int>& agroup)
{
  SubPtcl.resize(agroup.size() + 1);
  SubPtcl[0] = 0;
  for (int is = 0; is < agroup.size(); is++)
    SubPtcl[is + 1] = SubPtcl[is] + agroup[is];
  size_t nsum = SubPtcl[agroup.size()];
  resize(nsum);
  TotalNum = nsum;
  int loc  = 0;
  for (int i = 0; i < agroup.size(); i++)
  {
    for (int j = 0; j < agroup[i]; j++, loc++)
      GroupID[loc] = i;
  }
}

void ParticleSet::set_quantum_domain(quantum_domains qdomain)
{
  if (quantum_domain_valid(qdomain))
    quantum_domain = qdomain;
  else
    APP_ABORT("ParticleSet::set_quantum_domain\n  input quantum domain is not valid for particles");
}

void ParticleSet::resetGroups()
{
  int nspecies = mySpecies.getTotalNum();
  if (nspecies == 0)
  {
    APP_ABORT("ParticleSet::resetGroups() Failed. No species exisits");
  }
  int natt = mySpecies.numAttributes();
  int qind = mySpecies.addAttribute("charge");
  if (natt == qind)
  {
    app_log() << " Missing charge attribute of the SpeciesSet " << myName << " particleset" << std::endl;
    app_log() << " Assume neutral particles Z=0.0 " << std::endl;
    for (int ig = 0; ig < nspecies; ig++)
      mySpecies(qind, ig) = 0.0;
  }
  for (int iat = 0; iat < Z.size(); iat++)
    Z[iat] = mySpecies(qind, GroupID[iat]);
  natt        = mySpecies.numAttributes();
  int massind = mySpecies.addAttribute("mass");
  if (massind == natt)
  {
    for (int ig = 0; ig < nspecies; ig++)
      mySpecies(massind, ig) = 1.0;
  }
  SameMass  = true;
  double m0 = mySpecies(massind, 0);
  for (int ig = 1; ig < nspecies; ig++)
    SameMass &= (mySpecies(massind, ig) == m0);
  if (SameMass)
    app_log() << "  All the species have the same mass " << m0 << std::endl;
  else
    app_log() << "  Distinctive masses for each species " << std::endl;
  for (int iat = 0; iat < Mass.size(); iat++)
    Mass[iat] = mySpecies(massind, GroupID[iat]);
  std::vector<int> ng(nspecies, 0);
  for (int iat = 0; iat < GroupID.size(); iat++)
  {
    if (GroupID[iat] < nspecies)
      ng[GroupID[iat]]++;
    else
      APP_ABORT("ParticleSet::resetGroups() Failed. GroupID is out of bound.");
  }
  SubPtcl.resize(nspecies + 1);
  SubPtcl[0] = 0;
  for (int i = 0; i < nspecies; ++i)
    SubPtcl[i + 1] = SubPtcl[i] + ng[i];
  int membersize = mySpecies.addAttribute("membersize");
  for (int ig = 0; ig < nspecies; ++ig)
    mySpecies(membersize, ig) = ng[ig];
  //orgID=ID;
  //orgGroupID=GroupID;
  int new_id = 0;
  for (int i = 0; i < nspecies; ++i)
    for (int iat = 0; iat < GroupID.size(); ++iat)
      if (GroupID[iat] == i)
        IndirectID[new_id++] = ID[iat];
  IsGrouped = true;
  for (int iat = 0; iat < ID.size(); ++iat)
    IsGrouped &= (IndirectID[iat] == ID[iat]);
}

void ParticleSet::randomizeFromSource(ParticleSet& src)
{
  SpeciesSet& srcSpSet(src.getSpeciesSet());
  SpeciesSet& spSet(getSpeciesSet());
  int srcChargeIndx = srcSpSet.addAttribute("charge");
  int srcMemberIndx = srcSpSet.addAttribute("membersize");
  int ChargeIndex   = spSet.addAttribute("charge");
  int MemberIndx    = spSet.addAttribute("membersize");
  int Nsrc          = src.getTotalNum();
  int Nptcl         = getTotalNum();
  int NumSpecies    = spSet.TotalNum;
  int NumSrcSpecies = srcSpSet.TotalNum;
  //Store information about charges and number of each species
  std::vector<int> Zat, Zspec, NofSpecies, NofSrcSpecies, CurElec;
  Zat.resize(Nsrc);
  Zspec.resize(NumSrcSpecies);
  NofSpecies.resize(NumSpecies);
  CurElec.resize(NumSpecies);
  NofSrcSpecies.resize(NumSrcSpecies);
  for (int spec = 0; spec < NumSrcSpecies; spec++)
  {
    Zspec[spec]         = (int)round(srcSpSet(srcChargeIndx, spec));
    NofSrcSpecies[spec] = (int)round(srcSpSet(srcMemberIndx, spec));
  }
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    NofSpecies[spec] = (int)round(spSet(MemberIndx, spec));
    CurElec[spec]    = first(spec);
  }
  int totQ = 0;
  for (int iat = 0; iat < Nsrc; iat++)
    totQ += Zat[iat] = Zspec[src.GroupID[iat]];
  app_log() << "  Total ion charge    = " << totQ << std::endl;
  totQ -= Nptcl;
  app_log() << "  Total system charge = " << totQ << std::endl;
  // Now, loop over ions, attaching electrons to them to neutralize
  // charge
  int spToken = 0;
  // This is decremented when we run out of electrons in each species
  int spLeft = NumSpecies;
  std::vector<PosType> gaussRand(Nptcl);
  makeGaussRandom(gaussRand);
  for (int iat = 0; iat < Nsrc; iat++)
  {
    // Loop over electrons to add, selecting round-robin from the
    // electron species
    int z = Zat[iat];
    while (z > 0 && spLeft)
    {
      int sp = spToken++ % NumSpecies;
      if (NofSpecies[sp])
      {
        NofSpecies[sp]--;
        z--;
        int elec = CurElec[sp]++;
        app_log() << "  Assigning " << (sp ? "down" : "up  ") << " electron " << elec << " to ion " << iat
                  << " with charge " << z << std::endl;
        double radius = 0.5 * std::sqrt((double)Zat[iat]);
        R[elec]       = src.R[iat] + radius * gaussRand[elec];
      }
      else
        spLeft--;
    }
  }
  // Assign remaining electrons
  int ion = 0;
  for (int sp = 0; sp < NumSpecies; sp++)
  {
    for (int ie = 0; ie < NofSpecies[sp]; ie++)
    {
      int iat       = ion++ % Nsrc;
      double radius = std::sqrt((double)Zat[iat]);
      int elec      = CurElec[sp]++;
      R[elec]       = src.R[iat] + radius * gaussRand[elec];
    }
  }
}

///write to a std::ostream
bool ParticleSet::get(std::ostream& os) const
{
  os << "  ParticleSet '" << getName() << "' contains " << TotalNum << " particles : ";
  if (SubPtcl.size() > 0)
    for (int i = 0; i < SubPtcl.size() - 1; i++)
      os << " " << mySpecies.speciesName[i] << "(" << SubPtcl[i + 1] - SubPtcl[i] << ")";
  os << std::endl;
  if (!IsGrouped)
    os << "    Particles are not grouped by species in the input file. Algorithms may not be optimal!" << std::endl;
  os << std::endl;

  const size_t maxParticlesToPrint = 10;
  size_t numToPrint                = std::min(TotalNum, maxParticlesToPrint);

  for (int i = 0; i < numToPrint; i++)
  {
    os << "    " << mySpecies.speciesName[GroupID[i]] << R[i] << std::endl;
  }
  if (numToPrint < TotalNum)
  {
    os << "    (... and " << (TotalNum - numToPrint) << " more particle positions ...)" << std::endl;
  }
  os << std::endl;

  for (const std::string& description : distTableDescriptions)
    os << description;
  os << std::endl;
  return true;
}

///read from std::istream
bool ParticleSet::put(std::istream& is) { return true; }

///reset member data
void ParticleSet::reset() { app_log() << "<<<< going to set properties >>>> " << std::endl; }

///read the particleset
bool ParticleSet::put(xmlNodePtr cur) { return true; }

int ParticleSet::addTable(const ParticleSet& psrc, int dt_type, bool need_full_table_loadWalker)
{
  if (myName == "none" || psrc.getName() == "none")
    APP_ABORT("ParticleSet::addTable needs proper names for both source and target particle sets.");

  if (DistTables.size() > 0 && dt_type != DT_SOA_PREFERRED && !DistTables[0]->is_same_type(dt_type))
    APP_ABORT("ParticleSet::addTable Cannot mix AoS and SoA distance tables.\n");

  int tid;
  std::map<std::string, int>::iterator tit(myDistTableMap.find(psrc.getName()));
  if (tit == myDistTableMap.end())
  {
    std::ostringstream description;
    tid                = DistTables.size();
    int dt_type_in_use = (tid == 0 ? dt_type : DistTables[0]->DTType);
    if (myName == psrc.getName())
      DistTables.push_back(createDistanceTable(*this, dt_type_in_use, description));
    else
      DistTables.push_back(createDistanceTable(psrc, *this, dt_type_in_use, description));
    distTableDescriptions.push_back(description.str());
    myDistTableMap[psrc.getName()] = tid;
    app_debug() << "  ... ParticleSet::addTable Create Table #" << tid << " " << DistTables[tid]->Name << std::endl;
  }
  else
  {
    tid = (*tit).second;
    app_debug() << "  ... ParticleSet::addTable Reuse Table #" << tid << " " << DistTables[tid]->Name << std::endl;
  }
  DistTables[tid]->Need_full_table_loadWalker =
      (DistTables[tid]->Need_full_table_loadWalker || need_full_table_loadWalker);
  app_log().flush();
  return tid;
}

void ParticleSet::update(bool skipSK)
{
  ScopedTimer update_scope(myTimers[PS_update]);

  RSoA.copyIn(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (!skipSK && SK)
    SK->UpdateAllPart(*this);

  activePtcl = -1;
}

void ParticleSet::setActive(int iat)
{
  ScopedTimer set_active_scope(myTimers[PS_setActive]);

  for (size_t i = 0; i < DistTables.size(); i++)
    if (DistTables[i]->DTType == DT_SOA)
      DistTables[i]->evaluate(*this, iat);
}

void ParticleSet::flex_setActive(const RefVector<ParticleSet>& P_list, int iat)
{
  if (P_list.size() > 1)
  {
    ScopedTimer local_timer(P_list[0].get().myTimers[PS_setActive]);
    int dist_tables_size = P_list[0].get().DistTables.size();
#pragma omp parallel
    {
      for (size_t i = 0; i < dist_tables_size; i++)
      {
#pragma omp for
        for (int iw = 0; iw < P_list.size(); iw++)
        {
          P_list[iw].get().DistTables[i]->evaluate(P_list[iw], iat);
        }
      }
    }
  }
  else if (P_list.size() == 1)
    P_list[0].get().setActive(iat);
}

void ParticleSet::makeMove(Index_t iat, const SingleParticlePos_t& displ)
{
  activePtcl = iat;
  activePos  = R[iat] + displ;
  computeNewPosDistTablesAndSK(iat, activePos);
}

void ParticleSet::makeMoveWithSpin(Index_t iat, const SingleParticlePos_t& displ, const RealType& sdispl)
{
  makeMove(iat, displ);
  activeSpinVal = spins[iat] + sdispl;
}

void ParticleSet::flex_makeMove(const RefVector<ParticleSet>& P_list,
                                Index_t iat,
                                const std::vector<SingleParticlePos_t>& displs)
{
  if (P_list.size() > 1)
  {
    std::vector<SingleParticlePos_t> new_positions;
    new_positions.reserve(displs.size());

    for (int iw = 0; iw < P_list.size(); iw++)
    {
      P_list[iw].get().activePtcl = iat;
      P_list[iw].get().activePos  = P_list[iw].get().R[iat] + displs[iw];
      new_positions.push_back(P_list[iw].get().activePos);
    }

    mw_computeNewPosDistTablesAndSK(P_list, iat, new_positions);
  }
  else if (P_list.size() == 1)
    P_list[0].get().makeMove(iat, displs[0]);
}

bool ParticleSet::makeMoveAndCheck(Index_t iat, const SingleParticlePos_t& displ)
{
  activePtcl    = iat;
  activePos     = R[iat] + displ;
  bool is_valid = true;
  if (Lattice.explicitly_defined)
  {
    if (Lattice.outOfBound(Lattice.toUnit(displ)))
      is_valid = false;
    else
    {
      newRedPos = Lattice.toUnit(activePos);
      if (!Lattice.isValid(newRedPos))
        is_valid = false;
    }
  }
  computeNewPosDistTablesAndSK(iat, activePos);
  return is_valid;
}

bool ParticleSet::makeMoveAndCheckWithSpin(Index_t iat, const SingleParticlePos_t& displ, const RealType& sdispl)
{
  activeSpinVal = spins[iat] + sdispl;
  return makeMoveAndCheck(iat, displ);
}

void ParticleSet::computeNewPosDistTablesAndSK(Index_t iat, const SingleParticlePos_t& newpos)
{
  ScopedTimer compute_newpos_scope(myTimers[PS_newpos]);

  for (int i = 0; i < DistTables.size(); ++i)
    DistTables[i]->move(*this, newpos);
  //Do not change SK: 2007-05-18
  //Change SK only if DoUpdate is true: 2008-09-12
  if (SK && SK->DoUpdate)
    SK->makeMove(iat, newpos);
}

void ParticleSet::mw_computeNewPosDistTablesAndSK(const RefVector<ParticleSet>& P_list,
                                                  Index_t iat,
                                                  const std::vector<SingleParticlePos_t>& new_positions)
{
  ScopedTimer compute_newpos_scope(P_list[0].get().myTimers[PS_newpos]);
  int dist_tables_size = P_list[0].get().DistTables.size();
#pragma omp parallel
  {
    for (int i = 0; i < dist_tables_size; ++i)
    {
#pragma omp for
      for (int iw = 0; iw < P_list.size(); iw++)
        P_list[iw].get().DistTables[i]->move(P_list[iw], new_positions[iw]);
    }

    StructFact* SK = P_list[0].get().SK;
    if (SK && SK->DoUpdate)
    {
#pragma omp for
      for (int iw = 0; iw < P_list.size(); iw++)
        P_list[iw].get().SK->makeMove(iat, new_positions[iw]);
    }
  }
}


bool ParticleSet::makeMoveAllParticles(const Walker_t& awalker, const ParticlePos_t& deltaR, RealType dt)
{
  activePtcl = -1;
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt * deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat] + displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat] = newpos;
    }
  }
  else
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
      R[iat] = awalker.R[iat] + dt * deltaR[iat];
  }
  RSoA.copyIn(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMoveAllParticles(const Walker_t& awalker,
                                       const ParticlePos_t& deltaR,
                                       const std::vector<RealType>& dt)
{
  activePtcl = -1;
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt[iat] * deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat] + displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat] = newpos;
    }
  }
  else
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
      R[iat] = awalker.R[iat] + dt[iat] * deltaR[iat];
  }
  RSoA.copyIn(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

/** move a walker by dt*deltaR + drift
 * @param awalker initial walker configuration
 * @param drift drift vector
 * @param deltaR random displacement
 * @param dt timestep
 * @return true, if all the particle moves are legal under the boundary conditions
 */
bool ParticleSet::makeMoveAllParticlesWithDrift(const Walker_t& awalker,
                                                const ParticlePos_t& drift,
                                                const ParticlePos_t& deltaR,
                                                RealType dt)
{
  activePtcl = -1;
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt * deltaR[iat] + drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat] + displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat] = newpos;
    }
  }
  else
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
      R[iat] = awalker.R[iat] + dt * deltaR[iat] + drift[iat];
  }
  RSoA.copyIn(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMoveAllParticlesWithDrift(const Walker_t& awalker,
                                                const ParticlePos_t& drift,
                                                const ParticlePos_t& deltaR,
                                                const std::vector<RealType>& dt)
{
  activePtcl = -1;
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt[iat] * deltaR[iat] + drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat] + displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat] = newpos;
    }
  }
  else
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
      R[iat] = awalker.R[iat] + dt[iat] * deltaR[iat] + drift[iat];
  }
  RSoA.copyIn(R);

  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

/** update the particle attribute by the proposed move
 *@param iat the particle index
 *
 *When the activePtcl is equal to iat, overwrite the position and update the
 *content of the distance tables.
 */
void ParticleSet::acceptMove(Index_t iat)
{
  if (iat == activePtcl)
  {
    //Update position + distance-table
    for (int i = 0, n = DistTables.size(); i < n; i++)
      DistTables[i]->update(iat);

    //Do not change SK: 2007-05-18
    if (SK && SK->DoUpdate)
      SK->acceptMove(iat, GroupID[iat], R[iat]);

    R[iat]     = activePos;
    RSoA(iat)  = activePos;
    spins[iat] = activeSpinVal;
    activePtcl = -1;
  }
  else
  {
    std::ostringstream o;
    o << "  Illegal acceptMove " << iat << " != " << activePtcl;
    APP_ABORT(o.str());
  }
}

void ParticleSet::flex_donePbyP(const RefVector<ParticleSet>& P_list)
{
  for (int iw = 0; iw < P_list.size(); iw++)
    P_list[iw].get().donePbyP();
}

void ParticleSet::donePbyP()
{
  ScopedTimer donePbyP_scope(myTimers[PS_donePbyP]);
  if (SK && !SK->DoUpdate)
    SK->UpdateAllPart(*this);
  activePtcl = -1;
}

void ParticleSet::makeVirtualMoves(const SingleParticlePos_t& newpos)
{
  activePtcl = -1;
  activePos  = newpos;
  for (size_t i = 0; i < DistTables.size(); ++i)
    DistTables[i]->move(*this, newpos);
}

void ParticleSet::loadWalker(Walker_t& awalker, bool pbyp)
{
  R = awalker.R;
  RSoA.copyIn(R);
#if !defined(SOA_MEMORY_OPTIMIZED)
  G = awalker.G;
  L = awalker.L;
#endif
  if (pbyp)
  {
    // in certain cases, full tables must be ready
    for (int i = 0; i < DistTables.size(); i++)
      if (DistTables[i]->DTType == DT_AOS || DistTables[i]->Need_full_table_loadWalker)
        DistTables[i]->evaluate(*this);
    //computed so that other objects can use them, e.g., kSpaceJastrow
    if (SK && SK->DoUpdate)
      SK->UpdateAllPart(*this);
  }

  activePtcl = -1;
}

void ParticleSet::saveWalker(Walker_t& awalker)
{
  awalker.R = R;
#if !defined(SOA_MEMORY_OPTIMIZED)
  awalker.G = G;
  awalker.L = L;
#endif
  //PAOps<RealType,OHMMS_DIM>::copy(G,awalker.Drift);
  //if (SK)
  //  SK->UpdateAllPart(*this);
  //awalker.DataSet.rewind();
}

void ParticleSet::flex_saveWalker(RefVector<ParticleSet>& psets, RefVector<Walker_t>& walkers)
{
  int num_sets    = psets.size();
  auto saveWalker = [](ParticleSet& pset, Walker_t& walker) {
    walker.R = pset.R;
#if !defined(SOA_MEMORY_OPTIMIZED)
    walker.G = pset.G;
    walker.L = pset.L;
#endif
  };
  for (int iw = 0; iw < num_sets; ++iw)
    saveWalker(psets[iw], walkers[iw]);
}


void ParticleSet::initPropertyList()
{
  PropertyList.clear();
  //Need to add the default Properties according to the enumeration
  PropertyList.add("LogPsi");
  PropertyList.add("SignPsi");
  PropertyList.add("UmbrellaWeight");
  PropertyList.add("R2Accepted");
  PropertyList.add("R2Proposed");
  PropertyList.add("DriftScale");
  PropertyList.add("AltEnergy");
  PropertyList.add("LocalEnergy");
  PropertyList.add("LocalPotential");
  if (PropertyList.size() != NUMPROPERTIES)
  {
    app_error() << "The number of default properties for walkers  is not consistent." << std::endl;
    app_error() << "NUMPROPERTIES " << NUMPROPERTIES << " size of PropertyList " << PropertyList.size() << std::endl;
    APP_ABORT("ParticleSet::initPropertyList");
  }
}

void ParticleSet::clearDistanceTables()
{
  //Physically remove the tables
  delete_iter(DistTables.begin(), DistTables.end());
  DistTables.clear();
  //for(int i=0; i< DistTables.size(); i++) DistanceTable::removeTable(DistTables[i]->getName());
  //DistTables.erase(DistTables.begin(),DistTables.end());
}

int ParticleSet::addPropertyHistory(int leng)
{
  int newL                                    = PropertyHistory.size();
  std::vector<FullPrecRealType> newVecHistory = std::vector<FullPrecRealType>(leng, 0.0);
  PropertyHistory.push_back(newVecHistory);
  PHindex.push_back(0);
  return newL;
}

//      void ParticleSet::resetPropertyHistory( )
//     {
//       for(int i=0;i<PropertyHistory.size();i++)
//       {
//         PHindex[i]=0;
//  for(int k=0;k<PropertyHistory[i].size();k++)
//  {
//    PropertyHistory[i][k]=0.0;
//  }
//       }
//     }

//      void ParticleSet::addPropertyHistoryPoint(int index, RealType data)
//     {
//       PropertyHistory[index][PHindex[index]]=(data);
//       PHindex[index]++;
//       if (PHindex[index]==PropertyHistory[index].size()) PHindex[index]=0;
// //       PropertyHistory[index].pop_back();
//     }

//      void ParticleSet::rejectedMove()
//     {
//       for(int dindex=0;dindex<PropertyHistory.size();dindex++){
//         int lastIndex=PHindex[dindex]-1;
//         if (lastIndex<0) lastIndex+=PropertyHistory[dindex].size();
//         PropertyHistory[dindex][PHindex[dindex]]=PropertyHistory[dindex][lastIndex];
//         PHindex[dindex]++;
//         if (PHindex[dindex]==PropertyHistory[dindex].size()) PHindex[dindex]=0;
// //       PropertyHistory[dindex].push_front(PropertyHistory[dindex].front());
// //       PropertyHistory[dindex].pop_back();
//       }
//     }

} // namespace qmcplusplus
