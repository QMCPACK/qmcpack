//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
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
#include "Particle/DynamicCoordinatesBuilder.h"
#include "Particle/DistanceTableData.h"
#include "Particle/createDistanceTable.h"
#include "LongRange/StructFact.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"

//#define PACK_DISTANCETABLES

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

//using namespace particle_info;

#ifdef QMC_CUDA
template<>
int ParticleSet::Walker_t::cuda_DataSize = 0;
#endif

const TimerNameList_t<ParticleSet::PSTimers> ParticleSet::PSTimerNames = {{PS_newpos, "ParticleSet::computeNewPosDT"},
                                                                          {PS_donePbyP, "ParticleSet::donePbyP"},
                                                                          {PS_accept, "ParticleSet::acceptMove"},
                                                                          {PS_update, "ParticleSet::update"}};

ParticleSet::ParticleSet(const DynamicCoordinateKind kind)
    : quantum_domain(classical),
      IsGrouped(true),
      SameMass(true),
      ThreadID(0),
      activePtcl(-1),
      SK(0),
      Properties(0, 0, 1, WP::MAXPROPERTIES),
      myTwist(0.0),
      ParentName("0"),
      TotalNum(0),
      coordinates_(std::move(createDynamicCoordinates(kind)))
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
      Properties(p.Properties),
      myTwist(0.0),
      ParentName(p.parentName()),
      coordinates_(std::move(p.coordinates_->makeClone()))
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
    addTable(p.DistTables[i]->origin(), p.DistTables[i]->DTType, p.DistTables[i]->getFullTableNeeds());
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

  G = p.G;
  L = p.L;
}

ParticleSet::~ParticleSet()
{
  DEBUG_MEMORY("ParticleSet::~ParticleSet");
  delete_iter(DistTables.begin(), DistTables.end());
  if (SK)
    delete SK;
}

void ParticleSet::create(int numPtcl)
{
  resize(numPtcl);
  SubPtcl.resize(2);
  SubPtcl[0] = 0;
  SubPtcl[1] = numPtcl;
}

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
  // Usually an empty ParticleSet indicates an error in the input file,
  // but in some cases it is useful.  Allow an empty ParticleSet if it
  // has the special name "empty".
  if (nspecies == 0 && getName() != "empty")
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
  // safety check if any group of particles has size 0, instruct users to fix the input.
  for (int group_id = 0; group_id < nspecies; group_id++)
    if (ng[group_id] == 0 && getName() != "empty")
    {
      std::ostringstream err_msg;
      err_msg << "ParticleSet::resetGroups() Failed. ParticleSet '" << myName << "' "
              << "has group '" << mySpecies.speciesName[group_id] << "' containing 0 particles. "
              << "Remove this group from input!" << std::endl;
      APP_ABORT(err_msg.str());
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

int ParticleSet::addTable(const ParticleSet& psrc, int dt_type, bool need_full_table)
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
    app_debug() << "  ... ParticleSet::addTable Create Table #" << tid << " " << DistTables[tid]->getName()
                << std::endl;
  }
  else
  {
    tid = (*tit).second;
    app_debug() << "  ... ParticleSet::addTable Reuse Table #" << tid << " " << DistTables[tid]->getName() << std::endl;
  }

  DistTables[tid]->setFullTableNeeds(DistTables[tid]->getFullTableNeeds() || need_full_table);

  app_log().flush();
  return tid;
}

void ParticleSet::update(bool skipSK)
{
  ScopedTimer update_scope(myTimers[PS_update]);

  coordinates_->setAllParticlePos(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (!skipSK && SK)
    SK->UpdateAllPart(*this);

  activePtcl = -1;
}

void ParticleSet::flex_update(const RefVector<ParticleSet>& p_list, bool skipSK)
{
  if (p_list.size() > 1)
  {
    ScopedTimer update_scope(p_list[0].get().myTimers[PS_update]);

    for (ParticleSet& pset : p_list)
      pset.setCoordinates(pset.R);

    auto& dts = p_list[0].get().DistTables;
    for (int i = 0; i < dts.size(); ++i)
    {
      const auto dt_list(extractDTRefList(p_list, i));
      dts[i]->mw_evaluate(dt_list, p_list);
    }

    if (!skipSK && p_list[0].get().SK)
    {
#pragma omp parallel for
      for (int iw = 0; iw < p_list.size(); iw++)
        p_list[iw].get().SK->UpdateAllPart(p_list[iw]);
    }
  }
  else if (p_list.size() == 1)
    p_list[0].get().update(skipSK);
}

void ParticleSet::makeMove(Index_t iat, const SingleParticlePos_t& displ, bool maybe_accept)
{
  activePtcl    = iat;
  activePos     = R[iat] + displ;
  activeSpinVal = spins[iat];
  computeNewPosDistTablesAndSK(iat, activePos, maybe_accept);
}

void ParticleSet::makeMoveWithSpin(Index_t iat, const SingleParticlePos_t& displ, const Scalar_t& sdispl)
{
  makeMove(iat, displ);
  activeSpinVal += sdispl;
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
  activeSpinVal = spins[iat];
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
  computeNewPosDistTablesAndSK(iat, activePos, true);
  return is_valid;
}

bool ParticleSet::makeMoveAndCheckWithSpin(Index_t iat, const SingleParticlePos_t& displ, const Scalar_t& sdispl)
{
  bool is_valid = makeMoveAndCheck(iat, displ);
  activeSpinVal += sdispl;
  return is_valid;
}

void ParticleSet::computeNewPosDistTablesAndSK(Index_t iat, const SingleParticlePos_t& newpos, bool maybe_accept)
{
  ScopedTimer compute_newpos_scope(myTimers[PS_newpos]);

  for (int i = 0; i < DistTables.size(); ++i)
    DistTables[i]->move(*this, newpos, iat, maybe_accept);
  //Do not change SK: 2007-05-18
  //Change SK only if DoUpdate is true: 2008-09-12
  if (SK && SK->DoUpdate)
    SK->makeMove(iat, newpos);
}

void ParticleSet::mw_computeNewPosDistTablesAndSK(const RefVector<ParticleSet>& P_list,
                                                  Index_t iat,
                                                  const std::vector<SingleParticlePos_t>& new_positions,
                                                  bool maybe_accept)
{
  ScopedTimer compute_newpos_scope(P_list[0].get().myTimers[PS_newpos]);
  const int dist_tables_size = P_list[0].get().DistTables.size();
#pragma omp parallel
  {
    for (int i = 0; i < dist_tables_size; ++i)
    {
#pragma omp for
      for (int iw = 0; iw < P_list.size(); iw++)
        P_list[iw].get().DistTables[i]->move(P_list[iw], new_positions[iw], iat, maybe_accept);
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
  coordinates_->setAllParticlePos(R);
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
  coordinates_->setAllParticlePos(R);
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
  coordinates_->setAllParticlePos(R);
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
  coordinates_->setAllParticlePos(R);

  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

/** update the particle attribute by the proposed move
 *
 * When the activePtcl is equal to iat, overwrite the position and update the
 * content of the distance tables.
 */
void ParticleSet::acceptMove(Index_t iat, bool partial_table_update)
{
  ScopedTimer update_scope(myTimers[PS_accept]);
  if (iat == activePtcl)
  {
    //Update position + distance-table
    for (int i = 0, n = DistTables.size(); i < n; i++)
      DistTables[i]->update(iat, partial_table_update);

    //Do not change SK: 2007-05-18
    if (SK && SK->DoUpdate)
      SK->acceptMove(iat, GroupID[iat], R[iat]);

    R[iat] = activePos;
    coordinates_->setOneParticlePos(activePos, iat);
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
  if (P_list.size() > 1)
  {
    // Leaving bare omp pragma here. It can potentially be improved with cleaner abstraction.
    #pragma omp parallel for
    for (int iw = 0; iw < P_list.size(); iw++)
      P_list[iw].get().donePbyP();
  }
  else if (P_list.size() == 1)
    P_list[0].get().donePbyP();
}

void ParticleSet::donePbyP()
{
  ScopedTimer donePbyP_scope(myTimers[PS_donePbyP]);
  coordinates_->donePbyP();
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
  coordinates_->setAllParticlePos(R);
#if !defined(SOA_MEMORY_OPTIMIZED)
  G = awalker.G;
  L = awalker.L;
#endif
  if (pbyp)
  {
    ScopedTimer update_scope(myTimers[PS_update]);

    // in certain cases, full tables must be ready
    for (int i = 0; i < DistTables.size(); i++)
      if (DistTables[i]->DTType == DT_AOS || DistTables[i]->getFullTableNeeds())
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
  awalker.spins = spins;
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

  // There is no point in checking this, its quickly not consistent as other objects update property list.
  // if (PropertyList.size() != WP::NUMPROPERTIES)
  // {
  //   app_error() << "The number of default properties for walkers  is not consistent." << std::endl;
  //   app_error() << "NUMPROPERTIES " << WP::NUMPROPERTIES << " size of PropertyList " << PropertyList.size() << std::endl;
  //   APP_ABORT("ParticleSet::initPropertyList");
  // }
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

RefVector<DistanceTableData> ParticleSet::extractDTRefList(
    const RefVector<ParticleSet>& p_list,
    int id)
{
  RefVector<DistanceTableData> dt_list;
  dt_list.reserve(p_list.size());
  for (ParticleSet& p : p_list)
    dt_list.push_back(*p.DistTables[id]);
  return dt_list;
}

} // namespace qmcplusplus
