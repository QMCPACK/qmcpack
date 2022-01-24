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
#include "ParticleSet.h"
#include "Particle/DynamicCoordinatesBuilder.h"
#include "Particle/DistanceTable.h"
#include "Particle/createDistanceTable.h"
#include "LongRange/StructFact.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGeneratorGlobal.h"

//#define PACK_DISTANCETABLES

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

//using namespace particle_info;

#ifdef QMC_CUDA
template<>
int ParticleSet::Walker_t::cuda_DataSize = 0;
#endif

enum PSetTimers
{
  PS_newpos,
  PS_donePbyP,
  PS_accept,
  PS_loadWalker,
  PS_update,
  PS_dt_move,
  PS_mw_copy
};

static const TimerNameList_t<PSetTimers> generatePSetTimerNames(std::string& obj_name)
{
  return {{PS_newpos, "ParticleSet:" + obj_name + "::computeNewPosDT"},
          {PS_donePbyP, "ParticleSet:" + obj_name + "::donePbyP"},
          {PS_accept, "ParticleSet:" + obj_name + "::acceptMove"},
          {PS_loadWalker, "ParticleSet:" + obj_name + "::loadWalker"},
          {PS_update, "ParticleSet:" + obj_name + "::update"},
          {PS_dt_move, "ParticleSet:" + obj_name + "::dt_move"},
          {PS_mw_copy, "ParticleSet:" + obj_name + "::mw_copy"}};
}

ParticleSet::ParticleSet(const SimulationCell& simulation_cell, const DynamicCoordinateKind kind)
    : quantum_domain(classical),
      Properties(0, 0, 1, WP::MAXPROPERTIES),
      simulation_cell_(simulation_cell),
      is_grouped_(true),
      same_mass_(true),
      is_spinor_(false),
      active_ptcl_(-1),
      myTwist(0.0),
      ParentName("0"),
      TotalNum(0),
      coordinates_(createDynamicCoordinates(kind))
{
  initPropertyList();
  setup_timers(myTimers, generatePSetTimerNames(myName), timer_level_medium);
}

ParticleSet::ParticleSet(const ParticleSet& p)
    : Properties(p.Properties),
      simulation_cell_(p.simulation_cell_),
      is_grouped_(p.is_grouped_),
      same_mass_(true),
      is_spinor_(false),
      active_ptcl_(-1),
      my_species_(p.getSpeciesSet()),
      myTwist(0.0),
      ParentName(p.parentName()),
      coordinates_(p.coordinates_->makeClone())
{
  setQuantumDomain(p.quantum_domain);
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
    addTable(p.DistTables[i]->get_origin(), p.DistTables[i]->getModes());

  if (p.structure_factor_)
    structure_factor_ = std::make_unique<StructFact>(*p.structure_factor_);
  setup_timers(myTimers, generatePSetTimerNames(myName), timer_level_medium);
  myTwist = p.myTwist;

  G = p.G;
  L = p.L;
}

ParticleSet::~ParticleSet() = default;

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

void ParticleSet::setQuantumDomain(quantum_domains qdomain)
{
  if (quantumDomainValid(qdomain))
    quantum_domain = qdomain;
  else
    throw std::runtime_error("ParticleSet::setQuantumDomain\n  input quantum domain is not valid for particles");
}

void ParticleSet::resetGroups()
{
  int nspecies = my_species_.getTotalNum();
  // Usually an empty ParticleSet indicates an error in the input file,
  // but in some cases it is useful.  Allow an empty ParticleSet if it
  // has the special name "empty".
  if (nspecies == 0 && getName() != "empty")
  {
    throw std::runtime_error("ParticleSet::resetGroups() Failed. No species exisits");
  }
  int natt = my_species_.numAttributes();
  int qind = my_species_.addAttribute("charge");
  if (natt == qind)
  {
    app_log() << " Missing charge attribute of the SpeciesSet " << myName << " particleset" << std::endl;
    app_log() << " Assume neutral particles Z=0.0 " << std::endl;
    for (int ig = 0; ig < nspecies; ig++)
      my_species_(qind, ig) = 0.0;
  }
  for (int iat = 0; iat < Z.size(); iat++)
    Z[iat] = my_species_(qind, GroupID[iat]);
  natt        = my_species_.numAttributes();
  int massind = my_species_.addAttribute("mass");
  if (massind == natt)
  {
    for (int ig = 0; ig < nspecies; ig++)
      my_species_(massind, ig) = 1.0;
  }
  same_mass_ = true;
  double m0  = my_species_(massind, 0);
  for (int ig = 1; ig < nspecies; ig++)
    same_mass_ &= (my_species_(massind, ig) == m0);
  if (same_mass_)
    app_log() << "  All the species have the same mass " << m0 << std::endl;
  else
    app_log() << "  Distinctive masses for each species " << std::endl;
  for (int iat = 0; iat < Mass.size(); iat++)
    Mass[iat] = my_species_(massind, GroupID[iat]);
  std::vector<int> ng(nspecies, 0);
  for (int iat = 0; iat < GroupID.size(); iat++)
  {
    if (GroupID[iat] < nspecies)
      ng[GroupID[iat]]++;
    else
      throw std::runtime_error("ParticleSet::resetGroups() Failed. GroupID is out of bound.");
  }
  // safety check if any group of particles has size 0, instruct users to fix the input.
  for (int group_id = 0; group_id < nspecies; group_id++)
    if (ng[group_id] == 0 && getName() != "empty")
    {
      std::ostringstream err_msg;
      err_msg << "ParticleSet::resetGroups() Failed. ParticleSet '" << myName << "' "
              << "has group '" << my_species_.speciesName[group_id] << "' containing 0 particles. "
              << "Remove this group from input!" << std::endl;
      throw std::runtime_error(err_msg.str());
    }
  SubPtcl.resize(nspecies + 1);
  SubPtcl[0] = 0;
  for (int i = 0; i < nspecies; ++i)
    SubPtcl[i + 1] = SubPtcl[i] + ng[i];
  int membersize = my_species_.addAttribute("membersize");
  for (int ig = 0; ig < nspecies; ++ig)
    my_species_(membersize, ig) = ng[ig];
  /** ID map that reflects species group
   *
   * IsGrouped=true, if ID==IndirectID
   */
  ParticleIndex IndirectID(R.size());
  int new_id = 0;
  for (int i = 0; i < nspecies; ++i)
    for (int iat = 0; iat < GroupID.size(); ++iat)
      if (GroupID[iat] == i)
        IndirectID[new_id++] = iat;
  is_grouped_ = true;
  for (int iat = 0; iat < IndirectID.size(); ++iat)
    is_grouped_ &= (IndirectID[iat] == iat);
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
      os << " " << my_species_.speciesName[i] << "(" << SubPtcl[i + 1] - SubPtcl[i] << ")";
  os << std::endl;
  if (!is_grouped_)
    os << "    Particles are not grouped by species in the input file. Algorithms may not be optimal!" << std::endl;
  os << std::endl;

  const size_t maxParticlesToPrint = 10;
  size_t numToPrint                = std::min(TotalNum, maxParticlesToPrint);

  for (int i = 0; i < numToPrint; i++)
  {
    os << "    " << my_species_.speciesName[GroupID[i]] << R[i] << std::endl;
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

int ParticleSet::addTable(const ParticleSet& psrc, DTModes modes)
{
  if (myName == "none" || psrc.getName() == "none")
    throw std::runtime_error("ParticleSet::addTable needs proper names for both source and target particle sets.");

  int tid;
  std::map<std::string, int>::iterator tit(myDistTableMap.find(psrc.getName()));
  if (tit == myDistTableMap.end())
  {
    std::ostringstream description;
    tid = DistTables.size();
    if (myName == psrc.getName())
      DistTables.push_back(createDistanceTable(*this, description));
    else
      DistTables.push_back(createDistanceTable(psrc, *this, description));
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

  DistTables[tid]->setModes(DistTables[tid]->getModes() | modes);

  app_log().flush();
  return tid;
}

const DistanceTableAA& ParticleSet::getDistTableAA(int table_ID) const
{
  return dynamic_cast<DistanceTableAA&>(*DistTables[table_ID]);
}

const DistanceTableAB& ParticleSet::getDistTableAB(int table_ID) const
{
  return dynamic_cast<DistanceTableAB&>(*DistTables[table_ID]);
}

void ParticleSet::update(bool skipSK)
{
  ScopedTimer update_scope(myTimers[PS_update]);

  coordinates_->setAllParticlePos(R);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (!skipSK && structure_factor_)
    structure_factor_->updateAllPart(*this);

  active_ptcl_ = -1;
}

void ParticleSet::mw_update(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK)
{
  auto& p_leader = p_list.getLeader();
  ScopedTimer update_scope(p_leader.myTimers[PS_update]);

  for (ParticleSet& pset : p_list)
    pset.setCoordinates(pset.R);

  auto& dts = p_leader.DistTables;
  for (int i = 0; i < dts.size(); ++i)
  {
    const auto dt_list(extractDTRefList(p_list, i));
    dts[i]->mw_evaluate(dt_list, p_list);
  }

  if (!skipSK && p_leader.structure_factor_)
  {
#pragma omp parallel for
    for (int iw = 0; iw < p_list.size(); iw++)
      p_list[iw].structure_factor_->updateAllPart(p_list[iw]);
  }
}

void ParticleSet::makeMove(Index_t iat, const SingleParticlePos& displ, bool maybe_accept)
{
  active_ptcl_     = iat;
  active_pos_      = R[iat] + displ;
  active_spin_val_ = spins[iat];
  computeNewPosDistTablesAndSK(iat, active_pos_, maybe_accept);
}

void ParticleSet::makeMoveWithSpin(Index_t iat, const SingleParticlePos& displ, const Scalar_t& sdispl)
{
  makeMove(iat, displ);
  active_spin_val_ += sdispl;
}

void ParticleSet::mw_makeMove(const RefVectorWithLeader<ParticleSet>& p_list,
                              Index_t iat,
                              const std::vector<SingleParticlePos>& displs)
{
  std::vector<SingleParticlePos> new_positions;
  new_positions.reserve(displs.size());

  for (int iw = 0; iw < p_list.size(); iw++)
  {
    p_list[iw].active_ptcl_     = iat;
    p_list[iw].active_pos_      = p_list[iw].R[iat] + displs[iw];
    p_list[iw].active_spin_val_ = p_list[iw].spins[iat];
    new_positions.push_back(p_list[iw].active_pos_);
  }

  mw_computeNewPosDistTablesAndSK(p_list, iat, new_positions);
}

void ParticleSet::mw_makeMoveWithSpin(const RefVectorWithLeader<ParticleSet>& p_list,
                                      Index_t iat,
                                      const std::vector<SingleParticlePos>& displs,
                                      const std::vector<Scalar_t>& sdispls)
{
  mw_makeMove(p_list, iat, displs);
  for (int iw = 0; iw < p_list.size(); iw++)
    p_list[iw].active_spin_val_ += sdispls[iw];
}


bool ParticleSet::makeMoveAndCheck(Index_t iat, const SingleParticlePos& displ)
{
  active_ptcl_     = iat;
  active_pos_      = R[iat] + displ;
  active_spin_val_ = spins[iat];
  bool is_valid    = true;
  auto& Lattice    = simulation_cell_.getLattice();
  if (Lattice.explicitly_defined)
  {
    if (Lattice.outOfBound(Lattice.toUnit(displ)))
      is_valid = false;
    else
    {
      SingleParticlePos newRedPos = Lattice.toUnit(active_pos_);
      if (!Lattice.isValid(newRedPos))
        is_valid = false;
    }
  }
  computeNewPosDistTablesAndSK(iat, active_pos_, true);
  return is_valid;
}

bool ParticleSet::makeMoveAndCheckWithSpin(Index_t iat, const SingleParticlePos& displ, const Scalar_t& sdispl)
{
  bool is_valid = makeMoveAndCheck(iat, displ);
  active_spin_val_ += sdispl;
  return is_valid;
}

void ParticleSet::computeNewPosDistTablesAndSK(Index_t iat, const SingleParticlePos& newpos, bool maybe_accept)
{
  ScopedTimer compute_newpos_scope(myTimers[PS_newpos]);

  for (int i = 0; i < DistTables.size(); ++i)
    DistTables[i]->move(*this, newpos, iat, maybe_accept);
  //Do not change SK: 2007-05-18
  //Change SK only if DoUpdate is true: 2008-09-12
  if (structure_factor_ && structure_factor_->DoUpdate)
    structure_factor_->makeMove(iat, newpos);
}

void ParticleSet::mw_computeNewPosDistTablesAndSK(const RefVectorWithLeader<ParticleSet>& p_list,
                                                  Index_t iat,
                                                  const std::vector<SingleParticlePos>& new_positions,
                                                  bool maybe_accept)
{
  ParticleSet& p_leader = p_list.getLeader();
  ScopedTimer compute_newpos_scope(p_leader.myTimers[PS_newpos]);

  {
    ScopedTimer copy_scope(p_leader.myTimers[PS_mw_copy]);
    const auto coords_list(extractCoordsRefList(p_list));
    p_leader.coordinates_->mw_copyActivePos(coords_list, iat, new_positions);
  }

  {
    ScopedTimer dt_scope(p_leader.myTimers[PS_dt_move]);
    const int dist_tables_size = p_leader.DistTables.size();
    for (int i = 0; i < dist_tables_size; ++i)
    {
      const auto dt_list(extractDTRefList(p_list, i));
      p_leader.DistTables[i]->mw_move(dt_list, p_list, new_positions, iat, maybe_accept);
    }

    // DistTables mw_move calls are asynchronous. Wait for them before return.
    PRAGMA_OFFLOAD("omp taskwait")
  }
  auto& SK = p_leader.structure_factor_;
  if (SK && SK->DoUpdate)
  {
#pragma omp parallel for
    for (int iw = 0; iw < p_list.size(); iw++)
      p_list[iw].structure_factor_->makeMove(iat, new_positions[iw]);
  }
}


bool ParticleSet::makeMoveAllParticles(const Walker_t& awalker, const ParticlePos& deltaR, RealType dt)
{
  active_ptcl_  = -1;
  auto& Lattice = simulation_cell_.getLattice();
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos displ(dt * deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos newpos(awalker.R[iat] + displ);
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
  if (structure_factor_)
    structure_factor_->updateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMoveAllParticles(const Walker_t& awalker,
                                       const ParticlePos& deltaR,
                                       const std::vector<RealType>& dt)
{
  active_ptcl_  = -1;
  auto& Lattice = simulation_cell_.getLattice();
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos displ(dt[iat] * deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos newpos(awalker.R[iat] + displ);
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
  if (structure_factor_)
    structure_factor_->updateAllPart(*this);
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
                                                const ParticlePos& drift,
                                                const ParticlePos& deltaR,
                                                RealType dt)
{
  active_ptcl_  = -1;
  auto& Lattice = simulation_cell_.getLattice();
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos displ(dt * deltaR[iat] + drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos newpos(awalker.R[iat] + displ);
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
  if (structure_factor_)
    structure_factor_->updateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMoveAllParticlesWithDrift(const Walker_t& awalker,
                                                const ParticlePos& drift,
                                                const ParticlePos& deltaR,
                                                const std::vector<RealType>& dt)
{
  active_ptcl_  = -1;
  auto& Lattice = simulation_cell_.getLattice();
  if (Lattice.explicitly_defined)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
    {
      SingleParticlePos displ(dt[iat] * deltaR[iat] + drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos newpos(awalker.R[iat] + displ);
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
  if (structure_factor_)
    structure_factor_->updateAllPart(*this);
  //every move is valid
  return true;
}

/** update the particle attribute by the proposed move
 *
 * When the active_ptcl_ is equal to iat, overwrite the position and update the
 * content of the distance tables.
 */
void ParticleSet::acceptMove(Index_t iat)
{
#ifndef NDEBUG
  if (iat != active_ptcl_)
    throw std::runtime_error("Bug detected by acceptMove! Request electron is not active!");
#endif
  ScopedTimer update_scope(myTimers[PS_accept]);
  //Update position + distance-table
  coordinates_->setOneParticlePos(active_pos_, iat);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->update(iat);

  //Do not change SK: 2007-05-18
  if (structure_factor_ && structure_factor_->DoUpdate)
    structure_factor_->acceptMove(iat, GroupID[iat], R[iat]);

  R[iat]       = active_pos_;
  spins[iat]   = active_spin_val_;
  active_ptcl_ = -1;
}

void ParticleSet::acceptMoveForwardMode(Index_t iat)
{
  assert(iat == active_ptcl_);
  ScopedTimer update_scope(myTimers[PS_accept]);
  //Update position + distance-table
  coordinates_->setOneParticlePos(active_pos_, iat);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->updatePartial(iat, true);

  //Do not change SK: 2007-05-18
  if (structure_factor_ && structure_factor_->DoUpdate)
    structure_factor_->acceptMove(iat, GroupID[iat], R[iat]);

  R[iat]       = active_pos_;
  spins[iat]   = active_spin_val_;
  active_ptcl_ = -1;
}

void ParticleSet::accept_rejectMove(Index_t iat, bool accepted, bool forward_mode)
{
  if (forward_mode)
    if (accepted)
      acceptMoveForwardMode(iat);
    else
      rejectMoveForwardMode(iat);
  else if (accepted)
    acceptMove(iat);
  else
    rejectMove(iat);
}

void ParticleSet::rejectMove(Index_t iat)
{
#ifndef NDEBUG
  if (iat != active_ptcl_)
    throw std::runtime_error("Bug detected by rejectMove! Request electron is not active!");
#endif
  active_ptcl_ = -1;
}

void ParticleSet::rejectMoveForwardMode(Index_t iat)
{
  assert(iat == active_ptcl_);
  //Update distance-table
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->updatePartial(iat, false);
  active_ptcl_ = -1;
}

void ParticleSet::mw_accept_rejectMove(const RefVectorWithLeader<ParticleSet>& p_list,
                                       Index_t iat,
                                       const std::vector<bool>& isAccepted,
                                       bool forward_mode)
{
  if (forward_mode)
  {
    ParticleSet& p_leader = p_list.getLeader();
    ScopedTimer update_scope(p_leader.myTimers[PS_accept]);

    const auto coords_list(extractCoordsRefList(p_list));
    std::vector<SingleParticlePos> new_positions;
    new_positions.reserve(p_list.size());
    for (const ParticleSet& pset : p_list)
      new_positions.push_back(pset.active_pos_);
    p_leader.coordinates_->mw_acceptParticlePos(coords_list, iat, new_positions, isAccepted);

    auto& dts = p_leader.DistTables;
    for (int i = 0; i < dts.size(); ++i)
    {
      const auto dt_list(extractDTRefList(p_list, i));
      dts[i]->mw_updatePartial(dt_list, iat, isAccepted);
    }

    for (int iw = 0; iw < p_list.size(); iw++)
    {
      assert(iat == p_list[iw].active_ptcl_);
      if (isAccepted[iw])
      {
        //Do not change SK: 2007-05-18
        auto& SK = p_leader.structure_factor_;
        if (SK && SK->DoUpdate)
          p_list[iw].structure_factor_->acceptMove(iat, p_list[iw].GroupID[iat], p_list[iw].R[iat]);

        p_list[iw].R[iat]     = p_list[iw].active_pos_;
        p_list[iw].spins[iat] = p_list[iw].active_spin_val_;
      }
      p_list[iw].active_ptcl_ = -1;
      assert(p_list[iw].R[iat] == p_list[iw].coordinates_->getAllParticlePos()[iat]);
    }
  }
  else
  {
    // loop over single walker acceptMove/rejectMove doesn't work safely.
    // need to code carefully for both coordinate and distance table updates
    // disable non-forward mode cases
    if (!forward_mode)
      throw std::runtime_error("BUG calling mw_accept_rejectMove in non-forward mode");
  }
}

void ParticleSet::donePbyP(bool skipSK)
{
  ScopedTimer donePbyP_scope(myTimers[PS_donePbyP]);
  coordinates_->donePbyP();
  if (!skipSK && structure_factor_ && !structure_factor_->DoUpdate)
    structure_factor_->updateAllPart(*this);
  for (size_t i = 0; i < DistTables.size(); ++i)
    DistTables[i]->finalizePbyP(*this);
  active_ptcl_ = -1;
}

void ParticleSet::mw_donePbyP(const RefVectorWithLeader<ParticleSet>& p_list, bool skipSK)
{
  ParticleSet& p_leader = p_list.getLeader();
  ScopedTimer donePbyP_scope(p_leader.myTimers[PS_donePbyP]);

  for (ParticleSet& pset : p_list)
  {
    pset.coordinates_->donePbyP();
    pset.active_ptcl_ = -1;
  }

  if (!skipSK && p_leader.structure_factor_ && !p_leader.structure_factor_->DoUpdate)
  {
    auto sk_list = extractSKRefList(p_list);
    StructFact::mw_updateAllPart(sk_list, p_list);
  }

  auto& dts = p_leader.DistTables;
  for (int i = 0; i < dts.size(); ++i)
  {
    const auto dt_list(extractDTRefList(p_list, i));
    dts[i]->mw_finalizePbyP(dt_list, p_list);
  }
}

void ParticleSet::makeVirtualMoves(const SingleParticlePos& newpos)
{
  active_ptcl_ = -1;
  active_pos_  = newpos;
  for (size_t i = 0; i < DistTables.size(); ++i)
    DistTables[i]->move(*this, newpos, active_ptcl_, false);
}

void ParticleSet::loadWalker(Walker_t& awalker, bool pbyp)
{
  ScopedTimer update_scope(myTimers[PS_loadWalker]);
  R     = awalker.R;
  spins = awalker.spins;
  coordinates_->setAllParticlePos(R);
#if !defined(SOA_MEMORY_OPTIMIZED)
  G = awalker.G;
  L = awalker.L;
#endif
  if (pbyp)
  {
    // in certain cases, full tables must be ready
    for (int i = 0; i < DistTables.size(); i++)
      if (DistTables[i]->getModes() & DTModes::NEED_FULL_TABLE_ANYTIME)
        DistTables[i]->evaluate(*this);
    //computed so that other objects can use them, e.g., kSpaceJastrow
    if (structure_factor_ && structure_factor_->DoUpdate)
      structure_factor_->updateAllPart(*this);
  }

  active_ptcl_ = -1;
}

void ParticleSet::mw_loadWalker(const RefVectorWithLeader<ParticleSet>& p_list,
                                const RefVector<Walker_t>& walkers,
                                const std::vector<bool>& recompute,
                                bool pbyp)
{
  auto& p_leader = p_list.getLeader();
  ScopedTimer load_scope(p_leader.myTimers[PS_loadWalker]);

  auto loadWalkerConfig = [](ParticleSet& pset, Walker_t& awalker) {
    pset.R     = awalker.R;
    pset.spins = awalker.spins;
    pset.coordinates_->setAllParticlePos(pset.R);
  };
#pragma omp parallel for
  for (int iw = 0; iw < p_list.size(); ++iw)
    if (recompute[iw])
      loadWalkerConfig(p_list[iw], walkers[iw]);

  if (pbyp)
  {
    auto& dts = p_leader.DistTables;
    for (int i = 0; i < dts.size(); ++i)
    {
      const auto dt_list(extractDTRefList(p_list, i));
      dts[i]->mw_recompute(dt_list, p_list, recompute);
    }

    if (p_leader.structure_factor_ && p_leader.structure_factor_->DoUpdate)
    {
#pragma omp parallel for
      for (int iw = 0; iw < p_list.size(); iw++)
        p_list[iw].structure_factor_->updateAllPart(p_list[iw]);
    }
  }
}

void ParticleSet::saveWalker(Walker_t& awalker)
{
  awalker.R     = R;
  awalker.spins = spins;
#if !defined(SOA_MEMORY_OPTIMIZED)
  awalker.G = G;
  awalker.L = L;
#endif
}

void ParticleSet::mw_saveWalker(const RefVectorWithLeader<ParticleSet>& psets, const RefVector<Walker_t>& walkers)
{
  for (int iw = 0; iw < psets.size(); ++iw)
    psets[iw].saveWalker(walkers[iw]);
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
  //   throw std::runtime_error("ParticleSet::initPropertyList");
  // }
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


void ParticleSet::createResource(ResourceCollection& collection) const
{
  coordinates_->createResource(collection);
  for (int i = 0; i < DistTables.size(); i++)
    DistTables[i]->createResource(collection);
}

void ParticleSet::acquireResource(ResourceCollection& collection, const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& ps_leader = p_list.getLeader();
  ps_leader.coordinates_->acquireResource(collection, extractCoordsRefList(p_list));
  for (int i = 0; i < ps_leader.DistTables.size(); i++)
    ps_leader.DistTables[i]->acquireResource(collection, extractDTRefList(p_list, i));
}

void ParticleSet::releaseResource(ResourceCollection& collection, const RefVectorWithLeader<ParticleSet>& p_list)
{
  auto& ps_leader = p_list.getLeader();
  ps_leader.coordinates_->releaseResource(collection, extractCoordsRefList(p_list));
  for (int i = 0; i < ps_leader.DistTables.size(); i++)
    ps_leader.DistTables[i]->releaseResource(collection, extractDTRefList(p_list, i));
}

RefVectorWithLeader<DistanceTable> ParticleSet::extractDTRefList(const RefVectorWithLeader<ParticleSet>& p_list, int id)
{
  RefVectorWithLeader<DistanceTable> dt_list(*p_list.getLeader().DistTables[id]);
  dt_list.reserve(p_list.size());
  for (ParticleSet& p : p_list)
    dt_list.push_back(*p.DistTables[id]);
  return dt_list;
}

RefVectorWithLeader<DynamicCoordinates> ParticleSet::extractCoordsRefList(
    const RefVectorWithLeader<ParticleSet>& p_list)
{
  RefVectorWithLeader<DynamicCoordinates> coords_list(*p_list.getLeader().coordinates_);
  coords_list.reserve(p_list.size());
  for (ParticleSet& p : p_list)
    coords_list.push_back(*p.coordinates_);
  return coords_list;
}

RefVectorWithLeader<StructFact> ParticleSet::extractSKRefList(const RefVectorWithLeader<ParticleSet>& p_list)
{
  RefVectorWithLeader<StructFact> sk_list(*p_list.getLeader().structure_factor_);
  sk_list.reserve(p_list.size());
  for (ParticleSet& p : p_list)
    sk_list.push_back(*p.structure_factor_);
  return sk_list;
}

} // namespace qmcplusplus
