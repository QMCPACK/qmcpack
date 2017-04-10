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
#include "Particle/DistanceTable.h"
#include "LongRange/StructFact.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"

//#define PACK_DISTANCETABLES

namespace qmcplusplus
{

//using namespace particle_info;

#ifdef QMC_CUDA
template<> int ParticleSet::Walker_t::cuda_DataSize = 0;
#endif

///object counter
int  ParticleSet::PtclObjectCounter = 0;

void add_p_timer(std::vector<NewTimer*>& timers)
{
  timers.push_back(new NewTimer("ParticleSet::makeMove",timer_level_fine)); //timer for MC, ratio etc
  timers.push_back(new NewTimer("ParticleSet::makeMoveOnSphere",timer_level_fine)); //timer for the walker loop
  TimerManager.addTimer(timers[0]);
  TimerManager.addTimer(timers[1]);
}

ParticleSet::ParticleSet()
  : UseBoundBox(true), UseSphereUpdate(true), IsGrouped(true)
  , ThreadID(0), SK(0), ParentTag(-1), ParentName("0")
  , quantum_domain(classical)
{
  initParticleSet();
  initPropertyList();
  add_p_timer(myTimers);
}

ParticleSet::ParticleSet(const ParticleSet& p)
  : UseBoundBox(p.UseBoundBox), UseSphereUpdate(p.UseSphereUpdate),IsGrouped(p.IsGrouped)
  , ThreadID(0), mySpecies(p.getSpeciesSet()),SK(0), ParentTag(p.tag()), ParentName(p.parentName())
{
  set_quantum_domain(p.quantum_domain);
  initBase();
  initParticleSet();
  assign(p); //only the base is copied, assumes that other properties are not assignable
  //need explicit copy:
  Mass=p.Mass;
  Z=p.Z;
  std::ostringstream o;
  o<<p.getName()<<ObjectTag;
  this->setName(o.str());
  app_log() << "  Copying a particle set " << p.getName() << " to " << this->getName() << " groups=" << groups() << std::endl;
  PropertyList.Names=p.PropertyList.Names;
  PropertyList.Values=p.PropertyList.Values;
  PropertyHistory=p.PropertyHistory;
  Collectables=p.Collectables;
  //construct the distance tables with the same order
  if(p.DistTables.size())
  {
    app_log() << "  Cloning distance tables. It has " << p.DistTables.size() << std::endl;
    addTable(*this,p.DistTables[0]->DTType); //first is always for this-this paier
    for (int i=1; i<p.DistTables.size(); ++i)
      addTable(p.DistTables[i]->origin(),p.DistTables[i]->DTType);
  }
  if(p.SK)
  {
    LRBox=p.LRBox; //copy LRBox
    SK=new StructFact(*p.SK); //safe to use the copy constructor
    //R.InUnit=p.R.InUnit;
    //createSK();
    //SK->DoUpdate=p.SK->DoUpdate;
  }
  if (p.Sphere.size())
    resizeSphere(p.Sphere.size());
  add_p_timer(myTimers);
  myTwist=p.myTwist;

  RSoA.resize(getLocalNum());
}

ParticleSet::~ParticleSet()
{
  DEBUG_MEMORY("ParticleSet::~ParticleSet");
  delete_iter(DistTables.begin(), DistTables.end());
  if (SK)
    delete SK;
  delete_iter(Sphere.begin(), Sphere.end());
}

void ParticleSet::create(unsigned n)
{
  createBase(n);
  RSoA.resize(n);
}

void ParticleSet::create(const std::vector<int>& agroup)
{
  createBase(agroup);
  RSoA.resize(getTotalNum());
}

void ParticleSet::set_quantum_domain(quantum_domains qdomain)
{
  if(quantum_domain_valid(qdomain))
    quantum_domain = qdomain;
  else
    APP_ABORT("ParticleSet::set_quantum_domain\n  input quantum domain is not valid for particles");
}

void ParticleSet::initParticleSet()
{
  #pragma omp critical (PtclObjectCounter)
  {
    ObjectTag = PtclObjectCounter;
    PtclObjectCounter++;
  }

  G.setTypeName(ParticleTags::gradtype_tag);
  L.setTypeName(ParticleTags::laptype_tag);
  dG.setTypeName(ParticleTags::gradtype_tag);
  dL.setTypeName(ParticleTags::laptype_tag);

  G.setObjName("grad");
  L.setObjName("lap");
  dG.setObjName("dgrad");
  dL.setObjName("dlap");

  addAttribute(G);
  addAttribute(L);
  addAttribute(dG);
  addAttribute(dL);

  //more particle attributes
  Mass.setTypeName(ParticleTags::scalartype_tag);
  Mass.setObjName("mass");
  SameMass=true; //default is SameMass
  addAttribute(Mass);

  Z.setTypeName(ParticleTags::scalartype_tag);
  Z.setObjName("charge");
  addAttribute(Z);

  PCID.setTypeName(ParticleTags::indextype_tag); //add PCID tags
  PCID.setObjName("pcid");
  addAttribute(PCID);

  IndirectID.setTypeName(ParticleTags::indextype_tag); //add IndirectID tags
  IndirectID.setObjName("id1");
  addAttribute(IndirectID);

  myTwist=0.0;

  activeWalker=nullptr;
}

void ParticleSet::resetGroups()
{
  int nspecies=mySpecies.getTotalNum();
  if(nspecies==0)
  {
    APP_ABORT("ParticleSet::resetGroups() Failed. No species exisits");
  }
  int natt=mySpecies.numAttributes();
  int qind=mySpecies.addAttribute("charge");
  if(natt==qind)
  {
    app_log() << " Missing charge attribute of the SpeciesSet " << myName << " particleset" << std::endl;
    app_log() << " Assume neutral particles Z=0.0 " << std::endl;
    for(int ig=0; ig<nspecies; ig++)
      mySpecies(qind,ig)=0.0;
  }
  for(int iat=0; iat<Z.size(); iat++)
    Z[iat]=mySpecies(qind,GroupID[iat]);
  natt=mySpecies.numAttributes();
  int massind=mySpecies.addAttribute("mass");
  if(massind==natt)
  {
    for(int ig=0; ig<nspecies; ig++)
      mySpecies(massind,ig)=1.0;
  }
  SameMass=true;
  double m0=mySpecies(massind,0);
  for(int ig=1; ig<nspecies; ig++)
    SameMass &= (mySpecies(massind,ig)== m0);
  if(SameMass)
    app_log() << "  All the species have the same mass " << m0 << std::endl;
  else
    app_log() << "  Distinctive masses for each species " << std::endl;
  for(int iat=0; iat<Mass.size(); iat++)
    Mass[iat]=mySpecies(massind,GroupID[iat]);
  std::vector<int> ng(nspecies,0);
  for(int iat=0; iat<GroupID.size(); iat++)
  {
    if(GroupID[iat]<nspecies)
      ng[GroupID[iat]]++;
    else
      APP_ABORT("ParticleSet::resetGroups() Failed. GroupID is out of bound.");
  }
  SubPtcl.resize(nspecies+1);
  SubPtcl[0]=0;
  for(int i=0; i<nspecies; ++i)
    SubPtcl[i+1]=SubPtcl[i]+ng[i];
  int membersize= mySpecies.addAttribute("membersize");
  for(int ig=0; ig<nspecies; ++ig)
    mySpecies(membersize,ig)=ng[ig];
  //orgID=ID;
  //orgGroupID=GroupID;
  int new_id=0;
  for(int i=0; i<nspecies; ++i)
    for(int iat=0; iat<GroupID.size(); ++iat)
      if(GroupID[iat]==i)
        IndirectID[new_id++]=ID[iat];
  IsGrouped=true;
  for(int iat=0; iat<ID.size(); ++iat)
    IsGrouped &= (IndirectID[iat]==ID[iat]);
  if(IsGrouped)
    app_log() << "Particles are grouped. Safe to use groups " << std::endl;
  else
    app_log() << "ID is not grouped. Need to use IndirectID for species-dependent operations " << std::endl;
}

void
ParticleSet::randomizeFromSource (ParticleSet &src)
{
  SpeciesSet& srcSpSet(src.getSpeciesSet());
  SpeciesSet& spSet(getSpeciesSet());
  int srcChargeIndx = srcSpSet.addAttribute("charge");
  int srcMemberIndx = srcSpSet.addAttribute("membersize");
  int ChargeIndex   = spSet.addAttribute("charge");
  int MemberIndx    = spSet.addAttribute("membersize");
  int Nsrc  = src.getTotalNum();
  int Nptcl = getTotalNum();
  int NumSpecies    = spSet.TotalNum;
  int NumSrcSpecies = srcSpSet.TotalNum;
  //Store information about charges and number of each species
  std::vector<int> Zat, Zspec, NofSpecies, NofSrcSpecies, CurElec;
  Zat.resize(Nsrc);
  Zspec.resize(NumSrcSpecies);
  NofSpecies.resize(NumSpecies);
  CurElec.resize(NumSpecies);
  NofSrcSpecies.resize(NumSrcSpecies);
  for(int spec=0; spec<NumSrcSpecies; spec++)
  {
    Zspec[spec] = (int)round(srcSpSet(srcChargeIndx,spec));
    NofSrcSpecies[spec] = (int)round(srcSpSet(srcMemberIndx,spec));
  }
  for(int spec=0; spec<NumSpecies; spec++)
  {
    NofSpecies[spec] = (int)round(spSet(MemberIndx,spec));
    CurElec[spec] = first(spec);
  }
  int totQ=0;
  for(int iat=0; iat<Nsrc; iat++)
    totQ+=Zat[iat] = Zspec[src.GroupID[iat]];
  app_log() << "  Total ion charge    = " << totQ << std::endl;
  totQ -= Nptcl;
  app_log() << "  Total system charge = " << totQ << std::endl;
  // Now, loop over ions, attaching electrons to them to neutralize
  // charge
  int spToken = 0;
  // This is decremented when we run out of electrons in each species
  int spLeft = NumSpecies;
  std::vector<PosType> gaussRand (Nptcl);
  makeGaussRandom (gaussRand);
  for (int iat=0; iat<Nsrc; iat++)
  {
    // Loop over electrons to add, selecting round-robin from the
    // electron species
    int z = Zat[iat];
    while (z > 0  && spLeft)
    {
      int sp = spToken++ % NumSpecies;
      if (NofSpecies[sp])
      {
        NofSpecies[sp]--;
        z--;
        int elec = CurElec[sp]++;
        app_log() << "  Assigning " << (sp ? "down" : "up  ")
                  << " electron " << elec << " to ion " << iat
                  << " with charge " << z << std::endl;
        double radius = 0.5* std::sqrt((double)Zat[iat]);
        R[elec] = src.R[iat] + radius * gaussRand[elec];
      }
      else
        spLeft--;
    }
  }
  // Assign remaining electrons
  int ion=0;
  for (int sp=0; sp < NumSpecies; sp++)
  {
    for (int ie=0; ie<NofSpecies[sp]; ie++)
    {
      int iat = ion++ % Nsrc;
      double radius = std::sqrt((double)Zat[iat]);
      int elec = CurElec[sp]++;
      R[elec] = src.R[iat] + radius * gaussRand[elec];
    }
  }
}

///write to a std::ostream
bool ParticleSet::get(std::ostream& os) const
{
  os << "  ParticleSet " << getName() << " : ";
  for (int i=0; i<SubPtcl.size(); i++)
    os << SubPtcl[i] << " ";
  os <<"\n\n    " << LocalNum << "\n\n";
  const int maxParticlesToPrint = 10;
  int numToPrint = std::min(LocalNum, maxParticlesToPrint);

  for (int i=0; i<numToPrint; i++)
  {
    os << "    " << mySpecies.speciesName[GroupID[i]]  << R[i] << std::endl;
  }

  if (numToPrint < LocalNum)
  {
    os << "    (... and " << (LocalNum-numToPrint) << " more particle positions ...)" << std::endl;
  }

  return true;
}

///read from std::istream
bool ParticleSet::put( std::istream& is)
{
  return true;
}

///reset member data
void ParticleSet::reset()
{
  app_log() << "<<<< going to set properties >>>> " << std::endl;
}

///read the particleset
bool ParticleSet::put(xmlNodePtr cur)
{
  return true;
}

void ParticleSet::setBoundBox(bool yes)
{
  UseBoundBox=yes;
}

void ParticleSet::checkBoundBox(RealType rb)
{
  if (UseBoundBox && rb>Lattice.SimulationCellRadius)
  {
    app_warning()
        << "ParticleSet::checkBoundBox "
        << rb << "> SimulationCellRadius=" << Lattice.SimulationCellRadius
        << "\n Using SLOW method for the sphere update. " << std::endl;
    UseSphereUpdate=false;
  }
}
//void ParticleSet::setUpdateMode(int updatemode) {
//  if(DistTables.empty()) {
//    DistanceTable::getTables(ObjectTag,DistTables);
//    DistanceTable::create(1);
//    LOGMSG("ParticleSet::setUpdateMode to create distance tables.")
//    LOGMSG("\t the number of distance tables for " << getName() << " " << DistTables.size())
//  }
//}

///** add a distance table to DistTables list
// * @param d_table pointer to a DistanceTableData to be added
// *
// * DistTables is a list of DistanceTables which are updated by MC moves.
// */
//void ParticleSet::addTable(DistanceTableData* d_table) {
//  int oid=d_table->origin().tag();
//  int i=0;
//  int dsize=DistTables.size();
//  while(i<dsize) {
//    if(oid == DistTables[i]->origin().tag()) //table already exists
//      return;
//    ++i;
//  }
//  DistTables.push_back(d_table);
//}
int ParticleSet::addTable(const ParticleSet& psrc, int dt_type)
{
  if (DistTables.empty())
  {
    DistTables.reserve(4);
#if defined(ENABLE_AA_SOA)
    DistTables.push_back(createDistanceTable(*this,DT_SOA));
#else
    //if(dt_type==DT_SOA_PREFERRED) dt_type=DT_AOS; //safety
    DistTables.push_back(createDistanceTable(*this,dt_type));
#endif
    //add  this-this pair
    myDistTableMap.clear();
    myDistTableMap[ObjectTag]=0;
    app_log() << "  ... ParticleSet::addTable Create Table #0 " << DistTables[0]->Name << std::endl;
    DistTables[0]->ID=0;
    if (psrc.tag() == ObjectTag)
      return 0;
  }
  if (psrc.tag() == ObjectTag)
  {
    app_log() << "  ... ParticleSet::addTable Reuse Table #" << 0 << " " << DistTables[0]->Name << std::endl;
    //if(!DistTables[0]->is_same_type(dt_type))
    //{//itself is special, cannot mix them: some of the users do not check the index
    //  APP_ABORT("ParticleSet::addTable for itself Cannot mix AoS and SoA distance tables.\n");
    //}
    return 0;
  }
  int tsize=DistTables.size(),tid;
  std::map<int,int>::iterator tit(myDistTableMap.find(psrc.tag()));
  if (tit == myDistTableMap.end())
  {
    tid=DistTables.size();
    DistTables.push_back(createDistanceTable(psrc,*this,dt_type));
    myDistTableMap[psrc.tag()]=tid;
    DistTables[tid]->ID=tid;
    app_log() << "  ... ParticleSet::addTable Create Table #" << tid << " " << DistTables[tid]->Name << std::endl;
  }
  else
  {
    tid = (*tit).second;
    if(dt_type == DT_SOA_PREFERRED || DistTables[tid]->is_same_type(dt_type))  //good to reuse
    {
      app_log() << "  ... ParticleSet::addTable Reuse Table #" << tid << " " << DistTables[tid]->Name << std::endl;
    }
    else
    {
      APP_ABORT("ParticleSet::addTable Cannot mix AoS and SoA distance tables.\n");
    }
    //if(dt_type == DT_SOA || dt_type == DT_AOS) //not compatible
    //{
    //}
    //for DT_SOA_PREFERRED or DT_AOS_PREFERRED, return the existing table
  }
  app_log().flush();
  return tid;
}

int ParticleSet::getTable(const ParticleSet& psrc)
{
  int tid;
  if (DistTables.empty())
    tid = -1;
  else
    if (psrc.tag() == ObjectTag)
      tid = 0;
    else
    {
      std::map<int,int>::iterator tit(myDistTableMap.find(psrc.tag()));
      if (tit == myDistTableMap.end())
        tid = -1;
      else
        tid = (*tit).second;
    }
  return tid;
}

void ParticleSet::update(int iflag)
{
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
  for (int i=0; i< DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);

  Ready4Measure=true;
}

void ParticleSet::update(const ParticlePos_t& pos)
{
  R = pos;
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
  for (int i=0; i< DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK && !SK->DoUpdate)
    SK->UpdateAllPart(*this);

  Ready4Measure=true;
}

/** move a particle iat
 * @param iat the index of the particle to be moved
 * @param displ the displacement of the iath-particle position
 * @return the proposed position
 *
 * Update activePtcl index and activePos position for the proposed move.
 * Evaluate the related distance table data DistanceTableData::Temp.
 */
ParticleSet::SingleParticlePos_t
ParticleSet::makeMove(Index_t iat, const SingleParticlePos_t& displ)
{
  activePtcl=iat;
  activePos=R[iat]; //save the current position
  SingleParticlePos_t newpos(activePos+displ);
  for (int i=0; i< DistTables.size(); ++i)
    DistTables[i]->move(*this,newpos,iat);
  R[iat]=newpos;
  //Do not change SK: 2007-05-18
  //Change SK only if DoUpdate is true: 2008-09-12
  if (SK && SK->DoUpdate)
    SK->makeMove(iat,newpos);
  return newpos;
}

void ParticleSet::setActive(int iat)
{
  for (size_t i=0,n=DistTables.size(); i< n; i++)
    DistTables[i]->evaluate(*this,iat);
}


/** move a particle iat
 * @param iat the index of the particle to be moved
 * @param displ the displacement of the iath-particle position
 * @return the proposed position
 *
 * Update activePtcl index and activePos position for the proposed move.
 * Evaluate the related distance table data DistanceTableData::Temp.
 */
bool
ParticleSet::makeMoveAndCheck(Index_t iat, const SingleParticlePos_t& displ)
{
  myTimers[0]->start();
  activePtcl=iat;
  //SingleParticlePos_t red_displ(Lattice.toUnit(displ));
  if (UseBoundBox)
  {
    if (Lattice.outOfBound(Lattice.toUnit(displ)))
    {
      myTimers[0]->stop();
      return false;
    }
    activePos=R[iat]; //save the current position
    SingleParticlePos_t newpos(activePos+displ);
    newRedPos=Lattice.toUnit(newpos);
    if (Lattice.isValid(newRedPos))
    {
      for (int i=0; i< DistTables.size(); ++i)
        DistTables[i]->move(*this,newpos,iat);
      R[iat]=newpos;
      if (SK && SK->DoUpdate)
        SK->makeMove(iat,newpos);
      myTimers[0]->stop();
      return true;
    }
    //out of bound
    myTimers[0]->stop();
    return false;
  }
  else
  {
    activePos=R[iat]; //save the current position
    SingleParticlePos_t newpos(activePos+displ);
    for (int i=0; i< DistTables.size(); ++i)
      DistTables[i]->move(*this,newpos,iat);
    R[iat]=newpos;
    myTimers[0]->stop();
    return true;
  }
}

bool ParticleSet::makeMove(const Walker_t& awalker
                           , const ParticlePos_t& deltaR, RealType dt)
{
  activePtcl=-1;
  if (UseBoundBox)
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt*deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat]+displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat]=newpos;
    }
  }
  else
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
      R[iat]=awalker.R[iat]+dt*deltaR[iat];
  }
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
  for (int i=0; i< DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMove(const Walker_t& awalker
                           , const ParticlePos_t& deltaR, const std::vector<RealType>& dt)
{
  Ready4Measure=false;
  activePtcl=-1;
  if (UseBoundBox)
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt[iat]*deltaR[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat]+displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat]=newpos;
    }
  }
  else
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
      R[iat]=awalker.R[iat]+dt[iat]*deltaR[iat];
  }
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
  for (int i=0; i< DistTables.size(); i++)
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
bool ParticleSet::makeMoveWithDrift(const Walker_t& awalker
                                    , const ParticlePos_t& drift , const ParticlePos_t& deltaR
                                    , RealType dt)
{
  Ready4Measure=false;
  activePtcl=-1;
  if (UseBoundBox)
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt*deltaR[iat]+drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat]+displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat]=newpos;
    }
  }
  else
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
      R[iat]=awalker.R[iat]+dt*deltaR[iat]+drift[iat];
  }
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
  for (int i=0; i< DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}

bool ParticleSet::makeMoveWithDrift(const Walker_t& awalker
                                    , const ParticlePos_t& drift , const ParticlePos_t& deltaR
                                    , const std::vector<RealType>& dt)
{
  Ready4Measure=false;
  activePtcl=-1;
  if (UseBoundBox)
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
    {
      SingleParticlePos_t displ(dt[iat]*deltaR[iat]+drift[iat]);
      if (Lattice.outOfBound(Lattice.toUnit(displ)))
        return false;
      SingleParticlePos_t newpos(awalker.R[iat]+displ);
      if (!Lattice.isValid(Lattice.toUnit(newpos)))
        return false;
      R[iat]=newpos;
    }
  }
  else
  {
    for (int iat=0; iat<deltaR.size(); ++iat)
      R[iat]=awalker.R[iat]+dt[iat]*deltaR[iat]+drift[iat];
  }

#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif

  for (int i=0; i< DistTables.size(); i++)
    DistTables[i]->evaluate(*this);
  if (SK)
    SK->UpdateAllPart(*this);
  //every move is valid
  return true;
}


/** move the iat-th particle by displ
 *
 * @param iat the particle that is moved on a sphere
 * @param displ displacement from the current position
 */
void
ParticleSet::makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ)
{
  activePtcl=iat;
  activePos=R[iat]; //save the current position
  SingleParticlePos_t newpos(activePos+displ);
  for (int i=0; i< DistTables.size(); ++i)
    DistTables[i]->moveOnSphere(*this,newpos,iat);
  R[iat]=newpos;
  if (SK && SK->DoUpdate)
    SK->makeMove(iat,R[iat]);
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
    for (int i=0,n=DistTables.size(); i< n; i++)
      DistTables[i]->update(iat);

    if(RSoA.size() != getLocalNum())
      std::cout << "Die here " << RSoA.size() << std::endl;

    RSoA(iat)=R[iat];

    //Do not change SK: 2007-05-18
    if (SK && SK->DoUpdate)
      SK->acceptMove(iat,GroupID[iat]);
  }
  else
  {
    std::ostringstream o;
    o << "  Illegal acceptMove " << iat << " != " << activePtcl;
    APP_ABORT(o.str());
  }
}

void ParticleSet::rejectMove(Index_t iat)
{
  //restore the position by the saved activePos
  R[iat]=activePos;
  for (int i=0; i< DistTables.size(); ++i)
    DistTables[i]->activePtcl=-1;
}

void ParticleSet::donePbyP()
{
  for (size_t i=0,nt=DistTables.size(); i< nt; i++)
    DistTables[i]->donePbyP();
  if (SK && !SK->DoUpdate)
    SK->UpdateAllPart(*this);
  Ready4Measure=true;
}

void ParticleSet::makeVirtualMoves(const SingleParticlePos_t& newpos)
{
  activePtcl=0;
  activePos=R[0];
  for (size_t i=0; i< DistTables.size(); ++i)
    DistTables[i]->move(*this,newpos,0);
  R[0]=newpos;
}


/** resize Sphere by the LocalNum
 * @param nc number of centers to which Spherical grid will be assigned.
 */
void ParticleSet::resizeSphere(int nc)
{
  int nsadd=nc-Sphere.size();
  while (nsadd>0)
  {
    Sphere.push_back(new ParticlePos_t);
    --nsadd;
  }
}

void ParticleSet::loadWalker(Walker_t& awalker, bool pbyp)
{
  R = awalker.R;
#if defined(ENABLE_AA_SOA)
  RSoA.copyIn(R); 
#endif
#if !defined(SOA_MEMORY_OPTIMIZED)
  G = awalker.G;
  L = awalker.L;
#endif
  if (pbyp)
  {
    //for (size_t i=0,nt=DistTables.sie(); i< nt; i++)
    //{
    //  if(DistTables[i]!= DT_SOA) DistTables[i]->evaluate(*this);
    //}
    //computed so that other objects can use them, e.g., kSpaceJastrow
    if(SK && SK->DoUpdate)
      SK->UpdateAllPart(*this);
  }

  Ready4Measure=false;
}

void ParticleSet::loadWalker(Walker_t* awalker)
{
  if(activeWalker != awalker)
  {
    activeWalker=awalker;
    R = awalker->R;
  }
}

void ParticleSet::saveWalker(Walker_t& awalker)
{
  awalker.R=R;
#if !defined(SOA_MEMORY_OPTIMIZED)
  awalker.G=G;
  awalker.L=L;
#endif
  //PAOps<RealType,OHMMS_DIM>::copy(G,awalker.Drift);
  //if (SK)
  //  SK->UpdateAllPart(*this);
  //awalker.DataSet.rewind();
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
  delete_iter(DistTables.begin(),DistTables.end());
  DistTables.clear();
  //for(int i=0; i< DistTables.size(); i++) DistanceTable::removeTable(DistTables[i]->getName());
  //DistTables.erase(DistTables.begin(),DistTables.end());
}

int ParticleSet::addPropertyHistory(int leng)
{
  int newL = PropertyHistory.size();
  std::vector<EstimatorRealType> newVecHistory=std::vector<EstimatorRealType>(leng,0.0);
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




}

