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

#ifdef QMC_CUDA
  template<> int ParticleSet::Walker_t::cuda_DataSize = 0;
#endif

  ///object counter
  int  ParticleSet::PtclObjectCounter = 0;

  void add_p_timer(vector<NewTimer*>& timers)
  {
    timers.push_back(new NewTimer("ParticleSet::makeMove")); //timer for MC, ratio etc
    timers.push_back(new NewTimer("ParticleSet::makeMoveOnSphere")); //timer for the walker loop
    TimerManager.addTimer(timers[0]);
    TimerManager.addTimer(timers[1]);
  }

  ParticleSet::ParticleSet()
    : UseBoundBox(true), UseSphereUpdate(true), sorted_ids(false), reordered_ids(false)
      , SK(0), ParentTag(-1)
    {
      initParticleSet();
      initPropertyList();

      add_p_timer(myTimers);
    }

  ParticleSet::ParticleSet(const ParticleSet& p)
      : UseBoundBox(p.UseBoundBox), UseSphereUpdate(p.UseSphereUpdate)
        ,sorted_ids(p.sorted_ids), reordered_ids(p.reordered_ids)
      , SK(0), mySpecies(p.getSpeciesSet()), ParentTag(p.tag())
  {
    initBase();
    initParticleSet();
    assign(p);

    ostringstream o;
    o<<p.getName()<<ObjectTag;
    this->setName(o.str());
    app_log() << "  Copying a particle set " << p.getName() << " to " << this->getName() << " groups=" << groups() << endl;

    PropertyList.Names=p.PropertyList.Names;
    PropertyList.Values=p.PropertyList.Values;

    PropertyHistory=  p.PropertyHistory;

    Collectables=p.Collectables;

    //construct the distance tables with the same order
    //first is always for this-this paier
    for (int i=1;i<p.DistTables.size(); ++i) addTable(p.DistTables[i]->origin());

    if(p.SK)
    {
      R.InUnit=p.R.InUnit;
      createSK();
      SK->DoUpdate=p.SK->DoUpdate;
    }

    if (p.Sphere.size()) resizeSphere(p.Sphere.size());

    add_p_timer(myTimers);
    myTwist=p.myTwist;
  }


  ParticleSet::~ParticleSet()
  {
    DEBUG_MEMORY("ParticleSet::~ParticleSet");
    delete_iter(DistTables.begin(), DistTables.end());
    if (SK) delete SK;
    delete_iter(Sphere.begin(), Sphere.end());
  }

  void ParticleSet::initParticleSet()
  {
    ObjectTag = PtclObjectCounter;
#pragma omp atomic
    PtclObjectCounter++;

#if defined(QMC_COMPLEX)
    G.setTypeName(ParticleTags::gradtype_tag);
    L.setTypeName(ParticleTags::laptype_tag);
    dG.setTypeName(ParticleTags::gradtype_tag);
    dL.setTypeName(ParticleTags::laptype_tag);
#else
    G.setTypeName(ParticleTags::postype_tag);
    L.setTypeName(ParticleTags::scalartype_tag);
    dG.setTypeName(ParticleTags::postype_tag);
    dL.setTypeName(ParticleTags::scalartype_tag);
#endif

    G.setObjName("grad");
    L.setObjName("lap");
    dG.setObjName("dgrad");
    dL.setObjName("dlap");
    addAttribute(G);
    addAttribute(L);
    addAttribute(dG);
    addAttribute(dL);

    orgID.setTypeName(ParticleTags::indextype_tag); //add ID tags
    orgID.setObjName("id0"); 
    addAttribute(orgID);
    
    orgGroupID.setTypeName(ParticleTags::indextype_tag); //add ID tags
    orgGroupID.setObjName("gid0"); 
    addAttribute(orgGroupID);

    myTwist=0.0;

    ////this has to be in unit
    //redR.setObjName("redpos");
    //redR.InUnit=PosUnit::LatticeUnit;
    //addAttribute(redR);
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
    vector<int> Zat, Zspec, NofSpecies, NofSrcSpecies, CurElec;
    Zat.resize(Nsrc); Zspec.resize(NumSrcSpecies);

    NofSpecies.resize(NumSpecies);
    CurElec.resize(NumSpecies);
    NofSrcSpecies.resize(NumSrcSpecies);
    

    for(int spec=0; spec<NumSrcSpecies; spec++) { 
      Zspec[spec] = (int)round(srcSpSet(srcChargeIndx,spec));
      NofSrcSpecies[spec] = (int)round(srcSpSet(srcMemberIndx,spec));
    }

    for(int spec=0; spec<NumSpecies; spec++) {
      NofSpecies[spec] = (int)round(spSet(MemberIndx,spec));
      CurElec[spec] = first(spec);
    }

    int totQ=0;
    for(int iat=0; iat<Nsrc; iat++) 
      totQ+=Zat[iat] = Zspec[src.GroupID[iat]];

    app_log() << "  Total ion charge    = " << totQ << endl;
    totQ -= Nptcl;
    app_log() << "  Total system charge = " << totQ << endl;

    // Now, loop over ions, attaching electrons to them to neutralize
    // charge 
    int spToken = 0;
    // This is decremented when we run out of electrons in each species
    int spLeft = NumSpecies;
    vector<PosType> gaussRand (Nptcl);
    makeGaussRandom (gaussRand);
    for (int iat=0; iat<Nsrc; iat++) {
      // Loop over electrons to add, selecting round-robin from the
      // electron species
      int z = Zat[iat];
      while (z > 0  && spLeft) {
	int sp = spToken++ % NumSpecies;
	if (NofSpecies[sp]) {
	  NofSpecies[sp]--;
	  z--;
	  int elec = CurElec[sp]++;
	  app_log() << "  Assigning " << (sp ? "down" : "up  ") 
		    << " electron " << elec << " to ion " << iat 
		    << " with charge " << z << endl;
	  double radius = 0.5* std::sqrt((double)Zat[iat]);
	  R[elec] = src.R[iat] + radius * gaussRand[elec];
	}
	else
	  spLeft--;
      }
    }
    // Assign remaining electrons
    int ion=0;
    for (int sp=0; sp < NumSpecies; sp++) {
      for (int ie=0; ie<NofSpecies[sp]; ie++) {
	int iat = ion++ % Nsrc;
	double radius = std::sqrt((double)Zat[iat]);
	int elec = CurElec[sp]++;
	R[elec] = src.R[iat] + radius * gaussRand[elec];
      }
    }

  }

  ///write to a ostream
  bool ParticleSet::get(ostream& os) const
    {

      os << "  ParticleSet " << getName() << " : ";
      for (int i=0; i<SubPtcl.size(); i++) os << SubPtcl[i] << " ";
      os <<"\n\n    " << LocalNum << "\n\n";

      for (int i=0; i<LocalNum; i++)
        {
          os << "    " << mySpecies.speciesName[GroupID[i]]  << R[i] << endl;
        }
      return true;
    }

  ///read from istream
  bool ParticleSet::put(istream& is)
  {
    return true;
  }

  ///reset member data
  void ParticleSet::reset()
  {
    app_log() << "<<<< going to set properties >>>> " << endl;
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
        << "\n Using SLOW method for the sphere update. " <<endl;
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
  int ParticleSet::addTable(const ParticleSet& psrc)
  {
    if (DistTables.empty())
      {
        DistTables.reserve(4);
        DistTables.push_back(createDistanceTable(*this));
        //add  this-this pair
        myDistTableMap.clear();
        myDistTableMap[ObjectTag]=0;
        app_log() << "  ... ParticleSet::addTable Create Table #0 " << DistTables[0]->Name << endl;

        if (psrc.tag() == ObjectTag) return 0;
      }

    if (psrc.tag() == ObjectTag)
      {
        app_log() << "  ... ParticleSet::addTable Reuse Table #" << 0 << " " << DistTables[0]->Name <<endl;
        return 0;
      }

    int tsize=DistTables.size(),tid;
    map<int,int>::iterator tit(myDistTableMap.find(psrc.tag()));
    if (tit == myDistTableMap.end())
      {
        tid=DistTables.size();
        DistTables.push_back(createDistanceTable(psrc,*this));
        myDistTableMap[psrc.tag()]=tid;
        app_log() << "  ... ParticleSet::addTable Create Table #" << tid << " " << DistTables[tid]->Name <<endl;
      }
    else
      {
        tid = (*tit).second;
        app_log() << "  ... ParticleSet::addTable Reuse Table #" << tid << " " << DistTables[tid]->Name << endl;
      }
    return tid;
  }

  void ParticleSet::update(int iflag)
  {
    for (int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if (SK) SK->UpdateAllPart();
  }

  void ParticleSet::update(const ParticlePos_t& pos)
  {
    R = pos;
    for (int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if (SK && !SK->DoUpdate) SK->UpdateAllPart();
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
    if (SK && SK->DoUpdate) SK->makeMove(iat,newpos);
    return newpos;
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
        if (Lattice.outOfBound(Lattice.toUnit(displ))) return false;
        activePos=R[iat]; //save the current position
        SingleParticlePos_t newpos(activePos+displ);
        newRedPos=Lattice.toUnit(newpos);
        if (Lattice.isValid(newRedPos))
          {
            for (int i=0; i< DistTables.size(); ++i)
              DistTables[i]->move(*this,newpos,iat);
            R[iat]=newpos;
            if (SK && SK->DoUpdate) SK->makeMove(iat,newpos);
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
    if (UseBoundBox)
      {
        for (int iat=0; iat<deltaR.size(); ++iat)
          {
            SingleParticlePos_t displ(dt*deltaR[iat]);
            if (Lattice.outOfBound(Lattice.toUnit(displ))) return false;
            SingleParticlePos_t newpos(awalker.R[iat]+displ);
            if (!Lattice.isValid(Lattice.toUnit(newpos)))  return false;
            R[iat]=newpos;
          }
      }
    else
      {
        for (int iat=0; iat<deltaR.size(); ++iat)
          R[iat]=awalker.R[iat]+dt*deltaR[iat];
      }

    for (int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if (SK) SK->UpdateAllPart();
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
    if (UseBoundBox)
      {
        for (int iat=0; iat<deltaR.size(); ++iat)
          {
            SingleParticlePos_t displ(dt*deltaR[iat]+drift[iat]);
            if (Lattice.outOfBound(Lattice.toUnit(displ))) return false;
            SingleParticlePos_t newpos(awalker.R[iat]+displ);
            if (!Lattice.isValid(Lattice.toUnit(newpos))) return false;
            R[iat]=newpos;
          }
      }
    else
      {
        for (int iat=0; iat<deltaR.size(); ++iat)
          R[iat]=awalker.R[iat]+dt*deltaR[iat]+drift[iat];
      }

    for (int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if (SK) SK->UpdateAllPart();

    //every move is valid
    return true;
  }

  void
  ParticleSet::makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ)
  {
    myTimers[1]->start();
    if (UseSphereUpdate)
      {
        for (int i=0; i< DistTables.size(); ++i) DistTables[i]->moveOnSphere(*this,displ,iat);
      }
    else
      {
        PosType newpos=activePos+displ;
        for (int i=0; i< DistTables.size(); ++i) DistTables[i]->move(*this,newpos,iat);
      }
    activePtcl=iat;
    activePos=R[iat];
    R[iat]=activePos+displ;

    if (SK && SK->DoUpdate) SK->makeMove(iat,R[iat]);
    //SingleParticlePos_t dum=makeMove(iat,displ);
    myTimers[1]->stop();
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
        for (int i=0; i< DistTables.size(); i++)
          {
            DistTables[i]->update(iat);
          }
        //Do not change SK: 2007-05-18
        if (SK && SK->DoUpdate) SK->acceptMove(iat);
      }
    else
    {
      ostringstream o;
      o << "  Illegal acceptMove " << iat << " != " << activePtcl;
      APP_ABORT(o.str());
    }
  }

  void ParticleSet::makeVirtualMoves(const SingleParticlePos_t& newpos)
  {
    activePtcl=0;
    activePos=R[0];
    for (int i=0; i< DistTables.size(); ++i) DistTables[i]->move(*this,newpos,0);
    R[0]=newpos;
  }

  void ParticleSet::rejectMove(Index_t iat)
  {
    //restore the position by the saved activePos
    R[iat]=activePos;
    //Do not change SK: 2007-05-18
    //if(SK && SK->DoUpdate) SK->rejectMove(iat);
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
    G = awalker.G;
    L = awalker.L;
    //awalker.DataSet.rewind();
    //identical to copyFromBuffer but forbid using buffer
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->copyFromBuffer(buf);
//    if(SK) SK->copyFromBuffer(buf);
//#else
    if (pbyp)
      {
        for (int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
        if (SK) SK->UpdateAllPart();
      }
//#endif
  }

  void ParticleSet::saveWalker(Walker_t& awalker)
  {
    awalker.R=R;
    awalker.G=G;
    awalker.L=L;
    //PAOps<RealType,OHMMS_DIM>::copy(G,awalker.Drift);
    if (SK) SK->UpdateAllPart();
    //awalker.DataSet.rewind();
  }


//  void
//  ParticleSet::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
//    R = awalker.R;
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++)
//    {
//      DistTables[i]->evaluate(*this);
//      DistTables[i]->registerData(buf);
//    }
//    if(SK)
//    {
//      SK->UpdateAllPart();
//      SK->registerData(buf);
//    }
//#else
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
//    if(SK) SK->UpdateAllPart();
//#endif
//  }
//
//  void
//  ParticleSet::registerData(PooledData<RealType>& buf) {
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++)
//    {
//      DistTables[i]->evaluate(*this);
//      DistTables[i]->registerData(buf);
//    }
//    if(SK)
//    {
//      SK->UpdateAllPart();
//      SK->registerData(buf);
//    }
//#else
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
//    if(SK) SK->UpdateAllPart();
//#endif
//  }
//
//  void
//  ParticleSet::updateBuffer(Walker_t& awalker, PooledData<RealType>& buf) {
//    R = awalker.R;
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++)
//    {
//      DistTables[i]->evaluate(*this);
//      DistTables[i]->updateBuffer(buf);
//    }
//    if(SK)
//    {
//      SK->UpdateAllPart();
//      SK->updateBuffer(buf);
//    }
//#else
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
//    if(SK) SK->UpdateAllPart();
//#endif
//  }
//
//  void
//  ParticleSet::updateBuffer(PooledData<RealType>& buf) {
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++)
//    {
//      DistTables[i]->evaluate(*this);
//      DistTables[i]->updateBuffer(buf);
//    }
//    if(SK)
//    {
//      SK->UpdateAllPart();
//      SK->updateBuffer(buf);
//    }
//#else
//    //for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
//    if(SK) SK->UpdateAllPart();
//#endif
//  }
//
//  void
//  ParticleSet::copyToBuffer(PooledData<RealType>& buf) {
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->copyToBuffer(buf);
//    //Do not change SK: 2007-05-18
//    //if(SK) SK->copyToBuffer(buf);
//    if(SK)
//    {//need to calculate the Sk with the current position
//      if(!SK->DoUpdate) SK->UpdateAllPart();
//      SK->copyToBuffer(buf);
//    }
//#else
//    if(SK && !SK->DoUpdate) SK->UpdateAllPart();
//#endif
//  }
//
//  void
//  ParticleSet::copyFromBuffer(PooledData<RealType>& buf)
//  {
//#if defined(PACK_DISTANCETABLES)
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->copyFromBuffer(buf);
//    if(SK) SK->copyFromBuffer(buf);
//#else
//    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
//    if(SK) SK->UpdateAllPart();
//#endif
//  }
//
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
    PropertyList.add("BranchingEnergy");
    PropertyList.add("LocalEnergy");
    PropertyList.add("LocalPotential");

    if (PropertyList.size() != NUMPROPERTIES)
      {
        app_error() << "The number of default properties for walkers  is not consistent." << endl;
        app_error() << "NUMPROPERTIES " << NUMPROPERTIES << " size of PropertyList " << PropertyList.size() << endl;
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
    vector<RealType> newVecHistory=vector<RealType>(leng,0.0);
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

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
