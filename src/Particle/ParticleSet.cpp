//////////////////////////////////////////////////////////////////
//W.create
// (c) Copyright 2003  by Jeongnim Kim
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

namespace qmcplusplus {

  ///object counter 
  int  ParticleSet::PtclObjectCounter = 0;

  ParticleSet::ParticleSet(): SK(0), ParentTag(-1){ 
    initParticleSet();

    initPropertyList();
  }

  ParticleSet::ParticleSet(const ParticleSet& p): SK(0), mySpecies(p.getSpeciesSet()), 
  ParentTag(p.tag())
  {
    initBase();
    initParticleSet();
    assign(p);

    ostringstream o;
    o<<p.getName()<<ObjectTag;
    this->setName(o.str());
    app_log() << "  Copying a particle set " << p.getName() << " to " << this->getName() << endl;
    PropertyList.Names=p.PropertyList.Names;
    PropertyList.Values=p.PropertyList.Values;

    PropertyHistory=  p.PropertyHistory;
    phLength= p.phLength;
    
    //construct the distance tables with the same order
    //first is always for this-this paier
    for(int i=1;i<p.DistTables.size(); ++i) addTable(p.DistTables[i]->origin());

    if(p.SK) 
    {
      R.InUnit=p.R.InUnit;
      createSK();
      SK->DoUpdate=p.SK->DoUpdate;
    }

    if(p.Sphere.size())
    {
      resizeSphere(p.Sphere.size());
    }
  }


  ParticleSet::~ParticleSet() 
  {
    DEBUG_MEMORY("ParticleSet::~ParticleSet");
    delete_iter(DistTables.begin(), DistTables.end());
    if(SK) delete SK;
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
  }

  ///write to a ostream
  bool ParticleSet::get(ostream& os) const {

    os << "  ParticleSet " << getName() << " : ";
    for(int i=0; i<SubPtcl.size(); i++) os << SubPtcl[i] << " ";
    os <<"\n\n    " << LocalNum << "\n\n";

    for(int i=0; i<LocalNum; i++) {
      os << "    " << mySpecies.speciesName[GroupID[i]]  << R[i] << endl;
    }
    return true;
  }
    
  ///read from istream
  bool ParticleSet::put(istream& is) { 
    return true;
  }

  ///reset member data
  void ParticleSet::reset() { 
  }

  ///read the particleset
  bool ParticleSet::put(xmlNodePtr cur){
    return true;
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
    if(DistTables.empty())
    {
      DistTables.reserve(4);
      DistTables.push_back(createDistanceTable(*this));
      //add  this-this pair
      myDistTableMap.clear();
      myDistTableMap[ObjectTag]=0;
      app_log() << "  ... ParticleSet::addTable Create Table #0 " << DistTables[0]->Name << endl;
      
      if(psrc.tag() == ObjectTag) return 0;
    } 
    
    if(psrc.tag() == ObjectTag) {
      app_log() << "  ... ParticleSet::addTable Reuse Table #" << 0 << " " << DistTables[0]->Name <<endl;
      return 0;
    }

    int tsize=DistTables.size(),tid;
    map<int,int>::iterator tit(myDistTableMap.find(psrc.tag()));
    if(tit == myDistTableMap.end())
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

  void ParticleSet::update(int iflag) {

    //apply Boundary condition
    //R.setUnit(0);
    //double xL=Lattice.R(0,0);
    //double yL=Lattice.R(1,1);
    //double zL=Lattice.R(2,2);
    //for(int iat=0; iat<LocalNum; iat++) {
    //  if(R[iat][0]<0) R[iat][0]+=xL;
    //  else if(R[iat][0]>xL) R[iat][0]-=xL;

    //  if(R[iat][1]<0) R[iat][1]+=yL;
    //  else if(R[iat][0]>yL) R[iat][1]-=yL;

    //  if(R[iat][2]<0) R[iat][2]+=zL;
    //  else if(R[iat][0]>zL) R[iat][2]-=zL;
    //}

    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if(SK) SK->UpdateAllPart();
  }  

  void ParticleSet::update(const ParticlePos_t& pos) { 
    R = pos;
    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);
    if(SK && !SK->DoUpdate) SK->UpdateAllPart();
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
  ParticleSet::makeMove(Index_t iat, const SingleParticlePos_t& displ) {
    activePtcl=iat;
    activePos=R[iat]; //save the current position
    SingleParticlePos_t newpos(activePos+displ);
    for(int i=0; i< DistTables.size(); ++i) DistTables[i]->move(*this,newpos,iat);
    
    R[iat]=newpos;
    //Do not change SK: 2007-05-18
    //Change SK only if DoUpdate is true: 2008-09-12
    if(SK && SK->DoUpdate) SK->makeMove(iat,newpos);
    return newpos;
  }

  void
  ParticleSet::makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ) {
    SingleParticlePos_t dum=makeMove(iat,displ);
    //activePtcl=iat;
    //activePos=R[iat]; //save the current position
    //R[iat]=activePos+displ;
    //for(int i=0; i< DistTables.size(); ++i) 
    //  DistTables[i]->moveOnSphere(*this,displ,iat);
  }
  
  /** update the particle attribute by the proposed move
   *@param iat the particle index
   *
   *When the activePtcl is equal to iat, overwrite the position and update the
   *content of the distance tables.
   */
  void ParticleSet::acceptMove(Index_t iat) {
    if(iat == activePtcl) {
      //Update position + distance-table
      for(int i=0; i< DistTables.size(); i++) {
        DistTables[i]->update(iat);
      }
      //Do not change SK: 2007-05-18
      if(SK && SK->DoUpdate) SK->acceptMove(iat);
    }
    else
    {
      APP_ABORT("  Illegal acceptMove ");
    }
  }

  void ParticleSet::rejectMove(Index_t iat) {
    //restore the position by the saved activePos
    R[iat]=activePos;
    //Do not change SK: 2007-05-18
    //if(SK && SK->DoUpdate) SK->rejectMove(iat);
  }

  /** resize Sphere by the LocalNum
   * @param nc number of centers to which Spherical grid will be assigned.
   */
  void ParticleSet::resizeSphere(int nc) {
    int nsadd=nc-Sphere.size();
    while(nsadd>0) {
      Sphere.push_back(new ParticlePos_t);
      --nsadd;
    }
  }

  void 
  ParticleSet::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
    R = awalker.R;
    for(int i=0; i< DistTables.size(); i++) 
    {
      DistTables[i]->evaluate(*this);
      DistTables[i]->registerData(buf);
    }
    if(SK)
    {
      SK->UpdateAllPart();
      SK->registerData(buf);
    }
  }

  void 
  ParticleSet::registerData(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) 
    {
      DistTables[i]->evaluate(*this);
      DistTables[i]->registerData(buf);
    }
    if(SK)
    {
      SK->UpdateAllPart();
      SK->registerData(buf);
    }
  }
  
  void 
  ParticleSet::updateBuffer(Walker_t& awalker, PooledData<RealType>& buf) {
    R = awalker.R;
    for(int i=0; i< DistTables.size(); i++) 
    {
      DistTables[i]->evaluate(*this);
      DistTables[i]->updateBuffer(buf);
    }
    if(SK)
    {
      SK->UpdateAllPart();
      SK->updateBuffer(buf);
    }
  }
    
  void 
  ParticleSet::updateBuffer(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) 
    {
      DistTables[i]->evaluate(*this);
      DistTables[i]->updateBuffer(buf);
    }
    if(SK)
    {
      SK->UpdateAllPart();
      SK->updateBuffer(buf);
    }
  }
    
  void 
  ParticleSet::copyToBuffer(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) DistTables[i]->copyToBuffer(buf);
    //Do not change SK: 2007-05-18
    //if(SK) SK->copyToBuffer(buf);
    if(SK)
    {//need to calculate the Sk with the current position
      if(!SK->DoUpdate) SK->UpdateAllPart();
      SK->copyToBuffer(buf);
    }
  }
  
  void 
  ParticleSet::copyFromBuffer(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) DistTables[i]->copyFromBuffer(buf);
    if(SK) SK->copyFromBuffer(buf);
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
    PropertyList.add("LocalEnergy");
    PropertyList.add("LocalPotential");

    if(PropertyList.size() != NUMPROPERTIES)
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
  
   void ParticleSet::setPropertyHistoryLength(int leng)
  {
    if (leng>0) phLength=leng;
  }
    
  int ParticleSet::addPropertyHistory(int leng)
  {
    int newL = PropertyHistory.size();
    vector<double> newVecHistory(leng,0.0);
    PropertyHistory.push_back(newVecHistory);
    return newL;
  }

  void ParticleSet::addPropertyHistoryPoint(int Rindex, double data)
  {
    PropertyHistory[Rindex].insert(PropertyHistory[Rindex].begin(),1,data);
    PropertyHistory[Rindex].pop_back();
  }
  
  double ParticleSet::getPropertyHistoryPoint(int Rindex, double Tindex)
  {
    return PropertyHistory[Rindex][Tindex];
  }
    
   double ParticleSet::getPropertyHistoryAvg(int index)
  {
    double mean=0.0;
    vector<double>::iterator phStart=PropertyHistory[index].begin();
    for(;phStart!=PropertyHistory[index].end();phStart++){
      mean+= (*phStart);
    }
    return (mean/PropertyHistory[index].size());
  }
    
  double ParticleSet::getPropertyHistorySum(int index, int endN)
  {
    double mean=0.0;
    vector<double>::iterator phStart=PropertyHistory[index].begin();
    for(int i=0;((phStart!=PropertyHistory[index].end())&(i<endN));phStart++,i++){
      mean+= (*phStart);
    }
    return mean ;
  }
  
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
