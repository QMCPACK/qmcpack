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

  int  ParticleSet::PtclObjectCounter = 0;

  ParticleSet::ParticleSet(): SK(0), ParentTag(-1){ 
    initParticleSet();

    initPropertyList();
  }

  ParticleSet::ParticleSet(const ParticleSet& p): SK(0), mySpecies(p.getSpeciesSet()), ParentTag(p.tag()){
    initBase();
    initParticleSet();
    assign(p);

    if(p.SK) 
    {
      createSK();
    }

    PropertyList.Names=p.PropertyList.Names;
    PropertyList.Values=p.PropertyList.Values;
  }


  ParticleSet::~ParticleSet() {
    if(SK) delete SK;
    delete_iter(Sphere.begin(), Sphere.end());
  }

  void ParticleSet::initParticleSet() {
    ObjectTag = PtclObjectCounter;
    PtclObjectCounter++;
    G.setTypeName(ParticleTags::postype_tag);
    G.setObjName("grad");
    L.setTypeName(ParticleTags::scalartype_tag);
    L.setObjName("laplacian");
    addAttribute(G);
    addAttribute(L);
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

  /** add a distance table to DistTables list
   * @param d_table pointer to a DistanceTableData to be added
   *
   * DistTables is a list of DistanceTables which are updated by MC moves.
   */
  void ParticleSet::addTable(DistanceTableData* d_table) {
    int oid=d_table->origin().tag();
    int i=0;
    int dsize=DistTables.size();
    while(i<dsize) {
      if(oid == DistTables[i]->origin().tag()) //table already exists
        return;
      ++i;
    }
    DistTables.push_back(d_table);
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
    if(SK) SK->UpdateAllPart();
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
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->move(*this,newpos,iat);
    } 
    
    R[iat]=newpos;
    //Do not change SK: 2007-05-18
    //if(SK) SK->makeMove(iat,newpos);
    return newpos;
  }

  void
  ParticleSet::makeMoveOnSphere(Index_t iat, const SingleParticlePos_t& displ) {
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->moveOnSphere(*this,displ,iat);
    }
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
      //if(SK) SK->acceptMove(iat);
    }
  }

  void ParticleSet::rejectMove(Index_t iat) {
    //restore the position by the saved activePos
    R[iat]=activePos;
    //Do not change SK: 2007-05-18
    //if(SK) SK->rejectMove(iat);
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
      SK->UpdateAllPart();
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

  void ParticleSet::clearDistanceTables() {
    //Physically remove the tables
    for(int i=0; i< DistTables.size(); i++) DistanceTable::removeTable(DistTables[i]->getName());
    DistTables.erase(DistTables.begin(),DistTables.end());
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
