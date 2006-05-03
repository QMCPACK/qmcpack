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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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

namespace qmcplusplus {

  int  ParticleSet::PtclObjectCounter = 0;

  ParticleSet::ParticleSet(): SK(NULL), SKOld(NULL) {
    initParticleSet();

    initPropertyList();
  }

  ParticleSet::ParticleSet(const ParticleSet& p): SK(NULL), SKOld(NULL) {
    initBase();
    initParticleSet();
    assign(p);

    PropertyList.Name=p.PropertyList.Name;
    PropertyList.Values=p.PropertyList.Values;
  }


  ParticleSet::~ParticleSet() {
    if(SK) delete SK;
    if(SKOld) delete SKOld;
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

  void ParticleSet::addTable(DistanceTableData* d_table) {
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

    //Update SK if used. 
    //Swap pointers for old and current so that old is accurate.
    //This is for all-particle moves, so the contents of old before this point are irrelevant.
    if(SK){	
      void* tmp = SK;
      SK = SKOld;
      SKOld = static_cast<StructFact*>(tmp);
      //Evaluate current structurefactor.
      SK->UpdateAllPart();
    }
  }  

  void ParticleSet::update(const ParticlePos_t& pos) { 
    R = pos;
    for(int i=0; i< DistTables.size(); i++) DistTables[i]->evaluate(*this);

    //Update SK if used. 
    //Swap pointers for old and current so that old is accurate.
    //This is for all-particle moves, so the contents of old before this point are irrelevant.
    if(SK){	
      void* tmp = SK;
      SK = SKOld;
      SKOld = static_cast<StructFact*>(tmp);
      //Evaluate current structurefactor.
      SK->UpdateAllPart();
    }
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

    //Update SK if used. Use smart-moving. 
    //In this case we can't simply swap pointers SK and SKOld because we evaluate the new
    //SK by adding the change due to particle move. For this reason, the data in SKOld
    //must be copied from SK before the update of SK

    //R is not changed by ->move. Use newpos and activepos as new and old.
    if(SK) {
      *SKOld = *SK; //Copy data from current using overloaded assignment operator.
      SK->Update1Part(activePos,newpos,iat,GroupID[iat]); //old = activePos, new=newpos
    }

    return R[iat]=newpos;
    //activePos=R[iat]+displ;
    //for(int i=0; i< DistTables.size(); i++) {
    //  DistTables[i]->move(*this,activePos,iat);
    //}
    //return activePos;
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
    }
  }

  void ParticleSet::rejectMove(Index_t iat) {
    //restore the position by the saved activePos
    R[iat]=activePos;
    //If the move is rejected then we must restore the SKOld data back to the current.
    //Swap the pointers for a fast way.
    //Accept the proposed structurefactor. Just swap the pointers.
    if(SK){
      void* tmp = SK;
      SK = SKOld;
      SKOld = static_cast<StructFact*>(tmp);
    }
  }

  /** resize Sphere by the LocalNum
   * @param nc number of centers to which Spherical grid will be assigned.
   */
  void ParticleSet::resizeSphere(int nc) {
    if(Sphere.size() != nc) {
      for(int i=0; i<nc; i++) Sphere.push_back(new ParticlePos_t);
    }
  }

  void 
  ParticleSet::registerData(Walker_t& awalker, PooledData<RealType>& buf) {
    R = awalker.R;
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->evaluate(*this);
      DistTables[i]->registerData(buf);
    }
    if(SK){
      SK->UpdateAllPart();
      SK->registerData(buf);
    }
  }

  void 
  ParticleSet::registerData(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->evaluate(*this);
      DistTables[i]->registerData(buf);
    }
    if(SK){
      SK->UpdateAllPart();
      SK->registerData(buf);
    }
  }
  
  void 
  ParticleSet::updateBuffer(Walker_t& awalker, PooledData<RealType>& buf) {
    R = awalker.R;
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->evaluate(*this);
      DistTables[i]->updateBuffer(buf);
    }
    if(SK){
      SK->UpdateAllPart();
      SK->updateBuffer(buf);
    }
  }
    
  void 
  ParticleSet::copyToBuffer(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->copyToBuffer(buf);
    }
    if(SK)SK->copyToBuffer(buf);
  }
  
  void 
  ParticleSet::copyFromBuffer(PooledData<RealType>& buf) {
    for(int i=0; i< DistTables.size(); i++) {
      DistTables[i]->copyFromBuffer(buf);
    }
    if(SK)SK->copyFromBuffer(buf);
  }

  void ParticleSet::initPropertyList() {
    //Need to add the default Properties according to the enumeration
    PropertyList.add("LogPsi");
    PropertyList.add("SignPsi");
    PropertyList.add("UmbrellaWeight");
    PropertyList.add("LocalEnergy");
    PropertyList.add("LocalPotential");
  }

  void ParticleSet::clearDistanceTables() {
    //Physically remove the tables
    for(int i=0; i< DistTables.size(); i++) {
      DistanceTable::removeTable(DistTables[i]->getName());
    }
    DistTables.erase(DistTables.begin(),DistTables.end());
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
