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
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/AsymmetricDistanceTableData.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Lattice/ParticleBConds.h"

namespace qmcplusplus
{

  /**@{instantiation of static data members*/
  //map<string,DistanceTableData*>  DistanceTable::TableMap;
  ParticleSet::ParticleLayout_t* DistanceTable::SimulationCell=0;
  /**@}*/

  /** Adding SymmetricDTD to the list, e.g., el-el distance table
   *\param s source/target particle set
   *\return index of the distance table with the name
   */
  DistanceTableData* createDistanceTable(ParticleSet& s) 
  {
    typedef OHMMS_PRECISION RealType;
    enum {DIM=OHMMS_DIM};
    DistanceTableData* dt=0;
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,s);
    else
    {
      if(s.Lattice.DiagonalOnly)
      {
        app_log() << "Distance table specialized for an Orthorhombic cell " << endl;
        dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,s);
      }
      else
      {
        app_log() << "Distance table specialized for a generic cell " << endl;
        dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,s);
      }
    }
    ostringstream o;
    o << s.getName() << "_" << s.getName();
    dt->Name=o.str();//assign the table name
    return dt;
  }

  /** Adding SymmetricDTD to the list, e.g., el-el distance table
   *\param s source/target particle set
   *\return index of the distance table with the name
   */
  DistanceTableData* createDistanceTable(const ParticleSet& s, ParticleSet& t) 
  {
    typedef OHMMS_PRECISION RealType;
    enum {DIM=OHMMS_DIM};
    DistanceTableData* dt=0;
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,t);
    else 
    {
      if(s.Lattice.DiagonalOnly)
      {
        app_log() << "Distance table specialized for an Orthorhombic cell " << endl;
        dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,t);
      }
      else
      {
        app_log() << "Distance table specialized for a generic cell " << endl;
        dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,t);
      }
    }
    ostringstream o;
    o << s.getName() << "_" << t.getName();
    dt->Name=o.str();//assign the table name
    return dt;
  }


  void DistanceTable::createSimulationCell(xmlNodePtr cur) 
  {
    if(cur != NULL) {
      if(SimulationCell == 0) {
        SimulationCell = new ParticleLayout_t;
      }
      LatticeParser a(*SimulationCell);
      a.put(cur);
    }
  }

  /** Adding SymmetricDTD to the list, e.g., el-el distance table
   *\param s source/target particle set
   *\return DistanceTableData*
   */
  DistanceTableData* DistanceTable::add(ParticleSet& s)
  {
    int tid=s.addTable(s);
    return s.DistTables[tid];
  }

  /** Adding AsymmetricDTD to the list, e.g., el-nuclei distance table
   *\param s source particle set
   *\param t target particle set
   *\return DistanceTableData*
   */
  DistanceTableData* DistanceTable::add(const ParticleSet& s, ParticleSet& t) 
  {
    int tid=t.addTable(s);
    return t.DistTables[tid];
  }


//void 
//DistanceTable::removeTable(const string& tname) {
//  map<string,DistanceTableData*>::iterator it = TableMap.find(tname);
//  if(it != TableMap.end()) {
//    delete (*it).second;
//    TableMap.erase(it);
//  }
//}
//
//void DistanceTable::create(int walkers) {
//  map<string,DistanceTableData*>::iterator it = TableMap.begin();
//  map<string,DistanceTableData*>::iterator it_end = TableMap.end();
//  while(it != it_end) {
//    (*it).second->create(walkers);
//    ++it;
//  }
//}
//
//void DistanceTable::reset() {
//}

//May need to make it singleton
// class DistanceTableSingleton {
// public:
//   typedef DistanceTable<TinyVector<double,3> > Table_t;
//   static void set(int np, int nc=1) {ref_.set(np,nc);}
//   static Table_t& getTable() {return ref_;}
//   ///static const Table_t& getTable() const { return ref_; }
// private:
//   static Table_t ref_;
// };
// ///instantiate the singleton
// DistanceTable<TinyVector<double,3> > DistanceTableSingleton::ref_;
} //namespace qmcplusplus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
