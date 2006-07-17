//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
//#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/AsymmetricDistanceTableData.h"
#include "ParticleIO/ParticleLayoutIO.h"
using namespace qmcplusplus;

/**@{instantiation of static data members*/
map<string,DistanceTableData*>  DistanceTable::TableMap;
ParticleSet::ParticleLayout_t* DistanceTable::SimulationCell=0;
/**@}*/

/** dummy class to detemine the supercell type **/
template<unsigned D>
struct SuperCellType {};

/** specialization of SuperCellType for 3-dimensional cell
 */
template<>
struct SuperCellType<3> {
  /** convert box to an integer
   * @param box 3-dimensional boolean vector
   *
   * When all the directions are open, returns 0.
   * - 7 (1,1,1) bulk system
   * - 3 (1,1,0) slab system
   * - 1 (1,0,0) wire system
   */
  inline static int apply(const TinyVector<int,3>& box) {
    return box[0]+2*(box[1]+box[2]*2);
  }
};

template<class T, unsigned D, int SC>
struct DTD_BConds {};

template<class T>
struct DTD_BConds<T,3,0> {
  inline static T apply(const CrystalLattice<T,3>& lat, TinyVector<T,3>& a) {
    return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  }
};

template<class T>
struct DTD_BConds<T,3,7> {
  inline static T apply(const CrystalLattice<T,3>& lat, TinyVector<T,3>& a) {
    TinyVector<T,3> ar(lat.toUnit(a));
    /*
    if(ar[0]<-0.5) ar[0]+=1.0; 
    else if(ar[0]>=0.5) ar[0]-=1.0;
    if(ar[1]<-0.5) ar[1]+=1.0; 
    else if(ar[1]>=0.5) ar[1]-=1.0;
    if(ar[2]<-0.5) ar[2]+=1.0; 
    else if(ar[2]>=0.5) ar[2]-=1.0;
    */
    T x=fmod(ar[0],1.0);
    T y=fmod(ar[1],1.0);
    T z=fmod(ar[2],1.0);
    ar[0]=x-static_cast<int>(x*2.0);
    ar[1]=y-static_cast<int>(y*2.0);
    ar[2]=z-static_cast<int>(z*2.0);
    a=lat.toCart(ar);
    return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  }
};

/**@ingroup nnlist
 * @brief class to apply No Boundary conditions for distance evaluation
 *
 *NoBConds stands for No Boundary Conditions and is intended for
 finite systems with open or vanishing boundary conditions.
 Use a simple dot product assuming cartesian coordinates
 */
template<class T, unsigned D>
struct NoBConds {
  inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
    return dot(a,a);
  }
};

/**@ingroup nnlist
 * @brief class to apply Periodic Boundary conditions for distance evaluation
 *
 *PeriodicBConds stands for Periodic Boundary Conditions and is intended 
 for periodic systems. The Cartesian distances are evaluated
 according to the minimum-image convention.
 */
template<class T, unsigned D>
struct PeriodicBConds {
  inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
    TinyVector<T,D> ar(lat.toUnit(a));
    for(int idim=0; idim<D; idim++) {
      if(lat.BoxBConds[idim]) {
        if(ar[idim]<-0.5) ar[idim]+=1.0; 
        else if(ar[idim]>=0.5) ar[idim]-=1.0;
      }
    }
    a=lat.toCart(ar);
    return dot(a,a);
  }
};

void DistanceTable::createSimulationCell(xmlNodePtr cur) {
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
 *\param aname of a new DistanceTableData
 *\return index of the distance table with the name
 */
DistanceTableData*
DistanceTable::add(ParticleSet& s, const char* aname) {

  string newname;

  if(aname) {
    newname = aname;
  } else {
    newname = s.getName(); newname.append(s.getName());
  }

  map<string,DistanceTableData*>::iterator it = TableMap.find(newname);

  //the named pair does not exist, add a new symmetric metrics
  if(it == TableMap.end()) {
    //LOGMSG("Distance table " << newname << " is created.")
    DistanceTableData* dt=0;
    if(SuperCellType<OHMMS_DIM>::apply(s.Lattice.BoxBConds) == SUPERCELL_OPEN) 
      dt = new SymmetricDTD<DTD_BConds<OHMMS_PRECISION,OHMMS_DIM,0> >(s,s);
    else
      dt = new SymmetricDTD<DTD_BConds<OHMMS_PRECISION,OHMMS_DIM,7> >(s,s);

    //set the name of the table
    dt->setName(newname);
    TableMap[newname] = dt;
    s.addTable(dt);
    //add to the list
    return dt;
  } else {
    //LOGMSG("Distance table " << newname << " is reused")
    return (*it).second;
  }
}

/** Adding AsymmetricDTD to the list, e.g., el-nuclei distance table
 *\param s source particle set
 *\param t target particle set
 *\param aname of a new DistanceTableData
 *\return index of the distance table with the name
 */
DistanceTableData*
DistanceTable::add(const ParticleSet& s, ParticleSet& t, const char* aname) {

  string newname;
  if(aname) {
    newname = aname;
  } else {
    newname = s.getName();
    newname.append(t.getName());
  }

  map<string,DistanceTableData*>::iterator it = TableMap.find(newname);

  ///the named pair does not exist, add a new asymmetric metrics
  if(it == TableMap.end()) {
    DistanceTableData* dt=0;
    if(SuperCellType<OHMMS_DIM>::apply(s.Lattice.BoxBConds) == SUPERCELL_OPEN) 
      dt = new AsymmetricDTD<DTD_BConds<OHMMS_PRECISION,OHMMS_DIM,0> >(s,t);
    else 
      dt = new AsymmetricDTD<DTD_BConds<OHMMS_PRECISION,OHMMS_DIM,7> >(s,t);

    //set the name of the table
    dt->setName(newname);

    t.addTable(dt);
    TableMap[newname] = dt;
    return dt;
  } else {
    //LOGMSG("Distance table " << newname << " is reused")
    return (*it).second;
  }
}

void 
DistanceTable::removeTable(const string& tname) {
  map<string,DistanceTableData*>::iterator it = TableMap.find(tname);
  if(it != TableMap.end()) {
    delete (*it).second;
    TableMap.erase(it);
  }
}

void DistanceTable::create(int walkers) {
  map<string,DistanceTableData*>::iterator it = TableMap.begin();
  map<string,DistanceTableData*>::iterator it_end = TableMap.end();
  while(it != it_end) {
    (*it).second->create(walkers);
    ++it;
  }
}

void DistanceTable::reset() {
}

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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
