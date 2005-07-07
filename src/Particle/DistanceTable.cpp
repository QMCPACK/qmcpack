//////////////////////////////////////////////////////////////////
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
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/AsymmetricDistanceTableData.h"
using namespace ohmmsqmc;

//uncomment to test PBC
//#define TEST_PBC_DIST

//@todo introduce enum to automatically generate the distance table data with BC
//enum { No=0, PNN, PPN, PPP};

/**class for boundary conditions for distance evaluation
 *@brief use a simple dot product assuming cartesian coordinates
 *@todo implements periodic boundary conditions and factory function
 */
template<class T, unsigned D>
struct NoBConds {
  inline static T apply(const CrystalLattice<T,D>& lat, 
			const TinyVector<T,D>& a) {
    return dot(a,a);
  }
};

template<class T, unsigned D>
struct PeriodicBConds {
  inline static T apply(const CrystalLattice<T,D>& lat, 
			const TinyVector<T,D>& a) {
    TinyVector<T,D> ar(lat.toUnit(a));
    for(int idim=0; idim<D; idim++) 
      if(ar[idim]<-0.5) ar[idim]+=1.0; 
      else if(ar[idim]>=0.5) ar[idim]-=1.0;
    return lat.Dot(ar,ar);
  }
};

template<class T>
struct PeriodicBConds<T,3>{
  inline static T apply(const CrystalLattice<T,3>& lat, 
			const TinyVector<T,3>& a) {
    TinyVector<T,3> ar(lat.toUnit(a));
    if(ar[0]<-0.5) ar[0]+=1.0; else if(ar[0]>0.5) ar[0]-=1.0;
    if(ar[1]<-0.5) ar[1]+=1.0; else if(ar[1]>0.5) ar[1]-=1.0;
    if(ar[2]<-0.5) ar[2]+=1.0; else if(ar[2]>0.5) ar[2]-=1.0;
    return ar[0]*ar[0]*lat.M(0)+ar[1]*ar[1]*lat.M(4)+ar[2]*ar[2]*lat.M(8);
    //return lat.Dot(ar,ar);
  }
};


/**@{instantiation of static data members*/
vector<bool>                DistanceTable::Updated;
vector<DistanceTableData*>  DistanceTable::TableList;
map<string,int>             DistanceTable::TableMap;
vector<int>                 DistanceTable::VisitorID;
/**@}*/

/*!\fn int DistanceTable::add(const ParticleSet& s, const char* aname) 
 *\param s source/target particle set
 *\param aname of a new DistanceTableData
 *\return index of the distance table with the name
 *\brief Adding SymmetricDTD to the list, e.g., el-el distance table
 */
int
DistanceTable::add(const ParticleSet& s, const char* aname) {

  string newname;

  if(aname) {
    newname = aname;
  } else {
    newname = s.getName(); newname.append(s.getName());
  }

  map<string,int>::iterator it = TableMap.find(newname);

  ///the named pair does not exist, add a new symmetric metrics
  if(it == TableMap.end()) {
    LOGMSG("Distance table " << newname << " is created.")
    int n = TableList.size();
#if defined(TEST_PBC_DIST)
    TableList.push_back(new SymmetricDTD<PeriodicBConds<double,3> >(s,s));
#else
    TableList.push_back(new SymmetricDTD<NoBConds<double,3> >(s,s));
#endif
    TableMap[newname] = n;
    VisitorID.push_back(s.tag());
    return n;
  } else {
    LOGMSG("Distance table " << newname << " is reused")
    return (*it).second;
  }
}

/*!\fn int DistanceTable::add(const ParticleSet& s,const ParticleSet& t, const char* aname) 
 *\param s source particle set
 *\param s target particle set
 *\param aname of a new DistanceTableData
 *\return index of the distance table with the name
 *\brief Adding AsymmetricDTD to the list, e.g., el-nuclei distance table
 */
int
DistanceTable::add(const ParticleSet& s, const ParticleSet& t,  
		   const char* aname) {

  string newname;
  if(aname) {
    newname = aname;
  } else {
    newname = s.getName();
    newname.append(t.getName());
  }

  map<string,int>::iterator it = TableMap.find(newname);

  ///the named pair does not exist, add a new asymmetric metrics
  if(it == TableMap.end()) {
    LOGMSG("Distance table " << newname << " is created.")
    int n = TableList.size();
#if defined(TEST_PBC_DIST)
    TableList.push_back(new AsymmetricDTD<PeriodicBConds<double,3> > (s,t));
#else
    TableList.push_back(new AsymmetricDTD<NoBConds<double,3> > (s,t));
#endif
    TableMap[newname] = n;
    VisitorID.push_back(t.tag());
    return n;
  } else {
    LOGMSG("Distance table " << newname << " is reused")
    return (*it).second;
  }
}

void 
DistanceTable::getTables(int ptag, vector<DistanceTableData*>&  tables) {
  ///add the table if ptag matches to the source or visitor
  for(int i=0; i<TableList.size(); i++) 
    if(ptag == VisitorID[i] || ptag ==  TableList[i]->origin().tag()) 
      tables.push_back(TableList[i]);
}

void DistanceTable::create(int walkers) {
  for(int i=0; i<TableList.size(); i++) TableList[i]->create(walkers);
}

void DistanceTable::reset() {
  for(int i=0; i<Updated.size(); i++) Updated[i] = false;
}

void DistanceTable::update(ParticleSet& t) {
  for(int i=0; i< TableList.size(); i++) 
    if(t.tag() == VisitorID[i]) TableList[i]->evaluate(t);
}

void DistanceTable::update(WalkerSetRef& w) {
  int id = w.tag();
  for(int i=0; i< TableList.size(); i++) 
    if(id == VisitorID[i]) TableList[i]->evaluate(w);
}

void DistanceTable::registerData(PooledData<RealType>& buf) {
  for(int i=0; i<TableList.size(); i++) TableList[i]->registerData(buf);
}

void DistanceTable::copyToBuffer(PooledData<RealType>& buf) {
  for(int i=0; i<TableList.size(); i++) TableList[i]->copyToBuffer(buf);
}

void DistanceTable::copyFromBuffer(PooledData<RealType>& buf) {
  for(int i=0; i<TableList.size(); i++) TableList[i]->copyFromBuffer(buf);
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
