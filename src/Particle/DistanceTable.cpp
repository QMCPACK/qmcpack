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
#include "Particle/DistanceTableData.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/AsymmetricDistanceTableData.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "Lattice/SuperCellTraits.h"

namespace qmcplusplus
{

/**@{instantiation of static data members*/
map<string,DistanceTableData*>  DistanceTable::TableMap;
ParticleSet::ParticleLayout_t* DistanceTable::SimulationCell=0;
/**@}*/

template<int N,unsigned D>
  struct PowerOfN
  {
    enum {value=N*PowerOfN<N,D-1>::value};
  };

template<int N>
  struct PowerOfN<N,0>
  {
    enum {value=1};
  };

const int TwoPowerD=PowerOfN<2,OHMMS_DIM>::value;

/** common definition of a distance 
 *
 * @param lat supercell
 * @param a in/out displacement
 * @return dot(a,a)
 */
template<class T, unsigned D, int SC>
struct DTD_BConds 
{
  inline static T apply(const CrystalLattice<T,D>& lat, TinyVector<T,D>& a) {
    return dot(a,a);
  }
};

/** specialization for a periodic 3D cell
 */
template<class T>
struct DTD_BConds<T,3,SUPERCELL_BULK> {
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
    //T x=fmod(ar[0],1.0); ar[0]=x-static_cast<int>(x*2.0);
    //T y=fmod(ar[1],1.0); ar[1]=y-static_cast<int>(y*2.0);
    //T z=fmod(ar[2],1.0); ar[2]=z-static_cast<int>(z*2.0);
#if defined(HAVE_STD_ROUND)
    ar[0]=ar[0]-round(ar[0]);
    ar[1]=ar[1]-round(ar[1]);
    ar[2]=ar[2]-round(ar[2]);
#else
    T dmy0,dmy1,dmy2;
    T x=modf(ar[0],&dmy0); ar[0]=x-static_cast<int>(x*2.0);
    T y=modf(ar[1],&dmy1); ar[1]=y-static_cast<int>(y*2.0);
    T z=modf(ar[2],&dmy2); ar[2]=z-static_cast<int>(z*2.0);
#endif
    a=lat.toCart(ar);
    return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  }
};

/** specialization for a periodic 3D orthorombic cell
 */
template<class T>
struct DTD_BConds<T,3,SUPERCELL_BULK+TwoPowerD> {
  inline static T apply(const CrystalLattice<T,3>& lat, TinyVector<T,3>& a) 
  {
#if defined(HAVE_STD_ROUND)
    T x=a[0]*lat.OneOverLength[0]; a[0]=lat.Length[0]*(x-round(x));
    T y=a[1]*lat.OneOverLength[1]; a[1]=lat.Length[1]*(y-round(y));
    T z=a[2]*lat.OneOverLength[2]; a[2]=lat.Length[2]*(z-round(z));
#else
    T dmy0,dmy1,dmy2;
    T x=modf(a[0]*lat.OneOverLength[0],&dmy0); a[0]=lat.Length[0]*(x-static_cast<int>(x*2.0));
    T y=modf(a[1]*lat.OneOverLength[1],&dmy1); a[1]=lat.Length[1]*(y-static_cast<int>(y*2.0));
    T z=modf(a[2]*lat.OneOverLength[2],&dmy2); a[2]=lat.Length[2]*(z-static_cast<int>(z*2.0));
#endif
    return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  }
};

/** specialization for a periodic 2D cell
 */
template<class T>
struct DTD_BConds<T,2,SUPERCELL_BULK> {
  inline static T apply(const CrystalLattice<T,2>& lat, TinyVector<T,2>& a) {
    TinyVector<T,2> ar(lat.toUnit(a));
#if defined(HAVE_STD_ROUND)
    ar[0]=ar[0]-round(ar[0]);
    ar[1]=ar[1]-round(ar[1]);
#else
    T dmy0,dmy1;
    T x=modf(ar[0],&dmy0); ar[0]=x-static_cast<int>(x*2.0);
    T y=modf(ar[1],&dmy1); ar[1]=y-static_cast<int>(y*2.0);
#endif
    a=lat.toCart(ar);
    return a[0]*a[0]+a[1]*a[1];
  }
};

/** specialization for a periodic 2D orthorombic cell
 */
template<class T>
struct DTD_BConds<T,2,SUPERCELL_BULK+TwoPowerD> {
  inline static T apply(const CrystalLattice<T,2>& lat, TinyVector<T,2>& a) 
  {
#if defined(HAVE_STD_ROUND)
    T x=a[0]*lat.OneOverLength[0]; a[0]=lat.Length[0]*(x-round(x));
    T y=a[1]*lat.OneOverLength[1]; a[1]=lat.Length[1]*(y-round(y));
#else
    T dmy0,dmy1;
    T x=modf(a[0]*lat.OneOverLength[0],&dmy0); a[0]=lat.Length[0]*(x-static_cast<int>(x*2.0));
    T y=modf(a[1]*lat.OneOverLength[1],&dmy1); a[1]=lat.Length[1]*(y-static_cast<int>(y*2.0));
#endif
    return a[0]*a[0]+a[1]*a[1];
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
    //if(SuperCellType<DIM>::apply(s.Lattice.BoxBConds) == SUPERCELL_OPEN) 
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,s);
    else
    {
      if(s.Lattice.DiagonalOnly)
        dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,s);
      else
        dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,s);
    }

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
    //if(SuperCellType<DIM>::apply(s.Lattice.BoxBConds) == SUPERCELL_OPEN) 
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,t);
    else 
    {
      if(s.Lattice.DiagonalOnly)
        dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,t);
      else
        dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,t);
    }

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
} //namespace qmcplusplus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
