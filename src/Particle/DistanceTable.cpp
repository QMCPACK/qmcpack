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
#include "Lattice/ParticleBConds.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Particle/AsymmetricDistanceTableData.h"
namespace qmcplusplus
{

  /** Adding SymmetricDTD to the list, e.g., el-el distance table
   *\param s source/target particle set
   *\return index of the distance table with the name
   */
  DistanceTableData* createDistanceTable(ParticleSet& s) 
  {
    typedef OHMMS_PRECISION RealType;
    enum {DIM=OHMMS_DIM};

    int sc=s.Lattice.SuperCellEnum;
    DistanceTableData* dt=0;

    ostringstream o;

    o << "  Distance table for AA: source/target = " << s.getName() << "\n";

    if(sc == SUPERCELL_BULK)
    {
      if(s.Lattice.DiagonalOnly)
      {
        o << "    PBC=bulk Orthorhombic=yes\n";
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_BULK+TwoPowerD>(s,s);
      }
      else
      {
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_BULK>(s,s);
        o << "    PBC=bulk Orthorhombic=no\n" 
          << "    Setting Rmax = " << s.Lattice.SimulationCellRadius;
      }
    }
    else if(sc == SUPERCELL_SLAB)
    {
      if(s.Lattice.DiagonalOnly)
      {
        o << "    PBC=slab Orthorhombic=yes\n";
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_SLAB+TwoPowerD>(s,s);
      }
      else
      {
        o << "    PBC=slab Orthorhombic=no\n" 
          << "    Setting Rmax = " << s.Lattice.SimulationCellRadius;
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_SLAB>(s,s);
      }
    }
    else if(sc == SUPERCELL_WIRE)
    {
      o << "    PBC=wire Orthorhombic=NA\n";
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_WIRE>(s,s);
    }
    else  //open boundary condition
    {
      o << "    PBC=open Orthorhombic=NA\n";
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,s);
    }

    dt->UseBoundBox= (sc!=SUPERCELL_OPEN);

    ostringstream p;
    p << s.getName() << "_" << s.getName();
    dt->Name=p.str();//assign the table name

    if(dt->UseBoundBox) 
      o << " Using bonding box/reduced coordinates ";
    else
      o << " using Cartesian coordinates ";

    app_log() << o.str() << endl;

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
    //int sc=s.Lattice.SuperCellEnum;
    int sc=t.Lattice.SuperCellEnum;

    ostringstream o;
    o << "  Distance table for AB: source = " << s.getName() << " target = " << t.getName() << "\n";

    if(sc == SUPERCELL_BULK)
    {
      if(s.Lattice.DiagonalOnly)
      {
        o << "    PBC=bulk Orthorhombic=yes\n";
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_BULK+TwoPowerD>(s,t);
      }
      else
      {
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_BULK>(s,t);
        o << "    PBC=bulk Orthorhombic=no\n" 
          << "    Setting Rmax = " << s.Lattice.SimulationCellRadius;
      }
    }
    else if(sc == SUPERCELL_SLAB)
    {
      if(s.Lattice.DiagonalOnly)
      {
        o << "    PBC=slab Orthorhombic=yes\n";
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_SLAB+TwoPowerD>(s,t);
      }
      else
      {
        o << "    PBC=slab Orthorhombic=no\n" 
          << "    Setting Rmax = " << s.Lattice.SimulationCellRadius;
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_SLAB>(s,t);
      }
    }
    else if(sc == SUPERCELL_WIRE)
    {
      o << "    PBC=wire Orthorhombic=NA\n";
      dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_WIRE>(s,t);
    }
    else  //open boundary condition
    {
      o << "    PBC=open Orthorhombic=NA\n";
      dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,t);
      dt->UseBoundBox=false;
    }

    dt->UseBoundBox= (sc!=SUPERCELL_OPEN);

    ostringstream p;
    p << s.getName() << "_" << t.getName();
    dt->Name=p.str();//assign the table name

    if(dt->UseBoundBox) 
      o << " Using bonding box/reduced coordinates ";
    else
      o << " using Cartesian coordinates ";
    app_log() << o.str() << endl;
    return dt;
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

} //namespace qmcplusplus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
