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
    DistanceTableData* dt=0;
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
    {
      app_log() << "  Sym Distance table specialized for an open cell ";
      //dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,s);
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,s);
      dt->UseBoundBox=false;
    }
    else
    {
      if(s.Lattice.DiagonalOnly)
      {
        app_log() << "  Sym Distance table specialized for an Orthorhombic cell ";
        //dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,s);
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_BULK+TwoPowerD>(s,s);
      }
      else
      {
        //dt = new SymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,s);
        dt = new SymmetricDTD<RealType,DIM,SUPERCELL_BULK>(s,s);
        //dt->setRmax(s.Lattice.SimulationCellRadius);
        app_log() << "  Sym Distance table specialized for a generic cell ";
        app_log() << "  Setting Rmax = " << s.Lattice.SimulationCellRadius <<endl;
      }
      dt->UseBoundBox=true;
    }
    ostringstream o;
    o << s.getName() << "_" << s.getName();
    dt->Name=o.str();//assign the table name
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
    if(s.Lattice.SuperCellEnum == SUPERCELL_OPEN)
    {
      app_log() << "  Asymm Distance table specialized for an open cell ";
      //dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_OPEN> >(s,t);
      dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,t);
      dt->UseBoundBox=false;
    }
    else 
    {
      if(s.Lattice.DiagonalOnly)
      {
        app_log() << "  Asymm Distance table specialized for an Orthorhombic cell ";
        //dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK+TwoPowerD> >(s,t);
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_BULK+TwoPowerD>(s,t);
      }
      else
      {
        //dt = new AsymmetricDTD<DTD_BConds<RealType,DIM,SUPERCELL_BULK> >(s,t);
        dt = new AsymmetricDTD<RealType,DIM,SUPERCELL_BULK>(s,t);
        //dt->setRmax(s.Lattice.SimulationCellRadius);
        app_log() << "  Asymm Distance table specialized for a generic cell ";
        app_log() << "  Setting Rmax = " << s.Lattice.SimulationCellRadius <<endl;
      }
      dt->UseBoundBox=true;
    }
    ostringstream o;
    o << s.getName() << "_" << t.getName();
    dt->Name=o.str();//assign the table name
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
