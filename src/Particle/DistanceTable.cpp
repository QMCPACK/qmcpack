//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
  std::ostringstream o;
  o << "  Distance table for AA: source/target = " << s.getName() << "\n";
  if(sc == SUPERCELL_BULK)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    PBC=bulk Orthorhombic=yes Using SymmetricDTD<T,DIM,PPPO> " << PPPO << std::endl;
      dt = new SymmetricDTD<RealType,DIM,PPPO>(s,s);
    }
    else
    {
      o << "    PBC=bulk Orthorhombic=no ";
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        o << "  Using SymmetricDTD<T,D,PPPG> " << PPPG << std::endl;
        dt = new  SymmetricDTD<RealType,DIM,PPPG>(s,s);
        //o << "    PBC=bulk Orthorhombic=no SymmetricDTD<T,DIM,PPPX> " << PPPX << std::endl;
        //dt = new  SymmetricDTD<RealType,DIM,PPPX>(s,s);
      }
      else
      {
        o << "  Using SymmetricDTD<T,D,PPPS> " << PPPS << std::endl;
        dt = new  SymmetricDTD<RealType,DIM,PPPS>(s,s);
      }
      o << "\n    Setting Rmax = " << s.Lattice.SimulationCellRadius;
    }
  }
  else if(sc == SUPERCELL_SLAB)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    PBC=slab Orthorhombic=yes Using SymmetricDTD<T,D,PPNO> " << PPNO << std::endl;
      dt = new SymmetricDTD<RealType,DIM,PPNO>(s,s);
    }
    else
    {
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        o << "    PBC=slab Orthorhombic=no Using SymmetricDTD<T,D,PPNX> " << PPNX << std::endl;
        dt = new SymmetricDTD<RealType,DIM,PPNX>(s,s);
      }
      else
      {
        o << "    PBC=slab Orthorhombic=no Using SymmetricDTD<T,D,PPNS> " << PPNS << std::endl;
        dt = new SymmetricDTD<RealType,DIM,PPNS>(s,s);
      }
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
  dt->CellType=sc;
  std::ostringstream p;
  p << s.getName() << "_" << s.getName();
  dt->Name=p.str();//assign the table name
  if(sc != SUPERCELL_OPEN)
    o << " Using bounding box/reduced coordinates with ";
  else
    o << " using Cartesian coordinates with ";
  if(omp_get_thread_num()==0) 
  {
    app_log() << o.str() << std::endl;
    app_log().flush();
  }
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
  std::ostringstream o;
  o << "  Distance table for AB: source = " << s.getName() << " target = " << t.getName() << "\n";
  if(sc == SUPERCELL_BULK)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    PBC=bulk Orthorhombic=yes Using AsymmetricDTD<T,D,PPPO> " << PPPO << std::endl;
      dt = new AsymmetricDTD<RealType,DIM,PPPO>(s,t);
    }
    else
    {
      o << "    PBC=bulk Orthorhombic=no ";
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        o << " Using AsymmetricDTD<T,D,PPPG> " << PPPG << std::endl;
        dt = new AsymmetricDTD<RealType,DIM,PPPG>(s,t);
      }
      else
      {
        o << " Using AsymmetricDTD<T,D,PPPS> " << PPPS << std::endl;
        dt = new AsymmetricDTD<RealType,DIM,PPPS>(s,t);
      }
      o << "    Setting Rmax = " << s.Lattice.SimulationCellRadius;
    }
  }
  else if(sc == SUPERCELL_SLAB)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    PBC=slab Orthorhombic=yes Using AsymmetricDTD<T,D,PPNO> " << PPNO << std::endl;
      dt = new AsymmetricDTD<RealType,DIM,PPNO>(s,t);
    }
    else
    {
      o << "    PBC=slab Orthorhombic=no ";
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        o << " Using AsymmetricDTD<T,DIM,PPNX> " << PPNX << std::endl;
        dt = new AsymmetricDTD<RealType,DIM,PPNX>(s,t);
      }
      else
      {
        o << " Using AsymmetricDTD<T,DIM,PPNS> " << PPNS << std::endl;
        dt = new AsymmetricDTD<RealType,DIM,PPNS>(s,t);
      }
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
  }

  dt->CellType=sc;
  std::ostringstream p;
  p << s.getName() << "_" << t.getName();
  dt->Name=p.str();//assign the table name
  if(sc!=SUPERCELL_OPEN)
    o << " Using bonding box/reduced coordinates ";
  else
    o << " using Cartesian coordinates ";
  if(omp_get_thread_num()==0) 
  {
    app_log() << o.str() << std::endl;
    app_log().flush();
  }
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
