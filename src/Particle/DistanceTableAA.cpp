//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Lattice/ParticleBConds.h"
#include "Particle/SymmetricDistanceTableData.h"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/SoaDistanceTableAA.h"
namespace qmcplusplus
{

/** Adding SymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
DistanceTableData* createDistanceTable(ParticleSet& s, int dt_type)
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
        if(dt_type == DT_SOA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPPG> of SoA layout " << PPPG << std::endl;
          dt = new  SoaDistanceTableAA<RealType,DIM,PPPG+SOA_OFFSET>(s);
        }
        else
        {
          o << "  Using SymmetricDTD<T,D,PPPG> " << PPPG << std::endl;
          dt = new  SymmetricDTD<RealType,DIM,PPPG>(s,s);
        }
      }
      else
      {
        if(dt_type == DT_SOA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPPS> of SoA layout " << PPPS << std::endl;
          dt = new  SoaDistanceTableAA<RealType,DIM,PPPS+SOA_OFFSET>(s);
        }
        else
        {
          o << "  Using SymmetricDTD<T,D,PPPS> " << PPPS << std::endl;
          dt = new  SymmetricDTD<RealType,DIM,PPPS>(s,s);
        }
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
  dt->DTType=dt_type;
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

} //namespace qmcplusplus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
