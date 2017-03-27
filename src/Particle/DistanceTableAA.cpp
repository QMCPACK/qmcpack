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
  bool useSoA=(dt_type == DT_SOA || dt_type == DT_SOA_PREFERRED);
  o << "  Distance table for AA: source/target = " << s.getName() << " useSoA =" << useSoA << "\n";
  if(sc == SUPERCELL_BULK)
  {
    if(s.Lattice.DiagonalOnly)
    {
      if(useSoA)
      {
        o << "    PBC=bulk Orthorhombic=yes  Using SoaDistanceTableAA<T,D,PPPO> of SoA layout " << PPPO << std::endl;
        dt = new SoaDistanceTableAA<RealType,DIM,PPPO+SOA_OFFSET>(s);
      }
      else
      {
        o << "    PBC=bulk Orthorhombic=yes  Using SymmetricDTD<T,DIM,PPPO> " << PPPO << std::endl; 
        dt = new SymmetricDTD<RealType,DIM,PPPO>(s,s);
      }
    }
    else
    {
      o << "    PBC=bulk Orthorhombic=no";
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        if(useSoA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPPG> of SoA layout " << PPPG << std::endl;
          dt = new SoaDistanceTableAA<RealType,DIM,PPPG+SOA_OFFSET>(s);
        }
        else
        {
          o << "  Using SymmetricDTD<T,D,PPPG> " << PPPG << std::endl;
          dt = new  SymmetricDTD<RealType,DIM,PPPG>(s,s);
        }
      }
      else
      {
        if(useSoA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPPS> of SoA layout " << PPPS << std::endl;
          dt = new SoaDistanceTableAA<RealType,DIM,PPPS+SOA_OFFSET>(s);
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
      if(useSoA)
      {
        o << "    PBC=slab Orthorhombic=yes  Using SoaDistanceTableAA<T,D,PPNO> of SoA layout " << PPNO << std::endl;
        dt = new SoaDistanceTableAA<RealType,DIM,PPNO+SOA_OFFSET>(s);
      }
      else
      {
        o << "    PBC=slab Orthorhombic=yes  Using SymmetricDTD<T,D,PPNO> " << PPNO << std::endl;
        dt = new SymmetricDTD<RealType,DIM,PPNO>(s,s);
      }
    }
    else
    {
      o << "    PBC=slab Orthorhombic=no";
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        if(useSoA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPNX> of SoA layout " << PPNX << std::endl;
          dt = new SoaDistanceTableAA<RealType,DIM,PPNX+SOA_OFFSET>(s);
        }
        else
        {
          o << "  Using SymmetricDTD<T,D,PPNX> " << PPNX << std::endl;
          dt = new SymmetricDTD<RealType,DIM,PPNX>(s,s);
        }
      }
      else
      {
        if(useSoA)
        {
          o << "  Using SoaDistanceTableAA<T,D,PPNS> of SoA layout " << PPNS << std::endl;
          dt = new SoaDistanceTableAA<RealType,DIM,PPNS+SOA_OFFSET>(s);
        }
        else
        {
          o << "  Using SymmetricDTD<T,D,PPNS> " << PPNS << std::endl;
          dt = new SymmetricDTD<RealType,DIM,PPNS>(s,s);
        }
      }
    }
  }
  else if(sc == SUPERCELL_WIRE)
  {
    if(useSoA)
    {
      o << "    PBC=wire Orthorhombic=NA  Using SoA layout\n";
      dt = new SoaDistanceTableAA<RealType,DIM,SUPERCELL_WIRE+SOA_OFFSET>(s);
    }
    else
    {
      o << "    PBC=wire Orthorhombic=NA\n";
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_WIRE>(s,s);
    }
  }
  else  //open boundary condition
  {
    if(useSoA)
    {
      o << "    PBC=open Orthorhombic=NA  Using SoA layout\n";
      dt = new SoaDistanceTableAA<RealType,DIM,SUPERCELL_OPEN+SOA_OFFSET>(s);
    }
    else
    {
      o << "    PBC=open Orthorhombic=NA\n";
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,s);
    }
  }

  //set dt properties
  dt->CellType=sc;
  dt->DTType=(useSoA)? DT_SOA: DT_AOS;
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
