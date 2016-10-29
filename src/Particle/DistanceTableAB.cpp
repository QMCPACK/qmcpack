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
#include "Particle/AsymmetricDistanceTableData.h"
#include "simd/algorithm.hpp"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/SoaDistanceTableBA.h"
namespace qmcplusplus
{

/** Adding AsymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
DistanceTableData* createDistanceTable(const ParticleSet& s, ParticleSet& t, int dt_type)
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
        if(dt_type == DT_SOA)
        {
          o << "  Using SoaDistanceTableBA<T,D,PPPG> of SoA layout " << PPPG << std::endl;
          dt = new  SoaDistanceTableBA<RealType,DIM,PPPG+SOA_OFFSET>(s,t);
        }
        else
        {
          o << " Using AsymmetricDTD<T,D,PPPG> " << PPPG << std::endl;
          dt = new AsymmetricDTD<RealType,DIM,PPPG>(s,t);
        }
      }
      else
      {
        if(dt_type == DT_SOA)
        {
          o << " Using SoaDistanceTableBA<T,D,PPPS> with SoA layout " << PPPS << std::endl;
          dt = new SoaDistanceTableBA<RealType,DIM,PPPS+SOA_OFFSET>(s,t);
        }
        else
        {
          o << " Using AsymmetricDTD<T,D,PPPS> " << PPPS << std::endl;
          dt = new AsymmetricDTD<RealType,DIM,PPPS>(s,t);
        }
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
  dt->DTType=dt_type;
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


} //namespace qmcplusplus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
