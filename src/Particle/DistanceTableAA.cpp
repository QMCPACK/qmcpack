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
  o << "  Distance table for similar particles (A-A):" << std::endl;
  o << "    source/target: " << s.getName() << std::endl;
  if (useSoA) {
    o << "    Using structure-of-arrays (SoA) data layout" << std::endl;
  } else {
    o << "    Using array-of-structure (AoS) data layout (less efficient than SoA)" << std::endl;
  }

  if(sc == SUPERCELL_BULK)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic periodic cell in 3D." << std::endl;
      if(useSoA)
      {
        dt = new SoaDistanceTableAA<RealType,DIM,PPPO+SOA_OFFSET>(s);
      }
      else
      {
        dt = new SymmetricDTD<RealType,DIM,PPPO>(s,s);
      }
    }
    else
    {
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        o << "    Distance computations use general periodic cell in 3D with corner image checks." << std::endl;
        if(useSoA)
        {
          dt = new SoaDistanceTableAA<RealType,DIM,PPPG+SOA_OFFSET>(s);
        }
        else
        {
          dt = new  SymmetricDTD<RealType,DIM,PPPG>(s,s);
        }
      }
      else
      {
        o << "    Distance computations use general periodic cell in 3D without corner image checks." << std::endl;
        if(useSoA)
        {
          dt = new SoaDistanceTableAA<RealType,DIM,PPPS+SOA_OFFSET>(s);
        }
        else
        {
          dt = new  SymmetricDTD<RealType,DIM,PPPS>(s,s);
        }
      }
      o << "    Setting Rmax = " << s.Lattice.SimulationCellRadius << std::endl;
    }
  }
  else if(sc == SUPERCELL_SLAB)
  {
    if(s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic code for periodic cell in 2D." << std::endl;
      if(useSoA)
      {
        dt = new SoaDistanceTableAA<RealType,DIM,PPNO+SOA_OFFSET>(s);
      }
      else
      {
        dt = new SymmetricDTD<RealType,DIM,PPNO>(s,s);
      }
    }
    else
    {
      if(s.Lattice.WignerSeitzRadius>s.Lattice.SimulationCellRadius)
      {
        if(useSoA)
        {
          o << "    Distance computations use general periodic cell in 2D with corner image checks." << std::endl;
          dt = new SoaDistanceTableAA<RealType,DIM,PPNG+SOA_OFFSET>(s);
        }
        else
        {
          o << "    Distance computations use general periodic cell in 2D with all surrounding image checks." << std::endl;
          dt = new SymmetricDTD<RealType,DIM,PPNX>(s,s);
        }
      }
      else
      {
        o << "    Distance computations use general periodic cell in 2D without corner image checks." << std::endl;
        if(useSoA)
        {
          dt = new SoaDistanceTableAA<RealType,DIM,PPNS+SOA_OFFSET>(s);
        }
        else
        {
          dt = new SymmetricDTD<RealType,DIM,PPNS>(s,s);
        }
      }
    }
  }
  else if(sc == SUPERCELL_WIRE)
  {
    o << "    Distance computations use periodic cell in one dimension." << std::endl;
    if(useSoA)
    {
      dt = new SoaDistanceTableAA<RealType,DIM,SUPERCELL_WIRE+SOA_OFFSET>(s);
    }
    else
    {
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_WIRE>(s,s);
    }
  }
  else  //open boundary condition
  {
    o << "    Distance computations use open boundary conditions in 3D." << std::endl;
    if(useSoA)
    {
      dt = new SoaDistanceTableAA<RealType,DIM,SUPERCELL_OPEN+SOA_OFFSET>(s);
    }
    else
    {
      dt = new SymmetricDTD<RealType,DIM,SUPERCELL_OPEN>(s,s);
    }
  }


  //set dt properties
  dt->CellType=sc;
  dt->DTType=(useSoA)? DT_SOA: DT_AOS;
  std::ostringstream p;
  p << s.getName() << "_" << s.getName();
  dt->Name=p.str();//assign the table name

  if(omp_get_thread_num()==0) 
  {
    app_log() << std::endl;
    app_log() << o.str() << std::endl;
    app_log().flush();
  }

  return dt;
}

} //namespace qmcplusplus

