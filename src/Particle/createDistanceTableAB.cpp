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


#include "Particle/createDistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Lattice/ParticleBConds.h"
#ifndef ENABLE_SOA
#include "Particle/AsymmetricDistanceTableData.h"
#endif
#include "simd/algorithm.hpp"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/SoaDistanceTableAB.h"
namespace qmcplusplus
{
/** Adding AsymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
DistanceTableData* createDistanceTableAB(const ParticleSet& s, ParticleSet& t, int dt_type, std::ostream& description)
{
  using RealType = ParticleSet::RealType;
  enum
  {
    DIM = OHMMS_DIM
  };
  DistanceTableData* dt = 0;
  //int sc=s.Lattice.SuperCellEnum;
  int sc = t.Lattice.SuperCellEnum;
  std::ostringstream o;
  bool useSoA = (dt_type == DT_SOA || dt_type == DT_SOA_PREFERRED);
  o << "  Distance table for dissimilar particles (A-B):" << std::endl;
  o << "    source: " << s.getName() << "  target: " << t.getName() << std::endl;
  if (useSoA)
  {
    o << "    Using structure-of-arrays (SoA) data layout" << std::endl;
  }
  else
  {
#ifdef ENABLE_SOA
    APP_ABORT("createDistanceTable (AB). Using array-of-structure (AoS) data layout is no longer supported in builds "
              "with ENABLE_SOA=1.");
#else
    o << "    Using array-of-structure (AoS) data layout (less efficient than SoA)" << std::endl;
#endif
  }

  if (sc == SUPERCELL_BULK)
  {
    if (s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic periodic cell in 3D." << std::endl;
      if (useSoA)
      {
        dt = new SoaDistanceTableAB<RealType, DIM, PPPO + SOA_OFFSET>(s, t);
      }
      else
      {
#ifndef ENABLE_SOA
        dt = new AsymmetricDTD<RealType, DIM, PPPO>(s, t);
#endif
      }
    }
    else
    {
      if (s.Lattice.WignerSeitzRadius > s.Lattice.SimulationCellRadius)
      {
        o << "    Distance computations use general periodic cell in 3D with corner image checks." << std::endl;
        if (useSoA)
        {
          dt = new SoaDistanceTableAB<RealType, DIM, PPPG + SOA_OFFSET>(s, t);
        }
        else
        {
#ifndef ENABLE_SOA
          dt = new AsymmetricDTD<RealType, DIM, PPPG>(s, t);
#endif
        }
      }
      else
      {
        o << "    Distance computations use general periodic cell in 3D without corner image checks." << std::endl;
        if (useSoA)
        {
          dt = new SoaDistanceTableAB<RealType, DIM, PPPS + SOA_OFFSET>(s, t);
        }
        else
        {
#ifndef ENABLE_SOA
          dt = new AsymmetricDTD<RealType, DIM, PPPS>(s, t);
#endif
        }
      }
    }
  }
  else if (sc == SUPERCELL_SLAB)
  {
    if (s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic code for periodic cell in 2D." << std::endl;
      if (useSoA)
      {
        dt = new SoaDistanceTableAB<RealType, DIM, PPNO + SOA_OFFSET>(s, t);
      }
      else
      {
#ifndef ENABLE_SOA
        dt = new AsymmetricDTD<RealType, DIM, PPNO>(s, t);
#endif
      }
    }
    else
    {
      if (s.Lattice.WignerSeitzRadius > s.Lattice.SimulationCellRadius)
      {
        if (useSoA)
        {
          o << "    Distance computations use general periodic cell in 2D with corner image checks." << std::endl;
          dt = new SoaDistanceTableAB<RealType, DIM, PPNG + SOA_OFFSET>(s, t);
        }
        else
        {
          o << "    Distance computations use general periodic cell in 2D with all surrounding image checks."
            << std::endl;
#ifndef ENABLE_SOA
          dt = new AsymmetricDTD<RealType, DIM, PPNX>(s, t);
#endif
        }
      }
      else
      {
        o << "    Distance computations use general periodic cell in 2D without corner image checks." << std::endl;
        if (useSoA)
        {
          dt = new SoaDistanceTableAB<RealType, DIM, PPNS + SOA_OFFSET>(s, t);
        }
        else
        {
#ifndef ENABLE_SOA
          dt = new AsymmetricDTD<RealType, DIM, PPNS>(s, t);
#endif
        }
      }
    }
  }
  else if (sc == SUPERCELL_WIRE)
  {
    o << "    Distance computations use periodic cell in one dimension." << std::endl;
    if (useSoA)
    {
      dt = new SoaDistanceTableAB<RealType, DIM, SUPERCELL_WIRE + SOA_OFFSET>(s, t);
    }
    else
    {
#ifndef ENABLE_SOA
      dt = new AsymmetricDTD<RealType, DIM, SUPERCELL_WIRE>(s, t);
#endif
    }
  }
  else //open boundary condition
  {
    o << "    Distance computations use open boundary conditions in 3D." << std::endl;
    if (useSoA)
    {
      dt = new SoaDistanceTableAB<RealType, DIM, SUPERCELL_OPEN + SOA_OFFSET>(s, t);
    }
    else
    {
#ifndef ENABLE_SOA
      dt = new AsymmetricDTD<RealType, DIM, SUPERCELL_OPEN>(s, t);
#endif
    }
  }

  //set dt properties
  dt->DTType = (useSoA) ? DT_SOA : DT_AOS;
  std::ostringstream p;
  p << s.getName() << "_" << t.getName();
  dt->setName(p.str()); //assign the table name

  description << o.str() << std::endl;

  return dt;
}


} //namespace qmcplusplus
