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
#include "simd/algorithm.hpp"
#include "Lattice/ParticleBConds3DSoa.h"
#include "Particle/SoaDistanceTableABOMP.h"
namespace qmcplusplus
{
/** Adding AsymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return index of the distance table with the name
 */
DistanceTableData* createDistanceTableABOMP(const ParticleSet& s, ParticleSet& t, int dt_type, std::ostream& description)
{
  using RealType = ParticleSet::RealType;
  enum
  {
    DIM = OHMMS_DIM
  };
  DistanceTableData* dt = 0;
  int sc = t.Lattice.SuperCellEnum;
  std::ostringstream o;
  o << "  Distance table for dissimilar particles (A-B):" << std::endl;
  o << "    source: " << s.getName() << "  target: " << t.getName() << std::endl;
  o << "    Using structure-of-arrays (SoA) data layout and OpenMP offload" << std::endl;

  if (sc == SUPERCELL_BULK)
  {
    if (s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic periodic cell in 3D." << std::endl;
      dt = new SoaDistanceTableABOMP<RealType, DIM, PPPO + SOA_OFFSET>(s, t);
    }
    else
    {
      if (s.Lattice.WignerSeitzRadius > s.Lattice.SimulationCellRadius)
      {
        o << "    Distance computations use general periodic cell in 3D with corner image checks." << std::endl;
        dt = new SoaDistanceTableABOMP<RealType, DIM, PPPG + SOA_OFFSET>(s, t);
      }
      else
      {
        o << "    Distance computations use general periodic cell in 3D without corner image checks." << std::endl;
        dt = new SoaDistanceTableABOMP<RealType, DIM, PPPS + SOA_OFFSET>(s, t);
      }
    }
  }
  else if (sc == SUPERCELL_SLAB)
  {
    if (s.Lattice.DiagonalOnly)
    {
      o << "    Distance computations use orthorhombic code for periodic cell in 2D." << std::endl;
      dt = new SoaDistanceTableABOMP<RealType, DIM, PPNO + SOA_OFFSET>(s, t);
    }
    else
    {
      if (s.Lattice.WignerSeitzRadius > s.Lattice.SimulationCellRadius)
      {
        o << "    Distance computations use general periodic cell in 2D with corner image checks." << std::endl;
        dt = new SoaDistanceTableABOMP<RealType, DIM, PPNG + SOA_OFFSET>(s, t);
      }
      else
      {
        o << "    Distance computations use general periodic cell in 2D without corner image checks." << std::endl;
        dt = new SoaDistanceTableABOMP<RealType, DIM, PPNS + SOA_OFFSET>(s, t);
      }
    }
  }
  else if (sc == SUPERCELL_WIRE)
  {
    o << "    Distance computations use periodic cell in one dimension." << std::endl;
    dt = new SoaDistanceTableABOMP<RealType, DIM, SUPERCELL_WIRE + SOA_OFFSET>(s, t);
  }
  else //open boundary condition
  {
    o << "    Distance computations use open boundary conditions in 3D." << std::endl;
    dt = new SoaDistanceTableABOMP<RealType, DIM, SUPERCELL_OPEN + SOA_OFFSET>(s, t);
  }

  //set dt properties
  dt->DTType = DT_SOA;
  std::ostringstream p;
  p << s.getName() << "_" << t.getName();
  dt->setName(p.str()); //assign the table name

  description << o.str() << std::endl;

  return dt;
}


} //namespace qmcplusplus
