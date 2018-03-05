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
namespace qmcplusplus
{

/** Adding SymmetricDTD to the list, e.g., el-el distance table
 *\param s source/target particle set
 *\return DistanceTableData*
 */
DistanceTableData* DistanceTable::add(ParticleSet& s, int dt_type)
{
  int tid=s.addTable(s,dt_type);
  return s.DistTables[tid];
}

/** Adding AsymmetricDTD to the list, e.g., el-nuclei distance table
 *\param s source particle set
 *\param t target particle set
 *\return DistanceTableData*
 */
DistanceTableData* DistanceTable::add(const ParticleSet& s, ParticleSet& t, int dt_type)
{
  int tid=t.addTable(s,dt_type);
  return t.DistTables[tid];
}

} //namespace qmcplusplus
