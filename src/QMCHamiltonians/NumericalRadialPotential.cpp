//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/NumericalRadialPotential.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

NumericalRadialPotential::NumericalRadialPotential(ParticleSet& center): sourcePtcl(center),
  VofR(0)
{
  IsPairPotential=true;
  d_table = DistanceTable::add(center);
  Centers = center.getTotalNum();
}

NumericalRadialPotential::NumericalRadialPotential(ParticleSet& center, ParticleSet& visitor):
  sourcePtcl(center), VofR(0)
{
  IsPairPotential=false;
  d_table = DistanceTable::add(center,visitor);
  Centers = center.getTotalNum();
}

///destructor
NumericalRadialPotential::~NumericalRadialPotential()
{
  if(VofR)
    delete VofR;
}

void NumericalRadialPotential::resetTargetParticleSet(ParticleSet& P)
{
  if(IsPairPotential)
    d_table = DistanceTable::add(P);
  else
    d_table = DistanceTable::add(sourcePtcl,P);
}

QMCHamiltonianBase::Return_t
NumericalRadialPotential:: evaluate(ParticleSet& P)
{
  Value=0.0;
  for(int iat=0; iat<Centers; iat++)
  {
    RealType e = 0.0;
    for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++)
    {
      e += VofR->splint(d_table->r(nn));
    }
  }
  return Value;
}

/** read data of the numerical radial potential
 *
 * <pairpot type="numericaldata" source="e">
 *   <grid type="linear" ri="0" rf="10" npts="100"/> <!-- optional -->
 *   <data ri="0" rf="10" size="100"> <!-- overwrites grid rf,ri and npts -->
 *   data
 *   </data>
 * </pairpot>
 *
 *
 */
bool NumericalRadialPotential::put(xmlNodePtr cur)
{
  //cannot do it again
  if(VofR)
    return false;
  //define a grid spec
  RealType ri = 0.0;
  RealType rf = 0.0;
  int npts = -1;
  OhmmsAttributeSet radAttrib;
  radAttrib.add(npts,"npts");
  radAttrib.add(npts,"size");
  radAttrib.add(ri,"ri");
  radAttrib.add(rf,"rf");
  //read xml elemet
  std::vector<RealType> vr;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "grid")
    {
      radAttrib.put(cur);
    }
    else
      if(cname == "data")
      {
        radAttrib.put(cur);
        putContent(vr,cur);
      }
    cur=cur->next;
  }
  //need to die if the grid data
  if(npts<0 || rf<=ri)
    return false;
  GridType *agrid = new LinearGrid<RealType>;
  agrid->set(ri,rf,npts);
  VofR=new RadialPotentialType(agrid,vr);
  VofR->spline();
  return true;
}

bool NumericalRadialPotential::get(std::ostream& os) const
{
  os << "NumericalRadial potential: " << std::endl;
  return true;
}

QMCHamiltonianBase*
NumericalRadialPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  NumericalRadialPotential* apot=0;
  if(IsPairPotential)
    apot=new NumericalRadialPotential(qp);
  else
    apot=new NumericalRadialPotential(sourcePtcl,qp);
  apot->VofR = new RadialPotentialType(*VofR);
  return apot;
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: NumericalRadialPotential.h 1581 2007-01-04 16:02:14Z jnkim $
 ***************************************************************************/

