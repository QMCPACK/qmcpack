//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <math.h>
#include "SQD/SphericalPotential/StepPotential.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf
{
/*!
 *Add a step function.
 */
StepPotential::StepPotential()
{
  Rseg.resize(2);
  Vseg.resize(2);
  Vseg[0] = -10;
  Vseg[1] = 0;
  Rseg[0] = 4.0;
  Rseg[1] = 10;
  MinEigenValue = -10.0;
  MaxEigenValue = 0.0;
}

RadialPotentialBase::value_type
StepPotential::evaluate(const BasisSetType& psi,
                        RadialOrbitalSet_t& V,
                        int norb)
{
  if(!Vext)
  {
    integrand=new RadialOrbital_t(psi(0));
    Vext = new RadialOrbital_t(psi(0));
    int ilow=0, ihi=1;
    for(int ig=0; ig < psi.m_grid->size(); ig++)
    {
      (*Vext)(ig) = Vseg[ilow];
      if(psi.m_grid->r(ig)> Rseg[ilow])
      {
        ilow++;
        ihi++;
      }
    }
  }
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type v = (*Vext)(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += v;
      sum += pow(psi(o,ig),2);
    }
    (*integrand)(ig) = v*sum;
  }
  return integrate_RK2(*integrand);
}

int StepPotential::getNumOfNodes(int n, int l)
{
  return n;
}

bool StepPotential::put(xmlNodePtr cur)
{
  cur = cur->xmlChildrenNode;
  std::vector<double> r0, v0;
  double vmin=1000;
  double vmax=-1000;
  while(cur!= NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "Region")
    {
      xmlAttrPtr att = cur->properties;
      double rmin, rmax, v;
      while(att != NULL)
      {
        std::string aname((const char*)(att->name));
        const char* vname= (const char*)(att->children->content);
        if(aname == "rmin")
        {
          rmin = atof(vname);
        }
        else
          if(aname == "rmax")
          {
            rmax = atof(vname);
          }
          else
            if(aname == "value")
            {
              v = atof(vname);
              vmin = std::min(v,vmin);
              vmax = std::min(v,vmax);
            }
        att = att->next;
      }
      r0.push_back(rmax);
      v0.push_back(v);
    }
    cur = cur->next;
  }
  if(r0.size())
  {
    MinEigenValue = vmin;
    //This needs to be reconsidered
    //MaxEigenValue = vmax;
    Vseg=v0;
    Rseg=r0;
  }
  return true;
}
}
