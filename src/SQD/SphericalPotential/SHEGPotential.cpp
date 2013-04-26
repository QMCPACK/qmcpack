//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <math.h>
#include "SQD/SphericalPotential/SHEGPotential.h"
#include "Numerics/RadialFunctorUtility.h"
#include "OhmmsData/ParameterSet.h"

namespace ohmmshf
{
SHEGPotential::SHEGPotential(int nel, value_type rs):
  Ntot(nel), Zext(nel), Rs(rs), Zgauss(0.0), Sgauss(1.0)
{
  MaxEigenValue=0.0;
}

RadialPotentialBase::value_type
SHEGPotential::evaluate(const BasisSetType& psi,
                        RadialOrbitalSet_t& V,
                        int norb)
{
  //intialize the external potential
  if(!Vext)
  {
    integrand=new RadialOrbital_t(psi(0));
    //evaluate Rmax from the density
    Rmax=Rs*std::pow(Ntot,1.0/3.0);
    Qinfty=Zext;//assign Zext to Qinfty
    Rcut=Rmax;//assign Rmax to Rcut
    value_type normin=-Zext/2.0/Rmax;
    value_type r2=1.0/Rmax/Rmax;
    XMLReport("  Total charge            = " << Ntot);
    XMLReport("  Total background charge = " << Zext);
    XMLReport("  Rs (density)            = " << Rs);
    XMLReport("  Rmax (background)       = " << Rmax);
    XMLReport("  Zgauss (gaussian depth) = " << Zgauss);
    XMLReport("  Sgauss (gaussian width) = " << Sgauss);
    Vext = new RadialOrbital_t(psi(0));
    for(int ig=0; ig < psi.m_grid->size(); ++ig)
    {
      value_type r=(*psi.m_grid)(ig);
      if(r<=Rmax)
        (*Vext)(ig)=normin*(3.0-(r*r)*r2)-Zgauss * exp(-r*r/Sgauss/Sgauss);
      else
        (*Vext)(ig)=-Zext/r-Zgauss * exp(-r*r/Sgauss/Sgauss) ;
      // cout << r << " " << (*Vext)(ig) << endl;
    }
  }
  for(int ig=0; ig < psi.m_grid->size(); ++ig)
  {
    value_type t = (*Vext)(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += t;
      sum += psi(o,ig)*psi(o,ig);
      //sum += pow(psi(o,ig),2);
    }
    (*integrand)(ig) = t*sum;
  }
  return integrate_RK2(*integrand);
}

/** return the number of nodes for n and l
 * @param n principal quantum number
 * @param l angular number
 * @return number of nodes
 */
int SHEGPotential::getNumOfNodes(int n, int l)
{
  MinEigenValue=-Ntot*3./2./Rmax - Zgauss;
  //double check this
  return n;
}

bool SHEGPotential::put(xmlNodePtr cur)
{
  ParameterSet params;
  params.add(Rs,"rs","double");
  params.add(Zext,"Z","double");
  params.add(Zgauss,"Zgauss","double");
  params.add(Sgauss,"Sgauss","double");
  params.put(cur);
  if(Zext<Ntot)
  {
    XMLReport("Invalid background charge Z=" << Zext << "  Overwriting by " << Ntot);
    Zext=Ntot;
  }
  return true;
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 322 $   $Date: 2005-05-07 11:55:23 -0500 (Sat, 07 May 2005) $
 * $Id: SHEGPotential.cpp 322 2005-05-07 16:55:23Z jnkim $
 ***************************************************************************/
