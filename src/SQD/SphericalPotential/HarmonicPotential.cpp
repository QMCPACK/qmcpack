//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <math.h>
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf
{
/*!
 * \param omega The force constant of a harmonic oscillator.
 */
HarmonicPotential::HarmonicPotential(value_type omega): Omega(omega)
{
  MinEigenValue = 0.0;
}

RadialPotentialBase::value_type
HarmonicPotential::evaluate(const BasisSetType& psi,
                            RadialOrbitalSet_t& V,
                            int norb)
{
  if(!Vext)
  {
    Vext = new RadialOrbital_t(psi(0));
    integrand=new RadialOrbital_t(psi(0));
    value_type omega2=Omega*Omega*0.5;
    for(int ig=0; ig < psi.m_grid->size(); ig++)
    {
      (*Vext)(ig)= omega2*psi.m_grid->r(ig)*psi.m_grid->r(ig);
    }
  }
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type sum = 0.0;
    value_type v = (*Vext)(ig);
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += v;
      sum += pow(psi(o,ig),2);
    }
    (*integrand)(ig) = v*sum;
  }
  return integrate_RK2(*integrand);
}

int HarmonicPotential::getNumOfNodes(int n, int l)
{
  MaxEigenValue = (2.0*static_cast<value_type>(2*n+l)+3.0)*Omega;
  return n;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
