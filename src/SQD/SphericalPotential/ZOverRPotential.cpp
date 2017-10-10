//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include <math.h>
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf
{

/** constructor for the Nuclear Potential \f$V(r)=frac{Z}{r}\f$
 *
 * \param z the charge of the Nuclear Potential
 */
ZOverRPotential::ZOverRPotential(value_type z): Z(z)
{
  Qinfty=Z;
}

RadialPotentialBase::value_type
ZOverRPotential::evaluate(const BasisSetType& psi,
                          RadialOrbitalSet_t& V, int norb)
{
  if(!Vext)
  {
    Vext = new RadialOrbital_t(psi(0));
    integrand=new RadialOrbital_t(psi(0));
    for(int ig=0; ig < psi.m_grid->size(); ig++)
    {
      (*Vext)(ig) = -Z/psi.m_grid->r(ig);
    }
  }
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type t = (*Vext)(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += t;
      sum += pow(psi(o,ig),2);
    }
    (*integrand)(ig) = t*sum;
  }
  return integrate_RK2(*integrand);
}

int ZOverRPotential::getNumOfNodes(int n, int l)
{
  MinEigenValue = -2.0*Z*Z/static_cast<value_type>(n*n);
  return n-l-1;
}
}
