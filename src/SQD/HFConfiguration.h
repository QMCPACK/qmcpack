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
    
    


#ifndef OHMMS_MAIN_HFCONFIGURATION_H
#define OHMMS_MAIN_HFCONFIGURATION_H
/*!\mainpage AtomicHF
 *\author Jeongnim Kim
 *\author Jordan Vincent
 *
 *\section intro_sec Introduction
 *
 *Package to solve the Hartree-Fock equations for a centered (atomic-like)
 *spherically symmetric potential using the Numerov algorithm.
 *
 *\section main_sec Schrodinger Equation
 *
 *The single-particle Schrodinger Equation involving the mass
 \f$m\f$ and an external potential \f$V(\bf r)\f$,
 \f[
 \left[-\frac{1}{2m}\nabla^2 + V(\bf r)\right]\Psi(\bf r) = E \Psi(\bf r),
 \f]
 *is reduced to a Radial Schrodinger Equation
 \f[
 \left\{
 \frac{1}{2m}\frac{d^{2}}{dr^{2}}
 +\left[
 E-V(r)-\frac{l(l+1)}{2mr^{2}}
 \right]
 \right\} u_{l}(r)=0,
 \f]
 *for a spherically symmetric potentials \f$V(r)\f$ and the wave functions are
 \f[
 \Psi({\bf r}) = \sum_{l=0}^{\infty}\sum_{m=-l}^{l} A_{lm}
 \frac{u_{l}(r)}{r}Y_{l}^{m}(\theta,\phi).
 \f]
 *AtomicHF uses the Numerov algorithm to solve the Radial Schrodinger Equation.
 *Several external radial potentials are implemented which are defined
 *in \ref RadialPotential.
 *
 *Within Hartree-Fock approximation, the Schrodinger Equation is generalized
 *to include electron-electron Coulomb and exchange interactions (Hartree
 *and exchange terms). Hartree and Exchange terms are implemented in
 *HartreePotential and ExchangePotential, while the external potential and
 *kinetic energy terms are handled by ZOverRPotential,
 *HarmonicPotential, SJPseudopotential specific for the physical model.
 *
 *\section numerov_sec Numemov algorithm
 *
 *Numerov algorithm (B. Numerov, Publ. de l'observ. astrophsique
 *central de Russie, vol. 2, p 188, 1933.) is based on a two-step, fifth order
 *predictor-corrector method for a second-order ordinary differential equation
 *of type \f[\frac{d^2y}{dx^2} + k^2(x)y = S(x).\f]
 *Utilizing Taylor expansions for \f$x\f$, the solution is obtained by
 *recursively integrating forward in x as
 \f[
 (1+\frac{h^2}{12}k^2_{n+2})y_{n+2} -
 2(1-\frac{5h^2}{12}k^2_{n+1})y_{n+1} + (1+\frac{h^2}{12}k^2_n)y_n =
 \frac{h^2}{12}(S_{n+2}+10S_{n+1}+S_n) + \mathcal{O}(h^6).
 \f]
 *Numerov algorithm uses a uniform linear grid.
 Therefore, a transformation function is used to handle
 - transformations of variables between the physical grid to Numerov grid
 - cusp conditions or boundary conditions
 - range of the eigen values
 *
 *Any grid type can be used as far as the source functor can transform
 *the original grid variable to a uniform grid on \f$x\f$.
 *
 *\section hf_sec Hartree-Fock method
 *
 *Short description of HF method here.
 *
 *\section app_sec How to run AtomicHF code
 *
 */
#include "Configuration.h"
#include "SQD/YlmRnlSet.h"
/**@namespace ohmmshf
 *@brief  Define basic data types for the applications.
 *
 * In order to reduce complier-time complexity and to enable switching
 * between CPLUSPLUS libraries for array and expression template,
 * basic data types are defined.
 */
namespace ohmmshf
{

/**Trait class for spherical orbitals
 *
 *This class defines the data types to represent the spherical potential and wave functions.
 */
struct SphericalOrbitalTraits
{

  ///@typedef the basis set type
  typedef YlmRnlSet<double>                BasisSetType;

  ///@typedef the type of value: double
  typedef BasisSetType::value_type         value_type;

  ///@typedef the type of radial grids
  typedef BasisSetType::RadialGrid_t       RadialGrid_t;

  ///@typedef the type of radial orbitals
  typedef BasisSetType::RadialOrbital_t    RadialOrbital_t;

  ///@typedef the type of radial orbital set
  typedef BasisSetType::RadialOrbitalSet_t RadialOrbitalSet_t;
};
}

#endif

