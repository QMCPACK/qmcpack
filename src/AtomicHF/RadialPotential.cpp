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
#include "AtomicHF/RadialPotential.h"
#include "AtomicHF/Clebsch_Gordan.h"
#include "Numerics/RadialFunctorUtility.h"
using namespace ohmmshf;

/*!
 * \fn ZOverRFunctor::ZOverRFunctor(value_type z)
 * \param z the charge of the Nuclear Potential
 * \brief The Constructor for the Nuclear Potential
 *
 */

ZOverRFunctor::ZOverRFunctor(value_type z): Z(z) { }

/*! RadialPotentialBase::value_type
  PseudoPotential::evaluate(const HFAtomicOrbitals& psi,
  RadialOrbitalSet_t& V, int norb)
  * \param psi the wavefunction
  * \param V the potential
  * \param norb the number of orbitals
  * \return The sum of the Nuclear Potential matrix
  elements (energy): \f[
  \langle i|V(r)|i \rangle = \sum_{k=0}^{N_{orb}}
  \int_0^{\infty} dr\psi_k^*(r)V(r)\psi_k(r) \f]
  * \brief Calculates and assigns the values of the Nuclear
  Potential for each orbital.  The Nuclear Potential: \f[
  V_{Nuclear}(r) = -\frac{Z}{r} \f]
  *
  */

RadialPotentialBase::value_type
ZOverRFunctor::evaluate(const HFAtomicOrbitals& psi,
                        RadialOrbitalSet_t& V, int norb)
{
  RadialOrbital_t integrand(psi(0));
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type t = -Z/psi.m_grid->r(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += t;
      sum += pow(psi(o,ig),2);
    }
    integrand(ig) = t*sum;
  }
  return integrate_RK2(integrand);
}

/*!
 * \fn HarmonicFunctor::HarmonicFunctor(value_type omega)
 * \param omega
 * \brief The Constructor for the Harmonic Potential
 *
 */

HarmonicFunctor::HarmonicFunctor(value_type omega):
  Omega(omega) { }

/*! RadialPotentialBase::value_type
  PseudoPotential::evaluate(const HFAtomicOrbitals& psi,
  RadialOrbitalSet_t& V, int norb)
  * \param psi the wavefunction
  * \param V the potential
  * \param norb the number of orbitals
  * \return The sum of the Harmonic Potential matrix
  elements (energy): \f[
  \langle i|V(r)|i \rangle = \sum_{k=0}^{N_{orb}}
  \int_0^{\infty} dr\psi_k^*(r)V(r)\psi_k(r) \f]
  * \brief Calculates and assigns the values of the Harmonic
  Potential for each orbital.  The Harmonic Potential: \f[
  V_{Harmonic}(r) = \frac{1}{2}\omega^2 r^2 \f]
  *
  */

RadialPotentialBase::value_type
HarmonicFunctor::evaluate(const HFAtomicOrbitals& psi,
                          RadialOrbitalSet_t& V, int norb)
{
  RadialOrbital_t integrand(psi(0));
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type v = 0.5*Omega*Omega*psi.m_grid->r(ig)*psi.m_grid->r(ig);
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += v;
      sum += pow(psi(o,ig),2);
    }
    integrand(ig) = v*sum;
  }
  return integrate_RK2(integrand);
}

/*!
 * \fn PseudoPotential::PseudoPotential(VarRegistry<value_type>& vreg,
 value_type zeff, value_type r_c,
 value_type lambda)
 * \param vreg a registry for the optimizable variables
 \f$r_c\$f and \f$\lambda\f$
 * \param zeff \f$Z_{Eff}\f$ the effective charge of the core
 * \param r_c \f$r_c\$f the core-radius
 * \param lambda \f$\lambda\f$ the decay parameter
 * \brief The Constructor for the Joanappoulos-Starkloff
 pseudopotential
 *
 */

PseudoPotential::PseudoPotential(VarRegistry<value_type>& vreg,
                                 value_type zeff, value_type r_c,
                                 value_type lambda):
  Zeff(zeff), rc(r_c), SJ_lambda(lambda)
{
  vreg.add("SJ_lambda",&SJ_lambda);
  vreg.add("r_core",&rc);
}

PseudoPotential::PseudoPotential(value_type zeff, value_type r_c,
                                 value_type lambda):
  Zeff(zeff), rc(r_c), SJ_lambda(lambda) { }

/*! RadialPotentialBase::value_type
  PseudoPotential::evaluate(const HFAtomicOrbitals& psi,
  RadialOrbitalSet_t& V, int norb)
  * \param psi the wavefunction
  * \param V the potential
  * \param norb the number of orbitals
  * \return The sum of the PseudoPotential matrix elements: \f[
  \langle i|V_{PP}(r)|i \rangle = \sum_{k=0}^{N_{orb}}
  \int_0^{\infty} dr\psi_k^*(r)V_{PP}(r)\psi_k(r) \f]
  * \brief Calculates and assigns the values of the pseudopotential
  for each orbital.  The Starkloff-Joannapoulos pseudopotential: \f[
  V_{PP}(r) = -\frac{Z_{Eff}}{r}\frac{1-e^{-\lambda r}}
  {1-e^{-\lambda (r-r_c)}} \f]
  *
  */

RadialPotentialBase::value_type
PseudoPotential::evaluate(const HFAtomicOrbitals& psi,
                          RadialOrbitalSet_t& V, int norb)
{
  RadialOrbital_t integrand(psi(0));
  for(int ig=0; ig < psi.m_grid->size(); ig++)
  {
    value_type r = psi.m_grid->r(ig);
    value_type SJ_num = 1.0-exp(-SJ_lambda*r);
    value_type SJ_den = 1.0+exp(-SJ_lambda*(r-rc));
    value_type v = -Zeff/r*SJ_num/SJ_den;
    value_type sum = 0.0;
    for(int o=0; o < norb; o++)
    {
      V[o](ig) += v;
      sum += pow(psi(o,ig),2);
    }
    integrand(ig) = v*sum;
  }
  return integrate_RK2(integrand);
}

/*!
 * \fn HartreePotential::HartreePotential(Clebsch_Gordan* cg)
 * \param cg The Clebsch-Gordan matrix elements
 * \brief The Constructor for the HartreePotential
 *
 */

HartreePotential::HartreePotential(Clebsch_Gordan* cg):
  CG_coeff(cg) { }

/*!
 * \fn RadialPotentialBase::value_type
 HartreePotential::evaluate(const HFAtomicOrbitals& psi,
 RadialOrbitalSet_t& V, int norb)
 * \param psi The wavefuntion
 * \param V The potential
 * \param norb The number of orbitals
 * \return The sum of the Hartree matrix elements: \f[
 \langle ij|\frac{1}{r_{12}}|ij \rangle = \sum_{k=0}^{max(2l_i,2l_j)}
 c^k(l_i,m_i;l_i,m_i)c^k(l_j,m_j;l_j,m_j)R^k(i,j;i,j) \f]
 *
 * \brief Calculates and assigns the values for the Hartree potential
 for each orbital.  The Hartree potential: \f[ V_{Hartree}(r) =
 \sum_j\sum_{k=0}^{max(2l_i,2l_j)} (-1)^{m_i+m_j}
 \frac{(2l_i+1)(2l_j+1)}{(2k+1)^2}c_g(l_i,l_j,k,0,0)
 c_g(l_j,l_j,k,0,0)c_g(l_i,l_i,k,-m_i,m_i)
 c_g(l_i,l_i,k,-m_j,m_j)\frac{ Y_k(n_jl_j,n_jl_j/r)}{r} \f]
 *
 */

RadialPotentialBase::value_type
HartreePotential::evaluate(const HFAtomicOrbitals& psi,
                           RadialOrbitalSet_t& V, int norb)
{
  int kmax, k;
  int pt;
  value_type coeff, ith_orb_coeff, jth_orb_coeff, energy_coeff = 0;
  value_type Ehartree=0;
  RadialOrbital_t Ykii_r(psi(0));
  RadialOrbital_t Ykjj_r(psi(0));
  RadialOrbital_t Psisq_x_Yk(psi(0));
  int npts = psi.m_grid->size();
  for(int i=0; i < norb; i++)
  {
    int mi = psi.M[i];
    int li = psi.L[i];
    int two_li_plus_one = 2*li + 1;
    for(int j=i; j < norb; j++)
    {
      int mj = psi.M[j];
      int lj = psi.L[j];
      int two_lj_plus_one = 2*lj + 1;
      int kmax = (li > lj) ? 2*lj : 2*li;
      for(int k=kmax; k >= 0; k -= 2)
      {
        int two_k_plus_one = 2*k+1;
        int lmax = CG_coeff->Lmax;
        coeff = static_cast<value_type>(two_li_plus_one*two_lj_plus_one)/
                static_cast<value_type>(two_k_plus_one)/static_cast<value_type>(two_k_plus_one)
                * CG_coeff->cg(li,li,k,0+lmax,0+lmax) * CG_coeff->cg(lj,lj,k,0+lmax,0+lmax)
                * CG_coeff->cg(li,li,k,mi+lmax,-mi+lmax) * CG_coeff->cg(lj,lj,k,mj+lmax,-mj+lmax)
                * pow(-1.0, mi+mj);
        if(i == j)
          coeff /= 2.0;
        ith_orb_coeff = psi.Occ[i] * coeff;
        jth_orb_coeff = psi.Occ[j] * coeff;
        energy_coeff = psi.Occ[j] * psi.Occ[i] * coeff;
        Ykofr(Ykii_r, psi(i), psi(i), k);
        Ykofr(Ykjj_r, psi(j), psi(j), k);
        for(int gp=0; gp<npts; gp++)
        {
          V[i](gp) += jth_orb_coeff*Ykjj_r(gp);
          V[j](gp) += ith_orb_coeff*Ykii_r(gp);
        }
        Ehartree += Phisq_x_Yk(Ykjj_r, psi(i), psi(i), 0.5*energy_coeff);
        Ehartree += Phisq_x_Yk(Ykii_r, psi(j), psi(j), 0.5*energy_coeff);
      }
    }
  }
  return Ehartree;
}

/*!
 * \fn ExchangePotential::ExchangePotential(Clebsch_Gordan* cg)
 * \param cg The Clebsch-Gordan matrix elements
 * \brief The Constructor for the ExchangePotential
 *
 */


ExchangePotential::ExchangePotential(Clebsch_Gordan* cg):
  CG_coeff(cg) { }

/*!
 * \fn RadialPotentialBase::value_type
 ExchangePotential::evaluate(const HFAtomicOrbitals& psi,
 RadialOrbitalSet_t& V, int norb)
 * \param psi The wavefuntion
 * \param V The potential
 * \param norb The number of orbitals
 * \return The sum of the Exchange matrix elements: \f[ \langle
 ij|\frac{1}{r_{12}}|ji \rangle =
 \sum_{k=|l_i-l_j|}^{l_i+l_j}c^k(l_i,m_i;l_j,m_j)R^k(ij;ji) \f]
 * \brief Calculates and assigns the values for the exchange potential
 for each orbital.  The exchange potential: \f[ V_{exchange}(r) =
 \sum_j\delta(s_i,s_j)\sum_{k=|l_i-l_j|}^{l_i+l_j}
 \frac{(2l_i+1)(2l_j+1)}{(2k+1)^2}c_g(l_i,l_j,k,0,0)^2
 c_g(l_i,l_j,k,-m_i,m_j)^2\frac{ Y_k(n_il_i,n_jl_j/r)}{r} \f] is
 non-local.
*/

RadialPotentialBase::value_type
ExchangePotential::evaluate(const HFAtomicOrbitals& psi,
                            RadialOrbitalSet_t& V, int norb)
{
  value_type ith_orb_coeff, jth_orb_coeff, coeff;
  value_type energy_coeff=0;
  value_type Eexchange=0;
  ///    zero_all_orbitals(); V is reset before entering
  RadialOrbital_t Ykij_r(psi(0));
  // Loop over all pairs of electrons once
  for(int i=0; i < norb; i++)
  {
    int si = psi.S[i];
    int mi = psi.M[i];
    int li = psi.L[i];
    int two_li_plus_one = 2*li + 1;
    for(int j=i; j < norb; j++)
    {
      int sj = psi.S[j];
      int mj = psi.M[j];
      int lj = psi.L[j];
      int two_lj_plus_one = 2*lj + 1;
      int kmax = li + lj;
      int kmin = std::abs(li - lj);
      if( si == sj )
      {
        for(int k=kmax; k >= kmin; k-=2)
        {
          int two_k_plus_one = 2*k + 1;
          int lmax = CG_coeff->Lmax;
          coeff = static_cast<value_type>(two_li_plus_one * two_lj_plus_one) /
                  static_cast<value_type>(two_k_plus_one*two_k_plus_one)
                  * CG_coeff->cg(li,lj,k,0+lmax,0+lmax) * CG_coeff->cg(li,lj,k,-mi+lmax,mj+lmax)
                  * CG_coeff->cg(li,lj,k,0+lmax,0+lmax) * CG_coeff->cg(li,lj,k,-mi+lmax,mj+lmax);
          if(i == j)
            coeff /= 2.0;
          ith_orb_coeff = psi.Occ[i] * coeff;
          jth_orb_coeff = psi.Occ[j] * coeff;
          energy_coeff = psi.Occ[j] * psi.Occ[i] * coeff;
          Ykofr(Ykij_r, psi(i), psi(j), k); /// Ykofr_phi1_phi2
          Make_Loc_Pot(V[i], Ykij_r, psi(i), psi(j),jth_orb_coeff);
          Make_Loc_Pot(V[j], Ykij_r, psi(j), psi(i),ith_orb_coeff);
          Eexchange -= Phisq_x_Yk(Ykij_r, psi(i), psi(j), energy_coeff);
        }
      }
    }
  }
  return Eexchange;
}


