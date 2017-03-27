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
    
    


#ifndef OHMMS_RADIALPOTENTIALSET_H
#define OHMMS_RADIALPOTENTIALSET_H
#include "AtomicHF/RadialPotential.h"
#include "Numerics/OneDimIntegration.h"
namespace ohmmshf
{

struct RadialPotentialSet
{

  typedef RadialPotentialBase::value_type      value_type;
  typedef RadialPotentialBase::RadialOrbital_t RadialOrbital_t;
  typedef RadialPotentialBase::RadialGrid_t    RadialGrid_t;

  ///the grid potentials for each orbital
  RadialPotentialBase::RadialOrbitalSet_t V;
  ///the grid potentials for each orbital
  RadialPotentialBase::RadialOrbitalSet_t Vsave;
  ///the Self Consistent Fields
  std::vector<RadialPotentialBase*> SCF;

  ///constructor
  RadialPotentialSet() { }

  ///destructor
  ~RadialPotentialSet()
  {
    for(int i=0; i<SCF.size(); i++)
      delete SCF[i];
  }

  ///add a new Self Consistent Field to the list
  void add(RadialPotentialBase* a)
  {
    SCF.push_back(a);
  }

  ///return the number of Self Consistent Fields
  inline int size() const
  {
    return SCF.size();
  }

  ///create a grid potential for each orbital in the wavefunction
  void initialize(HFAtomicOrbitals& psi)
  {
    for(int i=0; i < psi.size(); i++)
      V.push_back(RadialOrbital_t(psi.m_grid));
    for(int i=0; i<psi.size(); i++)
      Vsave.push_back(RadialOrbital_t(psi.m_grid));
  }

  ///reset the grid potentials to zero
  void reset(HFAtomicOrbitals& psi)
  {
    for(int i=0; i<V.size(); i++)
    {
      V[i].m_Y = 0.0;
      Vsave[i].m_Y = 0.0;
    }
  }

  /*! \fn inline value_type evaluate(HFAtomicOrbitals& psi,
    std::vector<value_type>& Energy,
    int norb)
    * \param psi the wavefunction
    * \param Energy vector to store the sum of each SCF of
    the SCF matrix elements
    * \param norb number of orbitals
    * \brief Loop over all the SCFs to generate a new
    potential for each orbital.
  */

  inline value_type evaluate(HFAtomicOrbitals& psi,
                             std::vector<value_type>& Energy,
                             int norb)
  {
    value_type sum = 0.0;
    for(int i=0; i < V.size(); i++)
      V[i].m_Y = 0.0;
    for(int ip=0; ip<SCF.size(); ip++)
      sum += (Energy[ip] = SCF[ip]->evaluate(psi,V,norb));
    return sum;
  }

  /*! \fn inline void mix(value_type x)
   * \param x the mixing ratio
   * \brief Mix the old potential with the new
   potential: \f[ V_{New}(r) = (1-x)V_{New}(r) +
   x V_{Old}(r). \f]  This is for convergence purposes, so
   that the convergence is "smoother" and decreases the
   possibility of having dramatic changes in the potential.
   */

  inline void mix(value_type x)
  {
    value_type y = 1-x;
    for(int ob=0; ob < V.size(); ob++)
    {
      V[ob].m_Y = y*V[ob].m_Y + x*Vsave[ob].m_Y;
      Vsave[ob].m_Y = V[ob].m_Y;
    }
  }

  /*! \fn inline value_type calcKE(HFAtomicOrbitals& psi,
    value_type eigsum,
    int norb)
    * \param psi the wavefunction
    * \param eigsum the sum of the eigenvalues
    * \param norb the number of orbitals
    * \return the total kinetic energy
    * \brief Calculates the total kinetic energy:
    \f[ KE = \sum_i^{Norb} \left\{ \epsilon_i - \int_0^{\infty}
    dr \psi_i V(r) \psi_i \right\}. \f]  This is true because
    \f[ \langle \psi_i | -\frac{1}{2} \nabla^2 | \psi_i \rangle
    = -\frac{1}{2r^2}\frac{d}{dr}\left( r^2 \frac{d}{dr} -
    \frac{l(l+1)}{2r^2} = (\epsilon - V(r)) \f]
  */

  inline value_type calcKE(HFAtomicOrbitals& psi,
                           value_type eigsum,
                           int norb)
  {
    value_type sum=0.0;
    RadialOrbital_t integrand(psi(0));
    for(int ob=0; ob < norb; ob++)
    {
      for(int i=0; i < integrand.size(); i++)
      {
        integrand(i) = V[ob](i)*psi(ob,i)*psi(ob,i);
      }
      sum+=integrate_RK2(integrand);
    }
    return eigsum-sum;
  }

  /*! \fn void applyRestriction(HFAtomicOrbitals& psi)
   * \param psi the wavefunction
   * \brief Restrict the potential.  Normally each orbital
   \f$ \psi_i \f$ in the wavefunction has its own unique
   potential, but we want to restrict the potential to be
   the same for orbitals with the same quantum numbers,
   such as \f$ (n,l) \f$.  What this function does is
   assign the average potential to all the orbitals that
   are restricted to be the same.
  */

  void applyRestriction(HFAtomicOrbitals& psi)
  {
    static std::vector<value_type> sum;
    if(sum.empty())
    {
      sum.resize(psi.m_grid->size());
      for(int ig=0; ig < psi.m_grid->size(); ig++)
        sum[ig] = 0.0;
    }
    ///index of starting orbital index
    int o_start = 0;
    ///index of ending orbital index
    int o_end = 0;
    ///orbital index
    int orb = 0;
    while (orb < psi.size())
    {
      ///loop over unique orbitals
      for(int uorb=0; uorb < psi.NumUniqueOrb; uorb++)
      {
        ///for each unique orbital, loop over all
        ///identical orbitals
        for(int i=0; i < psi.IDcount[uorb]; i++)
        {
          ///add all the orbitals together for averaging
          for(int ig=0; ig < psi.m_grid->size(); ig++)
          {
            sum[ig] += V[orb](ig);
          }
          ///increment the orbital index
          orb++;
        }
        int o_end = o_start+psi.IDcount[uorb];
        ///assign the average back to the orbitals
        for(int o = o_start; o < o_end; o++)
        {
          for(int ig=0; ig < psi.m_grid->size(); ig++)
          {
            V[o](ig) = sum[ig]/psi.IDcount[uorb];
          }
        }
        o_start = o_end;
        ///reset the sum for the next average
        for(int ig=0; ig < psi.m_grid->size(); ig++)
          sum[ig] = 0.0;
      }
    }
  }

};

}
#endif
