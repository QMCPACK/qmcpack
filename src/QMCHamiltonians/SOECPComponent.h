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


#ifndef QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{
/** Contains a set of radial grid potentials around a center.
*/
class SOECPComponent : public QMCTraits
{
private:
  typedef OneDimGridBase<RealType> GridType;
  typedef OneDimCubicSpline<RealType> RadialPotentialType;

  ///Non Local part: angular momentum, potential and grid
  int lmax;
  ///the number of non-local channels
  int nchannel;
  ///Maximum cutoff the non-local pseudopotential
  RealType Rmax;
  ///Angular momentum map
  aligned_vector<int> angpp_m;
  ///Non-Local part of the pseudo-potential
  std::vector<RadialPotentialType*> sopp_m;

  ParticleSet::ParticleGradient_t dG;
  ParticleSet::ParticleLaplacian_t dL;

public:
  SOECPComponent(){};

  ///destructor
  ~SOECPComponent(){};

  SOECPComponent* makeClone(const ParticleSet& qp){};

  ///add a new Spin-Orbit component
  void add(int l, RadialPotentialType* pp);
  RealType test_splined_pot(int l, RealType r);
  /** @brief Evaluate the spin orbit pp contribution 
   * to total energy from ion "iat" and electron "iel".
   *
   * @param W electron particle set.
   * @param iat index of ion.
   * @param Psi trial wave function object
   * @param iel index of electron
   * @param r the distance between ion iat and electron iel.
   * @param dr displacement from ion iat to electron iel.
   *
   * @return RealType Contribution to $\frac{V\Psi_T}{\Psi_T}$ from ion iat and electron iel.
   */
  RealType evaluateOne(ParticleSet& W,
                       int iat,
                       TrialWaveFunction& Psi,
                       int iel,
                       RealType r,
                       const PosType& dr){return 0.0;};

  // This function needs to be updated to SoA. myTableIndex is introduced temporarily.
  RealType evaluateValueAndDerivatives(ParticleSet& P,
                                       int iat,
                                       TrialWaveFunction& psi,
                                       const opt_variables_type& optvars,
                                       const std::vector<RealType>& dlogpsi,
                                       std::vector<RealType>& dhpsioverpsi,
                                       const int myTableIndex){return 0.0;};

  void print(std::ostream& os){};

  //void initVirtualParticle(const ParticleSet& qp){};

  inline void setRmax(int rmax) { Rmax = rmax; }
  inline RealType getRmax() const { return Rmax; }
  inline void setLmax(int Lmax) { lmax = Lmax; }
  inline int getLmax() const { return lmax; }

  friend class ECPComponentBuilder;
}; 

} // namespace qmcplusplus
#endif
