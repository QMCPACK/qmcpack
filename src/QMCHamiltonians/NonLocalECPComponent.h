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


#ifndef QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_NONLOCAL_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/OperatorBase.h"
#include <ResourceCollection.h>
#include <TrialWaveFunction.h>
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"
#include "NLPPJob.h"

namespace qmcplusplus
{
/** Contains a set of radial grid potentials around a center.
*/
class NonLocalECPComponent : public QMCTraits
{
private:
  using SpherGridType       = std::vector<PosType>;
  using GridType            = OneDimGridBase<RealType>;
  using RadialPotentialType = OneDimCubicSpline<RealType>;

  ///Non Local part: angular momentum, potential and grid
  int lmax;
  ///the number of non-local channels
  int nchannel;
  ///the number of nknot
  int nknot;
  ///Maximum cutoff the non-local pseudopotential
  RealType Rmax;
  ///Angular momentum map
  aligned_vector<int> angpp_m;
  ///Weight of the angular momentum
  aligned_vector<RealType> wgt_angpp_m;
  /// Lfactor1[l]=(2*l+1)/(l+1)
  RealType Lfactor1[8];
  /// Lfactor1[l]=(l)/(l+1)
  RealType Lfactor2[8];
  ///Non-Local part of the pseudo-potential
  std::vector<RadialPotentialType*> nlpp_m;
  ///fixed Spherical Grid for species
  SpherGridType sgridxyz_m;
  ///randomized spherical grid
  SpherGridType rrotsgrid_m;
  ///weight of the spherical grid
  std::vector<RealType> sgridweight_m;
  ///Working arrays
  std::vector<ValueType> wvec, Amat, dAmat;

  //Position delta for virtual moves.
  std::vector<PosType> deltaV;
  //Array for P_l[cos(theta)].
  std::vector<RealType> lpol;
  //Array for P'_l[cos(theta)]
  std::vector<RealType> dlpol;
  //Array for v_l(r).
  std::vector<ValueType> vrad;
  //Array for (2l+1)*v'_l(r)/r.
  std::vector<RealType> dvrad;
  //$\Psi(...q...)/\Psi(...r...)$ for all quadrature points q.
  std::vector<ValueType> psiratio;
  //$\nabla \Psi(...q...)/\Psi(...r...)$ for all quadrature points q.
  //  $\nabla$ is w.r.t. the electron coordinates involved in the quadrature.
  std::vector<PosType> gradpsiratio;
  //This stores gradient of v(r):
  std::vector<PosType> vgrad;
  //This stores the gradient of the cos(theta) term in force expression.
  std::vector<PosType> cosgrad;
  //This stores grad psi/psi - dot(u,grad psi)
  std::vector<PosType> wfngrad;
  //This stores potential contribution per knot:
  std::vector<RealType> knot_pots;

  /// scratch spaces used by evaluateValueAndDerivatives
  Matrix<ValueType> dratio;
  std::vector<ValueType> dlogpsi_vp;

  // For Pulay correction to the force
  std::vector<RealType> WarpNorm;
  ParticleSet::ParticleGradient dG;
  ParticleSet::ParticleLaplacian dL;
  /// First index is knot, second is electron
  Matrix<PosType> Gnew;
  ///The gradient of the wave function w.r.t. the ion position
  ParticleSet::ParticleGradient Gion;

  ///virtual particle set: delayed initialization
  VirtualParticleSet* VP;

  /// build QP position deltas from the reference electron using internally stored random grid points
  void buildQuadraturePointDeltaPositions(RealType r, const PosType& dr, std::vector<PosType>& deltaV) const;

  /** finalize the calculation of $\frac{V\Psi_T}{\Psi_T}$
   */
  RealType calculateProjector(RealType r, const PosType& dr);

public:
  NonLocalECPComponent();

  ///destructor
  ~NonLocalECPComponent();

  NonLocalECPComponent* makeClone(const ParticleSet& qp);

  ///add a new Non Local component
  void add(int l, RadialPotentialType* pp);

  ///add knots to the spherical grid
  void addknot(const PosType& xyz, RealType weight)
  {
    sgridxyz_m.push_back(xyz);
    sgridweight_m.push_back(weight);
  }

  void resize_warrays(int n, int m, int l);

  void randomize_grid(RandomGenerator& myRNG);
  template<typename T>
  void randomize_grid(std::vector<T>& sphere, RandomGenerator& myRNG);

  /** contribute local non-local move data
   * @param iel reference electron id.
   * @param Txy nonlocal move data.
   */
  void contributeTxy(int iel, std::vector<NonLocalData>& Txy) const;

  /** @brief Evaluate the nonlocal pp contribution via randomized quadrature grid
   * to total energy from ion "iat" and electron "iel".
   *
   * @param W electron particle set.
   * @param iat index of ion.
   * @param Psi trial wave function object
   * @param iel index of electron
   * @param r the distance between ion iat and electron iel.
   * @param dr displacement from ion iat to electron iel.
   * @param use_DLA if ture, use determinant localization approximation (DLA).
   *
   * @return RealType Contribution to $\frac{V\Psi_T}{\Psi_T}$ from ion iat and electron iel.
   */
  RealType evaluateOne(ParticleSet& W,
                       int iat,
                       TrialWaveFunction& Psi,
                       int iel,
                       RealType r,
                       const PosType& dr,
                       bool use_DLA);

  /** @brief Evaluate the nonlocal pp contribution via randomized quadrature grid
   * to total energy from ion "iat" and electron "iel" for a batch of walkers.
   *
   * @param ecp_component_list a list of ECP components
   * @param p_list a list of electron particle set.
   * @param psi_list a list of trial wave function object
   * @param joblist a list of ion-electron pairs
   * @param pairpots a list of contribution to $\frac{V\Psi_T}{\Psi_T}$ from ion iat and electron iel.
   * @param use_DLA if ture, use determinant localization approximation (DLA).
   *
   * Note: ecp_component_list allows including different NLPP component for different walkers.
   * electrons in iel_list must be of the same group (spin)
   */
  static void mw_evaluateOne(const RefVectorWithLeader<NonLocalECPComponent>& ecp_component_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             const RefVectorWithLeader<TrialWaveFunction>& psi_list,
                             const RefVector<const NLPPJob<RealType>>& joblist,
                             std::vector<RealType>& pairpots,
                             ResourceCollection& collection,
                             bool use_DLA);

  /** @brief Evaluate the nonlocal pp contribution via randomized quadrature grid
   * to total energy from ion "iat" and electron "iel".
   *
   * @param W electron particle set.
   * @param iat index of ion.
   * @param Psi trial wave function object
   * @param iel index of electron
   * @param r the distance between ion iat and electron iel.
   * @param dr displacement from ion iat to electron iel.
   * @param force_iat 3d vector for Hellman-Feynman contribution.  This gets modified.
   *
   * @return RealType Contribution to $\frac{V\Psi_T}{\Psi_T}$ from ion iat and electron iel.
   */
  RealType evaluateOneWithForces(ParticleSet& W,
                                 int iat,
                                 TrialWaveFunction& Psi,
                                 int iel,
                                 RealType r,
                                 const PosType& dr,
                                 PosType& force_iat);

  /** @brief Evaluate the nonlocal pp energy, Hellman-Feynman force, and "Pulay" force contribution
   * via randomized quadrature grid from ion "iat" and electron "iel".
   *
   * @param W electron particle set.
   * @param ions ion particle set.
   * @param iat index of ion.
   * @param Psi trial wave function object
   * @param iel index of electron
   * @param r the distance between ion iat and electron iel.
   * @param dr displacement from ion iat to electron iel.
   * @param force_iat 3d vector for Hellman-Feynman contribution.  This gets modified.
   * @param pulay_terms Nion x 3 object, holding a contribution for each ionic gradient from \Psi_T.
   *
   * @return RealType Contribution to $\frac{V\Psi_T}{\Psi_T}$ from ion iat and electron iel.
   */
  RealType evaluateOneWithForces(ParticleSet& W,
                                 ParticleSet& ions,
                                 int iat,
                                 TrialWaveFunction& Psi,
                                 int iel,
                                 RealType r,
                                 const PosType& dr,
                                 PosType& force_iat,
                                 ParticleSet::ParticlePos& pulay_terms);

  // This function needs to be updated to SoA. myTableIndex is introduced temporarily.
  RealType evaluateValueAndDerivatives(ParticleSet& P,
                                       int iat,
                                       TrialWaveFunction& psi,
                                       int iel,
                                       RealType r,
                                       const PosType& dr,
                                       const opt_variables_type& optvars,
                                       const std::vector<ValueType>& dlogpsi,
                                       std::vector<ValueType>& dhpsioverpsi);

  void print(std::ostream& os);

  void initVirtualParticle(const ParticleSet& qp);

  inline void setRmax(int rmax) { Rmax = rmax; }
  inline RealType getRmax() const { return Rmax; }
  inline int getNknot() const { return nknot; }
  inline void setLmax(int Lmax) { lmax = Lmax; }
  inline int getLmax() const { return lmax; }
  const VirtualParticleSet* getVP() const { return VP; };

  // copy sgridxyz_m to rrotsgrid_m without rotation. For testing only.
  friend void copyGridUnrotatedForTest(NonLocalECPComponent& nlpp);

  friend struct ECPComponentBuilder;
  // a lazy temporal solution
  friend class NonLocalECPotential_CUDA;
}; //end of RadialPotentialSet

} // namespace qmcplusplus
#endif
