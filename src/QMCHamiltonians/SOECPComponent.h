//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCHamiltonians/RandomRotationMatrix.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include <ResourceCollection.h>
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{

namespace testing
{
class TestSOECPotential;
}
/** class SOECPComponent
 **  brief Computes the nonlocal spin-orbit interaction \f$\Delta V_SO(r) |ljm_j><ljm_j|\f$.
 **  details This computes the nonlocal spin-orbit interaction between a single ion species and 
 **           a given electron.  
 **           Currently, this class does nothing other than generate and store \f$\Delta V_SO(r)\f$
 **           for different orbital angular momenta.  Implementation coming soon!  
 **/

class SOECPComponent : public QMCTraits
{
private:
  using SpherGridType       = std::vector<PosType>;
  using GridType            = OneDimGridBase<RealType>;
  using RadialPotentialType = OneDimCubicSpline<RealType>;

  ///Non Local part: angular momentum, potential and grid
  int lmax_;
  ///the number of non-local channels
  int nchannel_;
  int nknot_;
  int sknot_;
  int total_knots_; //spin + spatial knots
  ///Maximum cutoff the non-local pseudopotential
  RealType rmax_;
  ///Angular momentum map
  aligned_vector<int> angpp_m_;
  ///Non-Local part of the pseudo-potential
  std::vector<RadialPotentialType*> sopp_m_;

  ComplexType sMatrixElements(RealType s1, RealType s2, int dim);
  ComplexType lmMatrixElements(int l, int m1, int m2, int dim);
  int kroneckerDelta(int x, int y);

  std::vector<PosType> deltaV_;
  std::vector<RealType> deltaS_;
  SpherGridType sgridxyz_m_;
  SpherGridType rrotsgrid_m_;
  std::vector<ValueType> psiratio_;
  std::vector<ValueType> vrad_;
  std::vector<RealType> sgridweight_m_;
  //total spin and quadrature weights
  std::vector<RealType> spin_quad_weights_;
  //work array
  std::vector<ValueType> wvec_;
  //scratch spaces used by evaluateValueAndDerivative
  Matrix<ValueType> dratio_;
  Vector<ValueType> dlogpsi_vp_;
  VirtualParticleSet* vp_;

  //This builds the full quadrature grid for the Simpsons rule used for spin integrals as well as
  //the spatial quadrature. In this function, it specifies the deltaS_ and deltaV_ for all the quadrature points and sets the interal weights
  //in spin_quad_weights
  //If there are s0,s1,...sN spin integral points and q0,q1,...qM spatial quadrature points, the order is
  // s0q0, s0q1, ..., s0qM, s1q0, ..., sNq0, ..., sNqM for each of the deltaS_, deltaV_, and spin_quad_weights_
  void buildTotalQuadrature(const RealType r, const PosType& dr, const RealType sold);

public:
  SOECPComponent();
  ~SOECPComponent();

  SOECPComponent* makeClone(const ParticleSet& qp);

  ///add a new Spin-Orbit component
  void add(int l, RadialPotentialType* pp);

  void resize_warrays(int n, int m, int s);

  void rotateQuadratureGrid(const TensorType& rmat);

  ///API for accessing the value of an SO radial potential at distance r.  For unit and other testing.
  friend RealType getSplinedSOPot(SOECPComponent* so_pp, int l, double r);
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
  RealType evaluateOne(ParticleSet& W, int iat, TrialWaveFunction& Psi, int iel, RealType r, const PosType& dr);

  RealType calculateProjector(RealType r, const PosType& dr, RealType sold);


  static void mw_evaluateOne(const RefVectorWithLeader<SOECPComponent>& soecp_component_list,
                             const RefVectorWithLeader<ParticleSet>& p_list,
                             const RefVectorWithLeader<TrialWaveFunction>& psi_list,
                             const RefVector<const NLPPJob<RealType>>& joblist,
                             std::vector<RealType>& pairpots,
                             ResourceCollection& collection);

  RealType evaluateValueAndDerivatives(ParticleSet& P,
                                       int iat,
                                       TrialWaveFunction& psi,
                                       int iel,
                                       RealType r,
                                       const PosType& dr,
                                       const opt_variables_type& optvars,
                                       const Vector<ValueType>& dlogpsi,
                                       Vector<ValueType>& dhpsioverpsi);

  void print(std::ostream& os);

  void initVirtualParticle(const ParticleSet& qp);
  void deleteVirtualParticle();

  inline void setRmax(RealType rmax) { rmax_ = rmax; }
  inline RealType getRmax() const { return rmax_; }
  inline void setLmax(int lmax) { lmax_ = lmax; }
  inline int getLmax() const { return lmax_; }
  inline int getNknot() const { return nknot_; }
  inline int getSknot() const { return sknot_; }

  const VirtualParticleSet* getVP() const { return vp_; };

  friend struct ECPComponentBuilder;
  friend void copyGridUnrotatedForTest(SOECPComponent& nlpp);

  friend class testing::TestSOECPotential;
};

} // namespace qmcplusplus
#endif
