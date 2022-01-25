//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratory
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#define QMCPLUSPLUS_SO_ECPOTENTIAL_COMPONENT_H
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimLinearSpline.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus
{
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
  int lmax;
  ///the number of non-local channels
  int nchannel;
  int nknot;
  int sknot;
  ///Maximum cutoff the non-local pseudopotential
  RealType Rmax;
  ///Angular momentum map
  aligned_vector<int> angpp_m;
  ///Non-Local part of the pseudo-potential
  std::vector<RadialPotentialType*> sopp_m;

  ComplexType sMatrixElements(RealType s1, RealType s2, int dim);
  ComplexType lmMatrixElements(int l, int m1, int m2, int dim);
  int kroneckerDelta(int x, int y);

  ComplexType getAngularIntegral(RealType sold,
                                 RealType snew,
                                 ParticleSet& W,
                                 TrialWaveFunction& Psi,
                                 int iel,
                                 RealType r,
                                 const PosType& dr);

  std::vector<PosType> deltaV;
  SpherGridType sgridxyz_m;
  SpherGridType rrotsgrid_m;
  std::vector<ValueType> psiratio;
  std::vector<ValueType> vrad;
  std::vector<RealType> sgridweight_m;


public:
  SOECPComponent();
  ~SOECPComponent();

  SOECPComponent* makeClone(const ParticleSet& qp);

  ///add a new Spin-Orbit component
  void add(int l, RadialPotentialType* pp);

  void resize_warrays(int n, int m, int s);

  void randomize_grid(RandomGenerator& myRNG);
  template<typename T>
  void randomize_grid(std::vector<T>& sphere, RandomGenerator& myRNG);

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

  // This function needs to be updated to SoA. myTableIndex is introduced temporarily.
  inline RealType evaluateValueAndDerivatives(ParticleSet& P,
                                              int iat,
                                              TrialWaveFunction& psi,
                                              const opt_variables_type& optvars,
                                              const std::vector<RealType>& dlogpsi,
                                              std::vector<RealType>& dhpsioverpsi,
                                              const int myTableIndex)
  {
    APP_ABORT("evaluateValueAndDerivatives not implemented yet\n");
    return 0.0;
  };

  void print(std::ostream& os);

  //void initVirtualParticle(const ParticleSet& qp){};

  inline void setRmax(int rmax) { Rmax = rmax; }
  inline RealType getRmax() const { return Rmax; }
  inline void setLmax(int Lmax) { lmax = Lmax; }
  inline int getLmax() const { return lmax; }
  inline int getNknot() const { return nknot; }
  inline int getSknot() const { return sknot; }

  friend struct ECPComponentBuilder;
  friend void copyGridUnrotatedForTest(SOECPComponent& nlpp);
};

} // namespace qmcplusplus
#endif
