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


#include "Particle/DistanceTable.h"
#include "NonLocalECPComponent.h"
#include "NLPPJob.h"
#include "NonLocalData.h"
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
NonLocalECPComponent::NonLocalECPComponent() : lmax(0), nchannel(0), nknot(0), Rmax(-1), VP(nullptr) {}

NonLocalECPComponent::~NonLocalECPComponent()
{
  for (int ip = 0; ip < nlpp_m.size(); ip++)
    delete nlpp_m[ip];
  if (VP)
    delete VP;
}

NonLocalECPComponent* NonLocalECPComponent::makeClone(const ParticleSet& qp)
{
  NonLocalECPComponent* myclone = new NonLocalECPComponent(*this);
  for (int i = 0; i < nlpp_m.size(); ++i)
    myclone->nlpp_m[i] = nlpp_m[i]->makeClone();
  if (VP)
    myclone->VP = new VirtualParticleSet(qp, nknot);
  return myclone;
}

void NonLocalECPComponent::initVirtualParticle(const ParticleSet& qp)
{
  assert(VP == 0);
  outputManager.pause();
  VP = new VirtualParticleSet(qp, nknot);
  outputManager.resume();
}

void NonLocalECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m.push_back(l);
  wgt_angpp_m.push_back(static_cast<RealType>(2 * l + 1));
  nlpp_m.push_back(pp);
}

void NonLocalECPComponent::resize_warrays(int n, int m, int l)
{
  psiratio.resize(n);
  gradpsiratio.resize(n);
  deltaV.resize(n);
  cosgrad.resize(n);
  wfngrad.resize(n);
  knot_pots.resize(n);
  vrad.resize(m);
  dvrad.resize(m);
  vgrad.resize(m);
  wvec.resize(m);
  Amat.resize(n * m);
  dAmat.resize(n * m);
  lpol.resize(l + 1, 1.0);
  dlpol.resize(l + 1, 0.0);
  rrotsgrid_m.resize(n);
  nchannel = nlpp_m.size();
  nknot    = sgridxyz_m.size();

  //Now we inititalize the quadrature grid rrotsgrid_m to the unrotated grid.
  rrotsgrid_m = sgridxyz_m;

  //This is just to check
  //for(int nl=1; nl<nlpp_m.size(); nl++) nlpp_m[nl]->setGridManager(false);
  if (lmax)
  {
    if (lmax > 7)
    {
      APP_ABORT("Increase the maximum angular momentum implemented.");
    }
    //Lfactor1.resize(lmax);
    //Lfactor2.resize(lmax);
    for (int nl = 0; nl < lmax; nl++)
    {
      Lfactor1[nl] = static_cast<RealType>(2 * nl + 1);
      Lfactor2[nl] = 1.0e0 / static_cast<RealType>(nl + 1);
    }
  }
}

void NonLocalECPComponent::print(std::ostream& os)
{
  os << "    Maximum angular mementum = " << lmax << std::endl;
  os << "    Number of non-local channels = " << nchannel << std::endl;
  for (int l = 0; l < nchannel; l++)
    os << "       l(" << l << ")=" << angpp_m[l] << std::endl;
  os << "    Cutoff radius = " << Rmax << std::endl;
  os << "    Spherical grids and weights: " << std::endl;
  for (int ik = 0; ik < nknot; ik++)
    os << "       " << sgridxyz_m[ik] << std::setw(20) << sgridweight_m[ik] << std::endl;
}

NonLocalECPComponent::RealType NonLocalECPComponent::evaluateOne(ParticleSet& W,
                                                                 int iat,
                                                                 TrialWaveFunction& psi,
                                                                 int iel,
                                                                 RealType r,
                                                                 const PosType& dr,
                                                                 bool use_DLA)
{
  buildQuadraturePointDeltaPositions(r, dr, deltaV);

  if (VP)
  {
    // Compute ratios with VP
    VP->makeMoves(iel, W.R[iel], deltaV, true, iat);
    if (use_DLA)
      psi.evaluateRatios(*VP, psiratio, TrialWaveFunction::ComputeType::FERMIONIC);
    else
      psi.evaluateRatios(*VP, psiratio);
  }
  else
  {
    // Compute ratio of wave functions
    for (int j = 0; j < nknot; j++)
    {
      W.makeMove(iel, deltaV[j], false);
      if (use_DLA)
        psiratio[j] = psi.calcRatio(W, iel, TrialWaveFunction::ComputeType::FERMIONIC);
      else
        psiratio[j] = psi.calcRatio(W, iel);
      W.rejectMove(iel);
      psi.resetPhaseDiff();
    }
  }

  return calculateProjector(r, dr);
}

NonLocalECPComponent::RealType NonLocalECPComponent::calculateProjector(RealType r, const PosType& dr)
{
  for (int j = 0; j < nknot; j++)
    psiratio[j] *= sgridweight_m[j];

  // Compute radial potential, multiplied by (2l+1) factor.
  for (int ip = 0; ip < nchannel; ip++)
    vrad[ip] = nlpp_m[ip]->splint(r) * wgt_angpp_m[ip];

  constexpr RealType czero(0);
  constexpr RealType cone(1);

  const RealType rinv = cone / r;
  RealType pairpot    = czero;
  // Compute spherical harmonics on grid
  for (int j = 0; j < nknot; j++)
  {
    RealType zz = dot(dr, rrotsgrid_m[j]) * rinv;
    // Forming the Legendre polynomials
    lpol[0]           = cone;
    RealType lpolprev = czero;
    for (int l = 0; l < lmax; l++)
    {
      //Not a big difference
      //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
      //lpol[l+1]/=(l+1);
      lpol[l + 1] = Lfactor1[l] * zz * lpol[l] - l * lpolprev;
      lpol[l + 1] *= Lfactor2[l];
      lpolprev = lpol[l];
    }

    ValueType lsum = 0.0;
    for (int l = 0; l < nchannel; l++)
      lsum += vrad[l] * lpol[angpp_m[l]];
    knot_pots[j] = std::real(lsum * psiratio[j]);
    pairpot += knot_pots[j];
  }

  return pairpot;
}

void NonLocalECPComponent::mw_evaluateOne(const RefVectorWithLeader<NonLocalECPComponent>& ecp_component_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          const RefVectorWithLeader<TrialWaveFunction>& psi_list,
                                          const RefVector<const NLPPJob<RealType>>& joblist,
                                          std::vector<RealType>& pairpots,
                                          ResourceCollection& collection,
                                          bool use_DLA)
{
  auto& ecp_component_leader = ecp_component_list.getLeader();
  if (ecp_component_leader.VP)
  {
    // Compute ratios with VP
    RefVectorWithLeader<VirtualParticleSet> vp_list(*ecp_component_leader.VP);
    RefVectorWithLeader<const VirtualParticleSet> const_vp_list(*ecp_component_leader.VP);
    RefVector<const std::vector<PosType>> deltaV_list;
    RefVector<std::vector<ValueType>> psiratios_list;
    vp_list.reserve(ecp_component_list.size());
    const_vp_list.reserve(ecp_component_list.size());
    deltaV_list.reserve(ecp_component_list.size());
    psiratios_list.reserve(ecp_component_list.size());

    for (size_t i = 0; i < ecp_component_list.size(); i++)
    {
      NonLocalECPComponent& component(ecp_component_list[i]);
      const NLPPJob<RealType>& job = joblist[i];

      component.buildQuadraturePointDeltaPositions(job.ion_elec_dist, job.ion_elec_displ, component.deltaV);

      vp_list.push_back(*component.VP);
      const_vp_list.push_back(*component.VP);
      deltaV_list.push_back(component.deltaV);
      psiratios_list.push_back(component.psiratio);
    }

    ResourceCollectionTeamLock<VirtualParticleSet> vp_res_lock(collection, vp_list);

    VirtualParticleSet::mw_makeMoves(vp_list, deltaV_list, joblist, true);

    if (use_DLA)
      TrialWaveFunction::mw_evaluateRatios(psi_list, const_vp_list, psiratios_list,
                                           TrialWaveFunction::ComputeType::FERMIONIC);
    else
      TrialWaveFunction::mw_evaluateRatios(psi_list, const_vp_list, psiratios_list);
  }
  else
  {
    // Compute ratios without VP. This is working but very slow code path.
#pragma omp parallel for
    for (size_t i = 0; i < p_list.size(); i++)
    {
      NonLocalECPComponent& component(ecp_component_list[i]);
      ParticleSet& W(p_list[i]);
      TrialWaveFunction& psi(psi_list[i]);
      const NLPPJob<RealType>& job = joblist[i];

      component.buildQuadraturePointDeltaPositions(job.ion_elec_dist, job.ion_elec_displ, component.deltaV);

      // Compute ratio of wave functions
      for (int j = 0; j < component.getNknot(); j++)
      {
        W.makeMove(job.electron_id, component.deltaV[j], false);
        if (use_DLA)
          component.psiratio[j] = psi.calcRatio(W, job.electron_id, TrialWaveFunction::ComputeType::FERMIONIC);
        else
          component.psiratio[j] = psi.calcRatio(W, job.electron_id);
        W.rejectMove(job.electron_id);
        psi.resetPhaseDiff();
      }
    }
  }

  for (size_t i = 0; i < p_list.size(); i++)
  {
    NonLocalECPComponent& component(ecp_component_list[i]);
    const NLPPJob<RealType>& job = joblist[i];
    pairpots[i]                  = component.calculateProjector(job.ion_elec_dist, job.ion_elec_displ);
  }
}

NonLocalECPComponent::RealType NonLocalECPComponent::evaluateOneWithForces(ParticleSet& W,
                                                                           int iat,
                                                                           TrialWaveFunction& psi,
                                                                           int iel,
                                                                           RealType r,
                                                                           const PosType& dr,
                                                                           PosType& force_iat)
{
  constexpr RealType czero(0);
  constexpr RealType cone(1);

  //We check that our quadrature grid is valid.  Namely, that all points lie on the unit sphere.
  //We check this by seeing if |r|^2 = 1 to machine precision.
  for (int j = 0; j < nknot; j++)
    assert(std::abs(std::sqrt(dot(rrotsgrid_m[j], rrotsgrid_m[j])) - 1) <
           100 * std::numeric_limits<RealType>::epsilon());


  for (int j = 0; j < nknot; j++)
    deltaV[j] = r * rrotsgrid_m[j] - dr;

  GradType gradtmp_(0);
  PosType realgradtmp_(0);

  //Pseudopotential derivative w.r.t. ions can be split up into 3 contributions:
  // term coming from the gradient of the radial potential
  PosType gradpotterm_(0);
  // term coming from gradient of legendre polynomial
  PosType gradlpolyterm_(0);
  // term coming from dependence of quadrature grid on ion position.
  PosType gradwfnterm_(0);

  if (VP)
  {
    APP_ABORT("NonLocalECPComponent::evaluateOneWithForces(...): Forces not implemented with virtual particle moves\n");
    // Compute ratios with VP
    VP->makeMoves(iel, W.R[iel], deltaV, true, iat);
    psi.evaluateRatios(*VP, psiratio);
  }
  else
  {
    // Compute ratio of wave functions
    for (int j = 0; j < nknot; j++)
    {
      W.makeMove(iel, deltaV[j], false);
      psiratio[j] = psi.calcRatioGrad(W, iel, gradtmp_);
      //QMCPACK spits out $\nabla\Psi(q)/\Psi(q)$.
      //Multiply times $\Psi(q)/\Psi(r)$ to get
      // $\nabla\Psi(q)/\Psi(r)
      gradtmp_ *= psiratio[j];
#if defined(QMC_COMPLEX)
      //And now we take the real part and save it.
      convertToReal(gradtmp_, gradpsiratio[j]);
#else
      //Real nonlocalpp forces seem to differ from those in the complex build.  Since
      //complex build has been validated against QE, that indicates there's a bug for the real build.
      gradpsiratio[j] = gradtmp_;
#endif
      W.rejectMove(iel);
      psi.resetPhaseDiff();
      //psi.rejectMove(iel);
    }
  }

  for (int j = 0; j < nknot; j++)
    psiratio[j] *= sgridweight_m[j];

  // This is just a temporary variable to dump d2/dr2 into for spline evaluation.
  RealType secondderiv(0);

  const RealType rinv = cone / r;

  // Compute radial potential and its derivative times (2l+1)
  for (int ip = 0; ip < nchannel; ip++)
  {
    //fun fact.  NLPComponent stores v(r) as v(r), and not as r*v(r) like in other places.
    vrad[ip]  = nlpp_m[ip]->splint(r, dvrad[ip], secondderiv) * wgt_angpp_m[ip];
    vgrad[ip] = dvrad[ip] * dr * wgt_angpp_m[ip] * rinv;
  }

  RealType pairpot = 0;
  // Compute spherical harmonics on grid
  for (int j = 0; j < nknot; j++)
  {
    RealType zz        = dot(dr, rrotsgrid_m[j]) * rinv;
    PosType uminusrvec = rrotsgrid_m[j] - zz * dr * rinv;

    cosgrad[j] = rinv * uminusrvec;

    RealType udotgradpsi = dot(gradpsiratio[j], rrotsgrid_m[j]);
    wfngrad[j]           = gradpsiratio[j] - dr * (udotgradpsi * rinv);
    wfngrad[j] *= sgridweight_m[j];

    // Forming the Legendre polynomials
    //P_0(x)=1; P'_0(x)=0.
    lpol[0]  = cone;
    dlpol[0] = czero;

    RealType lpolprev  = czero;
    RealType dlpolprev = czero;

    for (int l = 0; l < lmax; l++)
    {
      //Legendre polynomial recursion formula.
      lpol[l + 1] = Lfactor1[l] * zz * lpol[l] - l * lpolprev;
      lpol[l + 1] *= Lfactor2[l];

      //and for the derivative...
      dlpol[l + 1] = Lfactor1[l] * (zz * dlpol[l] + lpol[l]) - l * dlpolprev;
      dlpol[l + 1] *= Lfactor2[l];

      lpolprev  = lpol[l];
      dlpolprev = dlpol[l];
    }

    RealType lsum = czero;
    // Now to compute the forces:
    gradpotterm_   = 0;
    gradlpolyterm_ = 0;
    gradwfnterm_   = 0;

    for (int l = 0; l < nchannel; l++)
    {
      lsum += std::real(vrad[l]) * lpol[angpp_m[l]];
      gradpotterm_ += vgrad[l] * lpol[angpp_m[l]] * std::real(psiratio[j]);
      gradlpolyterm_ += std::real(vrad[l]) * dlpol[angpp_m[l]] * cosgrad[j] * std::real(psiratio[j]);
      gradwfnterm_ += std::real(vrad[l]) * lpol[angpp_m[l]] * wfngrad[j];
    }
    knot_pots[j] = std::real(lsum * psiratio[j]);
    pairpot += knot_pots[j];
    force_iat += gradpotterm_ + gradlpolyterm_ - gradwfnterm_;
  }

  return pairpot;
}

NonLocalECPComponent::RealType NonLocalECPComponent::evaluateOneWithForces(ParticleSet& W,
                                                                           ParticleSet& ions,
                                                                           int iat,
                                                                           TrialWaveFunction& psi,
                                                                           int iel,
                                                                           RealType r,
                                                                           const PosType& dr,
                                                                           PosType& force_iat,
                                                                           ParticleSet::ParticlePos& pulay_terms)
{
  constexpr RealType czero(0);
  constexpr RealType cone(1);

  //We check that our quadrature grid is valid.  Namely, that all points lie on the unit sphere.
  //We check this by seeing if |r|^2 = 1 to machine precision.
  for (int j = 0; j < nknot; j++)
    assert(std::abs(std::sqrt(dot(rrotsgrid_m[j], rrotsgrid_m[j])) - 1) <
           100 * std::numeric_limits<RealType>::epsilon());

  for (int j = 0; j < nknot; j++)
    deltaV[j] = r * rrotsgrid_m[j] - dr;

  GradType gradtmp_(0);
  PosType realgradtmp_(0);

  //Pseudopotential derivative w.r.t. ions can be split up into 3 contributions:
  // term coming from the gradient of the radial potential
  PosType gradpotterm_(0);
  // term coming from gradient of legendre polynomial
  PosType gradlpolyterm_(0);
  // term coming from dependence of quadrature grid on ion position.
  PosType gradwfnterm_(0);

  //Now for the Pulay specific stuff...
  // $\nabla_I \Psi(...r...)/\Psi(...r...)$
  ParticleSet::ParticlePos pulay_ref;
  ParticleSet::ParticlePos pulaytmp_;
  // $\nabla_I \Psi(...q...)/\Psi(...r...)$ for each quadrature point.
  std::vector<ParticleSet::ParticlePos> pulay_quad(nknot);

  //A working array for pulay stuff.
  GradType iongradtmp_(0);
  //resize everything.
  pulay_ref.resize(ions.getTotalNum());
  pulaytmp_.resize(ions.getTotalNum());
  for (size_t j = 0; j < nknot; j++)
    pulay_quad[j].resize(ions.getTotalNum());

  if (VP)
  {
    APP_ABORT("NonLocalECPComponent::evaluateOneWithForces(...): Forces not implemented with virtual particle moves\n");
    // Compute ratios with VP
    VP->makeMoves(iel, W.R[iel], deltaV, true, iat);
    psi.evaluateRatios(*VP, psiratio);
  }
  else
  {
    // Compute ratio of wave functions
    for (int j = 0; j < nknot; j++)
    {
      W.makeMove(iel, deltaV[j], false);
      psiratio[j] = psi.calcRatioGrad(W, iel, gradtmp_);
      //QMCPACK spits out $\nabla\Psi(q)/\Psi(q)$.
      //Multiply times $\Psi(q)/\Psi(r)$ to get
      // $\nabla\Psi(q)/\Psi(r)
      gradtmp_ *= psiratio[j];
#if defined(QMC_COMPLEX)
      //And now we take the real part and save it.
      convertToReal(gradtmp_, gradpsiratio[j]);
#else
      //Real nonlocalpp forces seem to differ from those in the complex build.  Since
      //complex build has been validated against QE, that indicates there's a bug for the real build.
      gradpsiratio[j] = gradtmp_;
#endif
      W.rejectMove(iel);
      psi.resetPhaseDiff();
      //psi.rejectMove(iel);
    }
  }

  for (int j = 0; j < nknot; j++)
    psiratio[j] *= sgridweight_m[j];

  // This is just a temporary variable to dump d2/dr2 into for spline evaluation.
  RealType secondderiv(0);

  const RealType rinv = cone / r;

  // Compute radial potential and its derivative times (2l+1)
  for (int ip = 0; ip < nchannel; ip++)
  {
    //fun fact.  NLPComponent stores v(r) as v(r), and not as r*v(r) like in other places.
    vrad[ip]  = nlpp_m[ip]->splint(r, dvrad[ip], secondderiv) * wgt_angpp_m[ip];
    vgrad[ip] = dvrad[ip] * dr * wgt_angpp_m[ip] * rinv;
  }

  //Now to construct the 3N dimensional ionic wfn derivatives for pulay terms.
  //This is going to be slow an painful for now.
  for (size_t jat = 0; jat < ions.getTotalNum(); jat++)
  {
    pulay_ref[jat] = psi.evalGradSource(W, ions, jat);
    gradpotterm_   = 0;
    for (size_t j = 0; j < nknot; j++)
    {
      deltaV[j] = r * rrotsgrid_m[j] - dr;
      //This sequence is necessary to update the distance tables and make the
      //inverse matrix available for force computation.  Move the particle to
      //quadrature point...
      W.makeMove(iel, deltaV[j]);
      psi.calcRatio(W, iel);
      psi.acceptMove(W, iel);
      W.acceptMove(iel); // it only updates the jel-th row of e-e table
      //Done with the move.  Ready for force computation.

      iongradtmp_ = psi.evalGradSource(W, ions, jat);
      iongradtmp_ *= psiratio[j];
#ifdef QMC_COMPLEX
      convertToReal(iongradtmp_, pulay_quad[j][jat]);
#endif
      pulay_quad[j][jat] = iongradtmp_;
      //And move the particle back.
      deltaV[j] = dr - r * rrotsgrid_m[j];

      // mirror the above in reverse order
      W.makeMove(iel, deltaV[j]);
      psi.calcRatio(W, iel);
      psi.acceptMove(W, iel);
      W.acceptMove(iel);
    }
  }

  RealType pairpot = 0;
  // Compute spherical harmonics on grid
  for (int j = 0; j < nknot; j++)
  {
    RealType zz        = dot(dr, rrotsgrid_m[j]) * rinv;
    PosType uminusrvec = rrotsgrid_m[j] - zz * dr * rinv;

    cosgrad[j] = rinv * uminusrvec;

    RealType udotgradpsi = dot(gradpsiratio[j], rrotsgrid_m[j]);
    wfngrad[j]           = gradpsiratio[j] - dr * (udotgradpsi * rinv);
    wfngrad[j] *= sgridweight_m[j];

    // Forming the Legendre polynomials
    //P_0(x)=1; P'_0(x)=0.
    lpol[0]  = cone;
    dlpol[0] = czero;

    RealType lpolprev  = czero;
    RealType dlpolprev = czero;

    for (int l = 0; l < lmax; l++)
    {
      //Legendre polynomial recursion formula.
      lpol[l + 1] = Lfactor1[l] * zz * lpol[l] - l * lpolprev;
      lpol[l + 1] *= Lfactor2[l];

      //and for the derivative...
      dlpol[l + 1] = Lfactor1[l] * (zz * dlpol[l] + lpol[l]) - l * dlpolprev;
      dlpol[l + 1] *= Lfactor2[l];

      lpolprev  = lpol[l];
      dlpolprev = dlpol[l];
    }

    RealType lsum = czero;
    // Now to compute the forces:
    gradpotterm_   = 0;
    gradlpolyterm_ = 0;
    gradwfnterm_   = 0;
    pulaytmp_      = 0;

    for (int l = 0; l < nchannel; l++)
    {
      //Note.  Because we are computing "forces", there's a -1 difference between this and
      //direct finite difference calculations.
      lsum += std::real(vrad[l]) * lpol[angpp_m[l]];
      gradpotterm_ += vgrad[l] * lpol[angpp_m[l]] * std::real(psiratio[j]);
      gradlpolyterm_ += std::real(vrad[l]) * dlpol[angpp_m[l]] * cosgrad[j] * std::real(psiratio[j]);
      gradwfnterm_ += std::real(vrad[l]) * lpol[angpp_m[l]] * wfngrad[j];
      pulaytmp_ -= std::real(vrad[l]) * lpol[angpp_m[l]] * pulay_quad[j];
    }
    knot_pots[j] = std::real(lsum * psiratio[j]);
    pulaytmp_ += knot_pots[j] * pulay_ref;
    pairpot += knot_pots[j];
    force_iat += gradpotterm_ + gradlpolyterm_ - gradwfnterm_;
    pulay_terms += pulaytmp_;
  }

  return pairpot;
}

///Randomly rotate sgrid_m
void NonLocalECPComponent::randomize_grid(RandomGenerator& myRNG)
{
  RealType phi(TWOPI * myRNG()), psi(TWOPI * myRNG()), cth(myRNG() - 0.5);
  RealType sph(std::sin(phi)), cph(std::cos(phi)), sth(std::sqrt(1.0 - cth * cth)), sps(std::sin(psi)),
      cps(std::cos(psi));
  TensorType rmat(cph * cth * cps - sph * sps, sph * cth * cps + cph * sps, -sth * cps, -cph * cth * sps - sph * cps,
                  -sph * cth * sps + cph * cps, sth * sps, cph * sth, sph * sth, cth);
  for (int i = 0; i < sgridxyz_m.size(); i++)
    rrotsgrid_m[i] = dot(rmat, sgridxyz_m[i]);
}

template<typename T>
void NonLocalECPComponent::randomize_grid(std::vector<T>& sphere, RandomGenerator& myRNG)
{
  RealType phi(TWOPI * myRNG()), psi(TWOPI * myRNG()), cth(myRNG() - 0.5);
  RealType sph(std::sin(phi)), cph(std::cos(phi)), sth(std::sqrt(1.0 - cth * cth)), sps(std::sin(psi)),
      cps(std::cos(psi));
  TensorType rmat(cph * cth * cps - sph * sps, sph * cth * cps + cph * sps, -sth * cps, -cph * cth * sps - sph * cps,
                  -sph * cth * sps + cph * cps, sth * sps, cph * sth, sph * sth, cth);
  SpherGridType::iterator it(sgridxyz_m.begin());
  SpherGridType::iterator it_end(sgridxyz_m.end());
  SpherGridType::iterator jt(rrotsgrid_m.begin());
  while (it != it_end)
  {
    *jt = dot(rmat, *it);
    ++it;
    ++jt;
  }
  //copy the randomized grid to sphere
  for (int i = 0; i < rrotsgrid_m.size(); i++)
    for (int j = 0; j < OHMMS_DIM; j++)
      sphere[OHMMS_DIM * i + j] = rrotsgrid_m[i][j];
}

void NonLocalECPComponent::buildQuadraturePointDeltaPositions(RealType r,
                                                              const PosType& dr,
                                                              std::vector<PosType>& deltaV) const
{
  for (int j = 0; j < nknot; j++)
    deltaV[j] = r * rrotsgrid_m[j] - dr;
}

void NonLocalECPComponent::contributeTxy(int iel, std::vector<NonLocalData>& Txy) const
{
  for (int j = 0; j < nknot; j++)
    Txy.push_back(NonLocalData(iel, knot_pots[j], deltaV[j]));
}

/// \relates NonLocalEcpComponent
template void NonLocalECPComponent::randomize_grid(std::vector<float>& sphere, RandomGenerator& myRNG);
template void NonLocalECPComponent::randomize_grid(std::vector<double>& sphere, RandomGenerator& myRNG);


} // namespace qmcplusplus
