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
NonLocalECPComponent::NonLocalECPComponent() : lmax(0), nchannel(0), nknot(0), Rmax(-1), VP(nullptr), do_randomize_grid_(true) {}

// unfortunately we continue the sloppy use of the default copy constructor followed by reassigning pointers.
// This prevents use of smart pointers and concievably sets us up for trouble with double frees and the destructor.
NonLocalECPComponent::NonLocalECPComponent(const NonLocalECPComponent& nl_ecpc, const ParticleSet& pset)
    : NonLocalECPComponent(nl_ecpc)
{
  for (int i = 0; i < nl_ecpc.nlpp_m.size(); ++i)
    nlpp_m[i] = nl_ecpc.nlpp_m[i]->makeClone();
  if (nl_ecpc.VP)
    VP = new VirtualParticleSet(pset, nknot, nl_ecpc.VP->getNumDistTables());
}

NonLocalECPComponent::~NonLocalECPComponent()
{
  for (int ip = 0; ip < nlpp_m.size(); ip++)
    delete nlpp_m[ip];
  if (VP)
    delete VP;
}

void NonLocalECPComponent::set_randomize_grid(bool do_randomize_grid)
{
  do_randomize_grid_ = do_randomize_grid;
}

void NonLocalECPComponent::initVirtualParticle(const ParticleSet& qp)
{
  assert(VP == 0);
  outputManager.pause();
  VP = new VirtualParticleSet(qp, nknot);
  outputManager.resume();
}

void NonLocalECPComponent::deleteVirtualParticle()
{
  if (VP)
    delete VP;
  VP = nullptr;
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
  wvec.resize(n);
  lpol.resize(l + 1, 1.0);
  //dlpol needs two data points to do a recursive construction.  Also is only nontrivial for l>1.
  //This +2 guards against l=0 case.
  dlpol.resize(l + 2, 0.0);
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
  os << "    Maximum angular momentum = " << lmax << std::endl;
  os << "    Number of non-local channels = " << nchannel << std::endl;
  for (int l = 0; l < nchannel; l++)
    os << "       l(" << l << ")=" << angpp_m[l] << std::endl;
  os << "    Cutoff radius = " << Rmax << std::endl;
  os << "    Number of spherical integration grid points = " << nknot << std::endl;
  if (outputManager.isActive(Verbosity::HIGH))
  {
    os << "    Spherical grid and weights: " << std::endl;
    for (int ik = 0; ik < nknot; ik++)
      os << "       " << sgridxyz_m[ik] << std::setw(20) << sgridweight_m[ik] << std::endl;
  }
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
    VP->makeMoves(W, iel, deltaV, true, iat);
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
      lpol[l + 1] = (Lfactor1[l] * zz * lpol[l] - l * lpolprev) * Lfactor2[l];
      lpolprev    = lpol[l];
    }

    RealType lsum = 0.0;
    for (int l = 0; l < nchannel; l++)
      lsum += vrad[l] * lpol[angpp_m[l]];
    knot_pots[j] = lsum * std::real(psiratio[j]);
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

    VirtualParticleSet::mw_makeMoves(vp_list, p_list, deltaV_list, joblist, true);

    if (use_DLA)
      TrialWaveFunction::mw_evaluateRatios(psi_list, const_vp_list, psiratios_list,
                                           TrialWaveFunction::ComputeType::FERMIONIC);
    else
      TrialWaveFunction::mw_evaluateRatios(psi_list, const_vp_list, psiratios_list);
  }
  else
  {
    // Compute ratios without VP. This is working but very slow code path.
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
    VP->makeMoves(W, iel, deltaV, true, iat);
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
      lpol[l + 1] = (Lfactor1[l] * zz * lpol[l] - l * lpolprev) * Lfactor2[l];

      //and for the derivative...
      dlpol[l + 1] = (Lfactor1[l] * (zz * dlpol[l] + lpol[l]) - l * dlpolprev) * Lfactor2[l];

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
      lsum += vrad[l] * lpol[angpp_m[l]];
      gradpotterm_ += vgrad[l] * lpol[angpp_m[l]] * std::real(psiratio[j]);
      gradlpolyterm_ += vrad[l] * dlpol[angpp_m[l]] * cosgrad[j] * std::real(psiratio[j]);
      gradwfnterm_ += vrad[l] * lpol[angpp_m[l]] * wfngrad[j];
    }
    knot_pots[j] = lsum * std::real(psiratio[j]);
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
    VP->makeMoves(W, iel, deltaV, true, iat);
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
    convertToReal(psi.evalGradSource(W, ions, jat), pulay_ref[jat]);
    gradpotterm_ = 0;
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
      convertToReal(iongradtmp_, pulay_quad[j][jat]);
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
      lpol[l + 1] = (Lfactor1[l] * zz * lpol[l] - l * lpolprev) * Lfactor2[l];

      //and for the derivative...
      dlpol[l + 1] = (Lfactor1[l] * (zz * dlpol[l] + lpol[l]) - l * dlpolprev) * Lfactor2[l];

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
      lsum += vrad[l] * lpol[angpp_m[l]];
      gradpotterm_ += vgrad[l] * lpol[angpp_m[l]] * std::real(psiratio[j]);
      gradlpolyterm_ += vrad[l] * dlpol[angpp_m[l]] * cosgrad[j] * std::real(psiratio[j]);
      gradwfnterm_ += vrad[l] * lpol[angpp_m[l]] * wfngrad[j];
      pulaytmp_ -= vrad[l] * lpol[angpp_m[l]] * pulay_quad[j];
    }
    knot_pots[j] = lsum * std::real(psiratio[j]);
    pulaytmp_ += knot_pots[j] * pulay_ref;
    pairpot += knot_pots[j];
    force_iat += gradpotterm_ + gradlpolyterm_ - gradwfnterm_;
    pulay_terms += pulaytmp_;
  }

  return pairpot;
}

void NonLocalECPComponent::evaluateOneBodyOpMatrixContribution(ParticleSet& W,
                                                               const int iat,
                                                               const TWFFastDerivWrapper& psi,
                                                               const int iel,
                                                               const RealType r,
                                                               const PosType& dr,
                                                               std::vector<ValueMatrix>& B)

{
  using ValueVector = SPOSet::ValueVector;
  //Some initial computation to find out the species and row number of electron.
  const IndexType gid        = W.getGroupID(iel);
  const IndexType sid        = psi.getTWFGroupIndex(gid);
  const IndexType firstIndex = W.first(gid);
  const IndexType thisIndex  = iel - firstIndex;

  const IndexType numOrbs = psi.numOrbitals(sid);
  ValueVector phi_row; //phi_0(r), phi_1(r), ...
  ValueVector temp_row;

  phi_row.resize(numOrbs);
  temp_row.resize(numOrbs);

  buildQuadraturePointDeltaPositions(r, dr, deltaV);


  constexpr RealType czero(0);
  constexpr RealType cone(1);

  const RealType rinv = cone / r;

  for (int ip = 0; ip < nchannel; ip++)
    vrad[ip] = nlpp_m[ip]->splint(r) * wgt_angpp_m[ip];

  for (int j = 0; j < nknot; j++)
  {
    W.makeMove(iel, deltaV[j], false); //Update distance tables.
    psi.getRowM(W, iel, phi_row);
    RealType jratio = psi.evaluateJastrowRatio(W, iel);
    W.rejectMove(iel);

    RealType zz = dot(dr, rrotsgrid_m[j]) * rinv;
    // Forming the Legendre polynomials
    lpol[0]           = cone;
    RealType lpolprev = czero;
    for (int l = 0; l < lmax; l++)
    {
      lpol[l + 1] = Lfactor2[l] * (Lfactor1[l] * zz * lpol[l] - l * lpolprev);
      lpolprev    = lpol[l];
    }

    ValueType lsum = 0.0;
    for (int l = 0; l < nchannel; l++)
    {
      temp_row = (vrad[l] * lpol[angpp_m[l]] * sgridweight_m[j]) * jratio * phi_row;
      for (int iorb = 0; iorb < numOrbs; iorb++)
        B[sid][thisIndex][iorb] += temp_row[iorb];
    }
  }
}

void NonLocalECPComponent::evaluateOneBodyOpMatrixdRContribution(ParticleSet& W,
                                                                 ParticleSet& ions,
                                                                 const int iat,
                                                                 const int iat_src,
                                                                 const TWFFastDerivWrapper& psi,
                                                                 const int iel,
                                                                 const RealType r,
                                                                 const PosType& dr,
                                                                 std::vector<std::vector<ValueMatrix>>& dB)
{
  using ValueVector = SPOSet::ValueVector;
  using GradVector  = SPOSet::GradVector;
  constexpr RealType czero(0);
  constexpr RealType cone(1);

  //We check that our quadrature grid is valid.  Namely, that all points lie on the unit sphere.
  //We check this by seeing if |r|^2 = 1 to machine precision.
  for (int j = 0; j < nknot; j++)
    assert(std::abs(std::sqrt(dot(rrotsgrid_m[j], rrotsgrid_m[j])) - 1) <
           100 * std::numeric_limits<RealType>::epsilon());

  for (int j = 0; j < nknot; j++)
    deltaV[j] = r * rrotsgrid_m[j] - dr;

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

  const IndexType gid        = W.getGroupID(iel);
  const IndexType sid        = psi.getTWFGroupIndex(gid);
  const IndexType thisEIndex = iel - W.first(gid);
  const IndexType numptcls   = W.last(gid) - W.first(gid);
  const IndexType norbs      = psi.numOrbitals(sid);
  SPOSet& spo(*psi.getSPOSet(sid));

  RealType pairpot = 0;
  // Compute spherical harmonics on grid
  GradMatrix iongrad_phimat;
  GradVector iongrad_phi;
  ValueVector phi, laplphi;
  ValueMatrix phimat, laplphimat;
  GradMatrix gradphimat;
  GradVector gradphi;
  GradMatrix wfgradmat;

  //This stores all ion gradients of the Jastrow, evaluated on all quadrature points.
  //Has length nknot.
  ParticleSet::ParticleGradient jgrad_quad;
  jgrad_quad.resize(nknot);

  wfgradmat.resize(nknot, norbs);

  iongrad_phimat.resize(nknot, norbs);
  laplphimat.resize(nknot, norbs);

  phimat.resize(nknot, norbs);
  gradphimat.resize(nknot, norbs);
  iongrad_phi.resize(norbs);
  phi.resize(norbs);
  gradphi.resize(norbs);
  laplphi.resize(norbs);

  GradVector udotgradpsi, wfgradrow;
  GradMatrix udotgradpsimat;
  GradVector gpot, glpoly, gwfn;
  ValueVector nlpp_prefactor;
  GradVector dlpoly_prefactor;
  GradVector dvdr_prefactor;

  nlpp_prefactor.resize(nknot);
  dlpoly_prefactor.resize(nknot);
  dvdr_prefactor.resize(nknot);

  udotgradpsimat.resize(nknot, norbs);
  wfgradrow.resize(norbs);
  udotgradpsi.resize(norbs);
  gpot.resize(norbs);
  glpoly.resize(norbs);
  gwfn.resize(norbs);

  buildQuadraturePointDeltaPositions(r, dr, deltaV);
  //This is the ion gradient of J at the original (non quadrature) coordinate.
  GradType jigradref(0.0);

  jigradref = psi.evaluateJastrowGradSource(W, ions, iat_src);

  //Until we have a more efficient routine, we move to a quadrature point,
  //update distance tables, compute the ion gradient of J, then move the particle back.
  //At cost of distance table updates.  Not good, but works.
  for (int j = 0; j < nknot; j++)
  {
    W.makeMove(iel, deltaV[j], false);
    W.acceptMove(iel);
    jgrad_quad[j] = psi.evaluateJastrowGradSource(W, ions, iat_src);
    W.makeMove(iel, -deltaV[j], false);
    W.acceptMove(iel);
  }

  for (int j = 0; j < nknot; j++)
  {
    W.makeMove(iel, deltaV[j], false);
    iongrad_phi = 0.0;
    spo.evaluateGradSourceRow(W, iel, ions, iat_src, iongrad_phi);
    GradType jegrad(0.0);
    GradType jigrad(0.0);

    RealType jratio = psi.calcJastrowRatioGrad(W, iel, jegrad);
    jigrad          = psi.evaluateJastrowGradSource(W, ions, iat_src);

    spo.evaluateVGL(W, iel, phi, gradphi, laplphi);

    //Quick comment on the matrix elements being computed below.
    //For the no jastrow implementation, phimat, gradphimat, iongrad_phimat were straightforward containers storing phi_j(r_i), grad(phi_j), etc.
    //Generalizing to jastrows is straightforward if we replace phi_j(q) with exp(J(q))/exp(J(r))*phi(q).  Storing these in the phimat, gradphimat, etc.
    //data structures allows us to not modify the rather complicated expressions we have already derived.
    if (iat == iat_src)
    {
      for (int iorb = 0; iorb < norbs; iorb++)
      {
        //Treating exp(J(q))/exp(J(r))phi_j(q) as the fundamental block.
        phimat[j][iorb] = jratio * phi[iorb];
        //This is the electron gradient of the above expression.
        gradphimat[j][iorb] = jratio * (gradphi[iorb] + GradType(jegrad) * phi[iorb]);
        laplphimat[j][iorb] = laplphi[iorb]; //this is not used, so not including jastrow contribution.
      }
    }
    for (int iorb = 0; iorb < norbs; iorb++)
    {
      //This is the ion gradient of exp(J(q))/exp(J(r))phi_j(q).
      iongrad_phimat[j][iorb] = jratio * (iongrad_phi[iorb] + phi[iorb] * (GradType(jgrad_quad[j]) - jigradref));
    }
    W.rejectMove(iel);
  }

  for (int j = 0; j < nknot; j++)
  {
    RealType zz        = dot(dr, rrotsgrid_m[j]) * rinv;
    PosType uminusrvec = rrotsgrid_m[j] - zz * dr * rinv;

    cosgrad[j] = rinv * uminusrvec;

    // Forming the Legendre polynomials
    //P_0(x)=1; P'_0(x)=0.
    lpol[0]  = cone;
    dlpol[0] = czero;
    dlpol[1] = cone;

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

    for (int l = 0; l < nchannel; l++)
    {
      //Note.  Because we are computing "forces", there's a -1 difference between this and
      //direct finite difference calculations.

      nlpp_prefactor[j] += sgridweight_m[j] * vrad[l] * lpol[angpp_m[l]];
      dvdr_prefactor[j] += sgridweight_m[j] * vgrad[l] * lpol[angpp_m[l]];
      dlpoly_prefactor[j] += sgridweight_m[j] * vrad[l] * dlpol[angpp_m[l]] * cosgrad[j];
    }
  }


  for (int j = 0; j < nknot; j++)
    for (int iorb = 0; iorb < norbs; iorb++)
      gwfn[iorb] += nlpp_prefactor[j] * (iongrad_phimat[j][iorb]);

  if (iat == iat_src)
  {
    for (int j = 0; j < nknot; j++)
      for (int iorb = 0; iorb < norbs; iorb++)
      {
        //this is for diagonal case.
        udotgradpsimat[j][iorb] = dot(gradphimat[j][iorb], rrotsgrid_m[j]);
        wfgradmat[j][iorb]      = gradphimat[j][iorb] - dr * (udotgradpsimat[j][iorb] * rinv);
        gpot[iorb] += dvdr_prefactor[j] * phimat[j][iorb];
        glpoly[iorb] += dlpoly_prefactor[j] * phimat[j][iorb];
        gwfn[iorb] += nlpp_prefactor[j] * (wfgradmat[j][iorb]);
      }
  }
  for (int idim = 0; idim < OHMMS_DIM; idim++)
    for (int iorb = 0; iorb < norbs; iorb++)
      dB[idim][sid][thisEIndex][iorb] += RealType(-1.0) * gpot[iorb][idim] - glpoly[iorb][idim] + gwfn[iorb][idim];
}

///Randomly rotate sgrid_m
void NonLocalECPComponent::rotateQuadratureGrid(const TensorType& rmat)
{
  for (int i = 0; i < sgridxyz_m.size(); i++)
    if (do_randomize_grid_)
      rrotsgrid_m[i] = dot(rmat, sgridxyz_m[i]);
    else
      rrotsgrid_m[i] = sgridxyz_m[i];
}

template<typename T>
void NonLocalECPComponent::rotateQuadratureGrid(std::vector<T>& sphere, const TensorType& rmat)
{
  SpherGridType::iterator it(sgridxyz_m.begin());
  SpherGridType::iterator it_end(sgridxyz_m.end());
  SpherGridType::iterator jt(rrotsgrid_m.begin());
  while (it != it_end)
  {
    if (do_randomize_grid_)
      *jt = dot(rmat, *it);
    else
      *jt = *it;
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
template void NonLocalECPComponent::rotateQuadratureGrid(std::vector<float>& sphere, const TensorType& rmat);
template void NonLocalECPComponent::rotateQuadratureGrid(std::vector<double>& sphere, const TensorType& rmat);


} // namespace qmcplusplus
