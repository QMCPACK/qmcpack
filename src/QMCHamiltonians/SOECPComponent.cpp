//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Raymond C. Clay, rclay@sandia.gov, Sandia National Laboratories
//                    Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond C. Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/DistanceTable.h"
#include "SOECPComponent.h"
#include "Numerics/Ylm.h"
#include "CPU/BLAS.hpp"
#include "NLPPJob.h"

namespace qmcplusplus
{
SOECPComponent::SOECPComponent()
    : lmax_(0), nchannel_(0), nknot_(0), sknot_(0), total_knots_(0), rmax_(-1), vp_(nullptr)
{}

SOECPComponent::~SOECPComponent()
{
  for (int i = 0; i < sopp_m_.size(); i++)
    delete sopp_m_[i];
  if (vp_)
    delete vp_;
}

void SOECPComponent::print(std::ostream& os) {}

void SOECPComponent::initVirtualParticle(const ParticleSet& qp)
{
  assert(vp_ == nullptr);
  outputManager.pause();
  vp_ = new VirtualParticleSet(qp, total_knots_);
  outputManager.resume();
}

void SOECPComponent::deleteVirtualParticle()
{
  if (vp_)
    delete vp_;
  vp_ = nullptr;
}

void SOECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m_.push_back(l);
  sopp_m_.push_back(pp);
}

SOECPComponent* SOECPComponent::makeClone(const ParticleSet& qp)
{
  SOECPComponent* myclone = new SOECPComponent(*this);
  for (int i = 0; i < sopp_m_.size(); i++)
    myclone->sopp_m_[i] = sopp_m_[i]->makeClone();
  if (vp_)
    myclone->vp_ = new VirtualParticleSet(qp, total_knots_);
  return myclone;
}

void SOECPComponent::resize_warrays(int n, int m, int s)
{
  vrad_.resize(m);
  rrotsgrid_m_.resize(n);
  nchannel_ = sopp_m_.size();
  nknot_    = sgridxyz_m_.size();
  sknot_    = s;
  if (sknot_ % 2 != 0 && sknot_ > 0)
    throw std::runtime_error("Spin integration uses Simpson's rule. Must have even number of knots in this case");

  //Need +1 for Simpsons rule to include both end points.
  //sknot here refers to the number of subintervals for integration
  total_knots_ = nknot_ * (sknot_ + 1);
  psiratio_.resize(total_knots_);
  deltaV_.resize(total_knots_);
  deltaS_.resize(total_knots_);
  spin_quad_weights_.resize(total_knots_);
  auto& up_row = spinor_multiplier_.first;
  auto& dn_row = spinor_multiplier_.second;
  up_row.resize(total_knots_);
  dn_row.resize(total_knots_);
  if (m != nchannel_)
  {
    APP_ABORT("SOECPComponent::resize_warrays has incorrect number of radial channels\n");
  }
}

int SOECPComponent::kroneckerDelta(int x, int y) { return (x == y) ? 1 : 0; }

SOECPComponent::ComplexType SOECPComponent::sMatrixElements(RealType s1, RealType s2, int dim)
{
  switch (dim)
  {
  case 0:
    return ComplexType(std::cos(s1 + s2), 0.0);
  case 1:
    return ComplexType(std::sin(s1 + s2), 0.0);
  case 2:
    return ComplexType(0.0, std::sin(s1 - s2));
  default:
    APP_ABORT("SOECPComponent::sMatrixElements invalid operator dimension\n");
    return 0;
  }
}

SOECPComponent::ComplexType SOECPComponent::lmMatrixElements(int l, int m1, int m2, int dim)
{
  ComplexType val;
  RealType onehalf = 0.5;
  RealType zero    = 0.0;
  switch (dim)
  {
  case 0: //x
    val = onehalf *
        ComplexType(std::sqrt(l * (l + 1) - m2 * (m2 + 1)) * kroneckerDelta(m1, m2 + 1) +
                        std::sqrt(l * (l + 1) - m2 * (m2 - 1)) * kroneckerDelta(m1, m2 - 1),
                    zero);
    return val;
  case 1:
    val = onehalf *
        ComplexType(0.0,
                    std::sqrt(l * (l + 1) - m2 * (m2 - 1)) * kroneckerDelta(m1, m2 - 1) -
                        std::sqrt(l * (l + 1) - m2 * (m2 + 1)) * kroneckerDelta(m1, m2 + 1));
    return val;
  case 2:
    val = ComplexType(m2 * kroneckerDelta(m1, m2), zero);
    return val;
  default:
    APP_ABORT("SOECPComponent::lMatrixElements Invalid operator dimension\n");
    return 0;
  }
}

SOECPComponent::ComplexType SOECPComponent::matrixElementDecomposed(int l, int m1, int m2, RealType spin, bool plus)
{
  RealType pm = plus ? RealType(1) : RealType(-1);
  RealType mp = -1 * pm;


  ComplexType lx = lmMatrixElements(l, m1, m2, 0);
  ComplexType ly = lmMatrixElements(l, m1, m2, 1);
  ComplexType lz = lmMatrixElements(l, m1, m2, 2);

  RealType cos2s = std::cos(2 * spin);
  RealType sin2s = std::sin(2 * spin);
  ComplexType phase(cos2s, mp * sin2s);
  ComplexType eye(0, 1);
  RealType onehalf = 0.5;
  ComplexType val  = onehalf * (pm * lz + phase * (lx + pm * eye * ly));
  return val;
}

SOECPComponent::RealType SOECPComponent::evaluateOne(ParticleSet& W,
                                                     int iat,
                                                     TrialWaveFunction& Psi,
                                                     int iel,
                                                     RealType r,
                                                     const PosType& dr)
{
  RealType sold = W.spins[iel];
  buildTotalQuadrature(r, dr, sold);

  if (vp_)
  {
    vp_->makeMovesWithSpin(W, iel, deltaV_, deltaS_, true, iat);
    Psi.evaluateRatios(*vp_, psiratio_);
  }
  else
    for (int iq = 0; iq < total_knots_; iq++)
    {
      W.makeMoveWithSpin(iel, deltaV_[iq], deltaS_[iq]);
      psiratio_[iq] = Psi.calcRatio(W, iel);
      W.rejectMove(iel);
      Psi.resetPhaseDiff();
    }

  return calculateProjector(r, dr, sold);
}

SOECPComponent::RealType SOECPComponent::calculateProjector(RealType r, const PosType& dr, RealType sold)
{
#ifdef QMC_COMPLEX
  wvec_.resize(total_knots_); //contribution from each quarature point
  ComplexType pairpot;
  for (int iq = 0; iq < total_knots_; iq++)
  {
    RealType snew = sold + deltaS_[iq];
    ComplexType lsum;
    for (int il = 0; il < nchannel_; il++)
    {
      int l = il + 1; //nchannels starts at l=1, so 0th element is p not s
      ComplexType msums;
      for (int m1 = -l; m1 <= l; m1++)
      {
        ComplexType Y = sphericalHarmonic(l, m1, dr);
        for (int m2 = -l; m2 <= l; m2++)
        {
          ComplexType ldots;
          for (int id = 0; id < 3; id++)
            ldots += lmMatrixElements(l, m1, m2, id) * sMatrixElements(sold, snew, id);
          ComplexType cY = std::conj(sphericalHarmonic(l, m2, rrotsgrid_m_[iq % nknot_]));
          msums += Y * cY * ldots;
        }
      }
      lsum += vrad_[il] * msums;
    }
    wvec_[iq] = psiratio_[iq] * lsum * spin_quad_weights_[iq];
    pairpot += wvec_[iq];
  }
  return std::real(pairpot);
#else
  throw std::runtime_error("SOECPComponent::calculateProjector only implemented in complex build.");
#endif
}

void SOECPComponent::setupExactSpinProjector(RealType r, const PosType& dr, RealType sold)
{
#ifdef QMC_COMPLEX
  auto& up_row = spinor_multiplier_.first;
  auto& dn_row = spinor_multiplier_.second;

  //calculate radial potentials
  for (int ip = 0; ip < nchannel_; ip++)
    vrad_[ip] = sopp_m_[ip]->splint(r);

  for (int iq = 0; iq < nknot_; iq++)
  {
    up_row[iq] = 0.0;
    dn_row[iq] = 0.0;
    for (int il = 0; il < nchannel_; il++)
    {
      int l = il + 1;
      for (int m1 = -l; m1 <= l; m1++)
      {
        ComplexType Y = sphericalHarmonic(l, m1, dr);
        for (int m2 = -l; m2 <= l; m2++)
        {
          ComplexType cY    = std::conj(sphericalHarmonic(l, m2, rrotsgrid_m_[iq]));
          ComplexType so_up = matrixElementDecomposed(l, m1, m2, sold, true);
          ComplexType so_dn = matrixElementDecomposed(l, m1, m2, sold, false);
          RealType fourpi   = 4.0 * M_PI;
          //Note: Need 4pi weight. In Normal NonLocalECP, 1/4Pi generated from transformation to legendre polynomials and gets absorbed into the
          // quadrature integration. We don't get the 1/4Pi from legendre here, so we need to scale by 4Pi.
          up_row[iq] += fourpi * vrad_[il] * Y * cY * so_up * sgridweight_m_[iq];
          dn_row[iq] += fourpi * vrad_[il] * Y * cY * so_dn * sgridweight_m_[iq];
        }
      }
    }
  }
#else
  throw std::runtime_error("SOECPComponent::setupExactSpinIntegration only implemented in complex build.");
#endif
}

SOECPComponent::RealType SOECPComponent::evaluateOneExactSpinIntegration(ParticleSet& W,
                                                                         const int iat,
                                                                         const TrialWaveFunction& psi,
                                                                         const int iel,
                                                                         const RealType r,
                                                                         const PosType& dr)
{
  RealType sold = W.spins[iel];

  for (int j = 0; j < nknot_; j++)
    deltaV_[j] = r * rrotsgrid_m_[j] - dr;
  vp_->makeMoves(W, iel, deltaV_, true, iat);

  setupExactSpinProjector(r, dr, sold);

  psi.evaluateSpinorRatios(*vp_, spinor_multiplier_, psiratio_);
  ComplexType pairpot;
  for (size_t iq = 0; iq < nknot_; iq++)
    pairpot += psiratio_[iq];
  return std::real(pairpot);
}

void SOECPComponent::mw_evaluateOne(const RefVectorWithLeader<SOECPComponent>& soecp_component_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    const RefVectorWithLeader<TrialWaveFunction>& psi_list,
                                    const RefVector<const NLPPJob<RealType>>& joblist,
                                    std::vector<RealType>& pairpots,
                                    ResourceCollection& collection)
{
  auto& soecp_component_leader = soecp_component_list.getLeader();
  if (soecp_component_leader.vp_)
  {
    // Compute ratios with VP
    RefVectorWithLeader<VirtualParticleSet> vp_list(*soecp_component_leader.vp_);
    RefVectorWithLeader<const VirtualParticleSet> const_vp_list(*soecp_component_leader.vp_);
    RefVector<const std::vector<PosType>> deltaV_list;
    RefVector<const std::vector<RealType>> deltaS_list;
    RefVector<std::vector<ValueType>> psiratios_list;
    vp_list.reserve(soecp_component_list.size());
    const_vp_list.reserve(soecp_component_list.size());
    deltaV_list.reserve(soecp_component_list.size());
    deltaS_list.reserve(soecp_component_list.size());
    psiratios_list.reserve(soecp_component_list.size());

    for (size_t i = 0; i < soecp_component_list.size(); i++)
    {
      SOECPComponent& component(soecp_component_list[i]);
      const NLPPJob<RealType>& job = joblist[i];
      const RealType sold          = p_list[i].spins[job.electron_id];

      component.buildTotalQuadrature(job.ion_elec_dist, job.ion_elec_displ, sold);

      vp_list.push_back(*component.vp_);
      const_vp_list.push_back(*component.vp_);
      deltaV_list.push_back(component.deltaV_);
      deltaS_list.push_back(component.deltaS_);
      psiratios_list.push_back(component.psiratio_);
    }

    ResourceCollectionTeamLock<VirtualParticleSet> vp_res_lock(collection, vp_list);

    VirtualParticleSet::mw_makeMovesWithSpin(vp_list, p_list, deltaV_list, deltaS_list, joblist, true);

    TrialWaveFunction::mw_evaluateRatios(psi_list, const_vp_list, psiratios_list);
  }
  else
  {
    // Compute ratios without VP. Slow
    for (size_t i = 0; i < p_list.size(); i++)
    {
      SOECPComponent& component(soecp_component_list[i]);
      ParticleSet& W(p_list[i]);
      TrialWaveFunction& psi(psi_list[i]);
      const NLPPJob<RealType>& job = joblist[i];

      const RealType sold = W.spins[job.electron_id];
      component.buildTotalQuadrature(job.ion_elec_dist, job.ion_elec_displ, sold);

      for (int j = 0; j < component.total_knots_; j++)
      {
        W.makeMoveWithSpin(job.electron_id, component.deltaV_[j], component.deltaS_[j]);
        component.psiratio_[j] = psi.calcRatio(W, job.electron_id);
        W.rejectMove(job.electron_id);
        psi.resetPhaseDiff();
      }
    }
  }

  for (size_t i = 0; i < p_list.size(); i++)
  {
    SOECPComponent& component(soecp_component_list[i]);
    const NLPPJob<RealType>& job = joblist[i];
    const RealType sold          = p_list[i].spins[job.electron_id];
    pairpots[i]                  = component.calculateProjector(job.ion_elec_dist, job.ion_elec_displ, sold);
  }
}

void SOECPComponent::mw_evaluateOneExactSpinIntegration(const RefVectorWithLeader<SOECPComponent>& soecp_component_list,
                                                        const RefVectorWithLeader<ParticleSet>& p_list,
                                                        const RefVectorWithLeader<TrialWaveFunction>& psi_list,
                                                        const RefVector<const NLPPJob<RealType>>& joblist,
                                                        std::vector<RealType>& pairpots,
                                                        ResourceCollection& collection)
{
#ifdef QMC_COMPLEX
  auto& soecp_component_leader = soecp_component_list.getLeader();
  // Compute ratios with VP
  RefVectorWithLeader<VirtualParticleSet> vp_list(*soecp_component_leader.vp_);
  RefVectorWithLeader<const VirtualParticleSet> const_vp_list(*soecp_component_leader.vp_);
  RefVector<const std::vector<PosType>> deltaV_list;
  RefVector<std::vector<ValueType>> psiratios_list;
  RefVector<std::pair<SPOSet::ValueVector, SPOSet::ValueVector>> spinor_multiplier_list;
  vp_list.reserve(soecp_component_list.size());
  const_vp_list.reserve(soecp_component_list.size());
  deltaV_list.reserve(soecp_component_list.size());
  psiratios_list.reserve(soecp_component_list.size());
  spinor_multiplier_list.reserve(soecp_component_list.size());

  for (size_t i = 0; i < soecp_component_list.size(); i++)
  {
    SOECPComponent& component(soecp_component_list[i]);
    const NLPPJob<RealType>& job = joblist[i];

    for (int j = 0; j < component.nknot_; j++)
      component.deltaV_[j] = job.ion_elec_dist * component.rrotsgrid_m_[j] - job.ion_elec_displ;

    vp_list.push_back(*component.vp_);
    const_vp_list.push_back(*component.vp_);
    deltaV_list.push_back(component.deltaV_);
    psiratios_list.push_back(component.psiratio_);
    spinor_multiplier_list.push_back(component.spinor_multiplier_);
  }

  ResourceCollectionTeamLock<VirtualParticleSet> vp_res_lock(collection, vp_list);

  VirtualParticleSet::mw_makeMoves(vp_list, p_list, deltaV_list, joblist, true);

  for (size_t i = 0; i < p_list.size(); i++)
  {
    SOECPComponent& component(soecp_component_list[i]);
    const NLPPJob<RealType>& job = joblist[i];
    const RealType sold          = p_list[i].spins[job.electron_id];
    component.setupExactSpinProjector(job.ion_elec_dist, job.ion_elec_displ, sold);
  }

  TrialWaveFunction::mw_evaluateSpinorRatios(psi_list, const_vp_list, spinor_multiplier_list, psiratios_list);

  for (size_t i = 0; i < p_list.size(); i++)
  {
    pairpots[i] = 0;
    SOECPComponent& component(soecp_component_list[i]);
    std::vector<ValueType>& psiratio(psiratios_list[i]);
    for (int iq = 0; iq < component.nknot_; iq++)
      pairpots[i] += std::real(psiratio[iq]);
  }
#else
  throw std::runtime_error("SOECPComponent::mw_evaluateOneExactSpinIntegration only implemented in complex build");
#endif
}

SOECPComponent::RealType SOECPComponent::evaluateValueAndDerivatives(ParticleSet& W,
                                                                     int iat,
                                                                     TrialWaveFunction& Psi,
                                                                     int iel,
                                                                     RealType r,
                                                                     const PosType& dr,
                                                                     const opt_variables_type& optvars,
                                                                     const Vector<ValueType>& dlogpsi,
                                                                     Vector<ValueType>& dhpsioverpsi)
{
#ifndef QMC_COMPLEX
  throw std::runtime_error("SOECPComponent::evaluateValueAndDerivatives should not be called in real build\n");
#else

  const size_t num_vars = optvars.num_active_vars;
  dratio_.resize(total_knots_, num_vars);
  dlogpsi_vp_.resize(dlogpsi.size());

  RealType sold = W.spins[iel];
  buildTotalQuadrature(r, dr, sold);

  //Now we have all the spin and spatial quadrature points acculated to use in evaluation
  //Now we need to obtain dlogpsi and dlogpsi_vp
  if (vp_)
  {
    // Compute ratios with VP
    vp_->makeMovesWithSpin(W, iel, deltaV_, deltaS_, true, iat);
    Psi.evaluateDerivRatios(*vp_, optvars, psiratio_, dratio_);
  }
  else
    for (int iq = 0; iq < total_knots_; iq++)
    {
      PosType posold = W.R[iel];
      W.makeMoveWithSpin(iel, deltaV_[iq], deltaS_[iq]);
      psiratio_[iq] = Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);

      std::fill(dlogpsi_vp_.begin(), dlogpsi_vp_.end(), 0.0);
      Psi.evaluateDerivativesWF(W, optvars, dlogpsi_vp_);
      for (int v = 0; v < dlogpsi_vp_.size(); ++v)
        dratio_(iq, v) = dlogpsi_vp_[v] - dlogpsi[v];

      W.makeMoveWithSpin(iel, -deltaV_[iq], -deltaS_[iq]);
      Psi.calcRatio(W, iel);
      Psi.acceptMove(W, iel);
      W.acceptMove(iel);
    }

  ComplexType pairpot = calculateProjector(r, dr, sold);
  //calculateProjector stores quad points in wvec_
  BLAS::gemv('N', num_vars, total_knots_, 1.0, dratio_.data(), num_vars, wvec_.data(), 1, 1.0, dhpsioverpsi.data(), 1);

  return std::real(pairpot);
#endif
}

SOECPComponent::RealType SOECPComponent::evaluateValueAndDerivativesExactSpinIntegration(
    ParticleSet& W,
    int iat,
    TrialWaveFunction& Psi,
    int iel,
    RealType r,
    const PosType& dr,
    const opt_variables_type& optvars,
    const Vector<ValueType>& dlogpsi,
    Vector<ValueType>& dhpsioverpsi)
{
#ifndef QMC_COMPLEX
  throw std::runtime_error("SOECPComponent::evaluateValueAndDerivatives should not be called in real build\n");
#else
  const size_t num_vars = optvars.num_active_vars;
  dratio_.resize(total_knots_, num_vars);
  dlogpsi_vp_.resize(dlogpsi.size());

  RealType sold = W.spins[iel];
  for (int j = 0; j < nknot_; j++)
    deltaV_[j] = r * rrotsgrid_m_[j] - dr;
  vp_->makeMoves(W, iel, deltaV_, true, iat);
  setupExactSpinProjector(r, dr, sold);
  Psi.evaluateSpinorDerivRatios(*vp_, spinor_multiplier_, optvars, psiratio_, dratio_);

  BLAS::gemv('N', num_vars, total_knots_, 1.0, dratio_.data(), num_vars, psiratio_.data(), 1, 1.0, dhpsioverpsi.data(), 1);

  ComplexType pairpot;
  for (size_t iq = 0; iq < nknot_; iq++)
    pairpot += psiratio_[iq];
  return std::real(pairpot);
#endif
}

void SOECPComponent::buildTotalQuadrature(const RealType r, const PosType& dr, const RealType sold)
{
  int count = 0;
  RealType smin(0.0);
  RealType smax(TWOPI);
  RealType dS = (smax - smin) / sknot_; //step size for spin


  //for a given spin point in the Simpsons integration, this will copy the spatial quadrature points into
  //the global quadrature arrays which includes spin and spatial quadrature.
  //Sets deltaS_, deltaV_, and spin_quad_weights_ accordingly.
  auto addSpatialQuadrature = [&](const int is, const RealType r, const PosType& dr, const RealType ds,
                                  const RealType spin_weight) {
    for (int iq = 0; iq < nknot_; iq++)
    {
      int offset      = is * nknot_ + iq;
      deltaV_[offset] = r * rrotsgrid_m_[iq] - dr;
      deltaS_[offset] = ds;
      //spin integral norm is 1/(2.0 * pi), spatial integral has 4.0 * pi, so overall weight is 2
      spin_quad_weights_[offset] = 2.0 * spin_weight * sgridweight_m_[iq];
      count++;
    }
  };

  //simpsons 1/3 rule for spin integral

  //odd points
  for (int is = 1; is <= sknot_ - 1; is += 2)
  {
    RealType snew = smin + is * dS;
    addSpatialQuadrature(is, r, dr, snew - sold, RealType(4. / 3.) * dS);
  }

  //even points
  for (int is = 2; is <= sknot_ - 2; is += 2)
  {
    RealType snew = smin + is * dS;
    addSpatialQuadrature(is, r, dr, snew - sold, RealType(2. / 3.) * dS);
  }

  //end points
  addSpatialQuadrature(0, r, dr, smin - sold, RealType(1. / 3.) * dS);
  addSpatialQuadrature(sknot_, r, dr, smax - sold, RealType(1. / 3.) * dS);

  assert(count == total_knots_);

  //also set the radial function
  for (int ip = 0; ip < nchannel_; ip++)
    vrad_[ip] = sopp_m_[ip]->splint(r);
}

void SOECPComponent::rotateQuadratureGrid(const TensorType& rmat)
{
  for (int i = 0; i < sgridxyz_m_.size(); i++)
    rrotsgrid_m_[i] = dot(rmat, sgridxyz_m_[i]);
}

} // namespace qmcplusplus
