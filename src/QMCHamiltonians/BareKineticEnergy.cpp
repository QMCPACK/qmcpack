//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "BareKineticEnergy.h"
#include "Particle/ParticleSet.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "BareKineticHelper.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
using WP       = WalkerProperties::Indexes;
using Return_t = BareKineticEnergy::Return_t;

/** constructor with particleset
   * @param target particleset
   *
   * Store mass per species and use SameMass to choose the methods.
   * if SameMass, probably faster and easy to vectorize but no impact on the performance.
   */
BareKineticEnergy::BareKineticEnergy(ParticleSet& p) : Ps(p)
{
  setEnergyDomain(KINETIC);
  oneBodyQuantumDomain(p);
  streaming_kinetic      = false;
  streaming_kinetic_comp = false;
  streaming_momentum     = false;
  update_mode_.set(OPTIMIZABLE, 1);
  SpeciesSet& tspecies(p.getSpeciesSet());
  MinusOver2M.resize(tspecies.size());
  int massind = tspecies.addAttribute("mass");
  SameMass    = true;
  M           = tspecies(massind, 0);
  OneOver2M   = 0.5 / M;
  for (int i = 0; i < tspecies.size(); ++i)
  {
    SameMass &= (std::abs(tspecies(massind, i) - M) < 1e-6);
    MinusOver2M[i] = -1.0 / (2.0 * tspecies(massind, i));
  }
}

///destructor
BareKineticEnergy::~BareKineticEnergy() = default;

#if !defined(REMOVE_TRACEMANAGER)
void BareKineticEnergy::contributeParticleQuantities()
{
  request_.contribute_array(name_);
  request_.contribute_array(name_ + "_complex");
  request_.contribute_array("momentum");
}

void BareKineticEnergy::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_) || request_.streaming_array(name_ + "_complex") ||
      request_.streaming_array("momentum");
  if (streaming_particles_)
  {
    T_sample      = tm.checkout_real<1>(name_, Ps);
    T_sample_comp = tm.checkout_complex<1>(name_ + "_complex", Ps);
    p_sample      = tm.checkout_complex<2>("momentum", Ps, DIM);
  }
}

void BareKineticEnergy::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete T_sample;
    delete T_sample_comp;
    delete p_sample;
  }
}
#endif


Return_t BareKineticEnergy::evaluate(ParticleSet& P)
{
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles_)
  {
    value_ = evaluate_sp(P);
  }
  else
#endif
      if (SameMass)
  {
//app_log() << "Here" << std::endl;
#ifdef QMC_COMPLEX
    value_ = std::real(CplxDot(P.G, P.G) + CplxSum(P.L));
    value_ *= -OneOver2M;
#else
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -OneOver2M;
#endif
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < MinusOver2M.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * MinusOver2M[i];
    }
  }
  return value_;
}

/**@brief Function to compute the value, direct ionic gradient terms, and pulay terms for the local kinetic energy.
 *  
 *  This general function represents the OperatorBase interface for computing.  For an operator \hat{O}, this
 *  function will return \frac{\hat{O}\Psi_T}{\Psi_T},  \frac{\partial(\hat{O})\Psi_T}{\Psi_T}, and 
 *  \frac{\hat{O}\partial\Psi_T}{\Psi_T} - \frac{\hat{O}\Psi_T}{\Psi_T}\frac{\partial \Psi_T}{\Psi_T}.  These are 
 *  referred to as Value, HF term, and pulay term respectively.
 *
 * @param P electron particle set.
 * @param ions ion particle set
 * @param psi Trial wave function object.
 * @param hf_terms 3Nion dimensional object. All direct force terms, or ionic gradient of operator itself.
 *                 Contribution of this operator is ADDED onto hf_terms.
 * @param pulay_terms The terms coming from ionic gradients of trial wavefunction.  Contribution of this operator is
 *                 ADDED onto pulay_terms.
 * @return Value of kinetic energy operator at electron/ion positions given by P and ions.  The force contributions from
 *          this operator are added into hf_terms and pulay_terms.
 */
Return_t BareKineticEnergy::evaluateWithIonDerivs(ParticleSet& P,
                                                  ParticleSet& ions,
                                                  TrialWaveFunction& psi,
                                                  ParticleSet::ParticlePos& hf_terms,
                                                  ParticleSet::ParticlePos& pulay_terms)
{
  using ParticlePos       = ParticleSet::ParticlePos;
  using ParticleGradient  = ParticleSet::ParticleGradient;
  using ParticleLaplacian = ParticleSet::ParticleLaplacian;

  int Nions = ions.getTotalNum();
  int Nelec = P.getTotalNum();

  //These are intermediate arrays for potentially complex math.
  ParticleLaplacian term2_(Nelec);
  ParticleGradient term4_(Nelec);

  //Potentially complex temporary array for \partial \psi/\psi and \nabla^2 \partial \psi / \psi
  ParticleGradient iongradpsi_(Nions), pulaytmp_(Nions);
  //temporary arrays that will be explicitly real.
  ParticlePos pulaytmpreal_(Nions), iongradpsireal_(Nions);


  TinyVector<ParticleGradient, OHMMS_DIM> iongrad_grad_;
  TinyVector<ParticleLaplacian, OHMMS_DIM> iongrad_lapl_;

  for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
  {
    iongrad_grad_[iondim].resize(Nelec);
    iongrad_lapl_[iondim].resize(Nelec);
  }

  iongradpsi_     = 0;
  iongradpsireal_ = 0;
  pulaytmpreal_   = 0;
  pulaytmp_       = 0;

  RealType logpsi_ = psi.evaluateLog(P);
  for (int iat = 0; iat < Nions; iat++)
  {
    //reset the iongrad_X containers.
    for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
    {
      iongrad_grad_[iondim] = 0;
      iongrad_lapl_[iondim] = 0;
    }
    iongradpsi_[iat] = psi.evalGradSource(P, ions, iat, iongrad_grad_, iongrad_lapl_);
    //conversion from potentially complex to definitely real.
    convertToReal(iongradpsi_[iat], iongradpsireal_[iat]);
    if (SameMass)
    {
      for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
      {
        //These term[24]_ variables exist because I want to do complex math first, and then take the real part at the
        //end.  Sum() and Dot() perform the needed functions and spit out the real part at the end.
        term2_                 = P.L * iongradpsi_[iat][iondim];
        term4_                 = P.G * iongradpsi_[iat][iondim];
        pulaytmp_[iat][iondim] = -OneOver2M *
            (Sum(iongrad_lapl_[iondim]) + Sum(term2_) + 2.0 * Dot(iongrad_grad_[iondim], P.G) + Dot(P.G, term4_));
      }
    }
    else
    {
      for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
      {
        for (int g = 0; g < MinusOver2M.size(); g++)
        {
          for (int iel = P.first(g); iel < P.last(g); iel++)
          {
            pulaytmp_[iat][iondim] += MinusOver2M[g] *
                (dlaplacian(P.G[iel], P.L[iel], iongrad_grad_[iondim][iel], iongrad_lapl_[iondim][iel],
                            iongradpsi_[iat][iondim]));
          }
        }
      }
    }
    //convert to real.
    convertToReal(pulaytmp_[iat], pulaytmpreal_[iat]);
  }

  if (SameMass)
  {
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -OneOver2M;
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < MinusOver2M.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * MinusOver2M[i];
    }
  }
  pulaytmpreal_ -= value_ * iongradpsireal_;


  pulay_terms += pulaytmpreal_;
  return value_;
}

void BareKineticEnergy::evaluateOneBodyOpMatrix(ParticleSet& P,
                                                const TWFFastDerivWrapper& psi,
                                                std::vector<ValueMatrix>& B)
{
  ParticleSet::ParticleGradient G;
  ParticleSet::ParticleLaplacian L;

  const IndexType nelec = P.getTotalNum();
  G.resize(nelec);
  L.resize(nelec);

  const IndexType ngroups = P.groups();
  assert(B.size() == ngroups);
  std::vector<ValueMatrix> M;
  std::vector<GradMatrix> grad_M;
  std::vector<ValueMatrix> lapl_M;
  std::vector<ValueMatrix> gradJdotgradPhi;
  for (int ig = 0; ig < ngroups; ig++)
  {
    const IndexType sid    = psi.getTWFGroupIndex(ig);
    const IndexType norbs  = psi.numOrbitals(sid);
    const IndexType first  = P.first(ig);
    const IndexType last   = P.last(ig);
    const IndexType nptcls = last - first;
    ValueMatrix zeromat;
    GradMatrix zerogradmat;

    zeromat.resize(nptcls, norbs);
    zerogradmat.resize(nptcls, norbs);

    M.push_back(zeromat);
    grad_M.push_back(zerogradmat);
    lapl_M.push_back(zeromat);
    gradJdotgradPhi.push_back(zeromat);
  }

  psi.getEGradELaplM(P, M, grad_M, lapl_M);
  psi.evaluateJastrowVGL(P, G, L);

  for (int ig = 0; ig < ngroups; ig++)
  {
    const IndexType sid    = psi.getTWFGroupIndex(ig);
    const IndexType norbs  = psi.numOrbitals(sid);
    const IndexType first  = P.first(ig);
    const IndexType last   = P.last(ig);
    const IndexType nptcls = last - first;
    for (int iel = first; iel < last; iel++)
    {
      for (int iorb = 0; iorb < norbs; iorb++)
      {
        gradJdotgradPhi[sid][iel - first][iorb] = RealType(2.0) * dot(GradType(G[iel]), grad_M[sid][iel - first][iorb]);
        B[sid][iel - first][iorb] += RealType(MinusOver2M[ig]) *
            (lapl_M[sid][iel - first][iorb] + gradJdotgradPhi[sid][iel - first][iorb] +
             ValueType(L[iel] + dot(G[iel], G[iel])) * M[sid][iel - first][iorb]);
      }
    }
  }
}

void BareKineticEnergy::evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                                          ParticleSet& source,
                                                          const TWFFastDerivWrapper& psi,
                                                          const int iat,
                                                          std::vector<std::vector<ValueMatrix>>& Bforce)
{
  const IndexType ngroups = P.groups();
  const IndexType nelec   = P.getTotalNum();

  ParticleSet::ParticleGradient Gtmp, G;
  ParticleSet::ParticleLaplacian Ltmp, L;
  Gtmp.resize(nelec);
  G.resize(nelec);
  Ltmp.resize(nelec);
  L.resize(nelec);

  std::vector<ValueMatrix> M;
  std::vector<GradMatrix> grad_M;
  std::vector<ValueMatrix> lapl_M;

  TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM> dG;
  TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM> dL;

  for (int dim = 0; dim < OHMMS_DIM; dim++)
  {
    dG[dim] = Gtmp;
    dL[dim] = Ltmp;
  }

  assert(Bforce.size() == OHMMS_DIM);
  assert(Bforce[0].size() == ngroups);
  std::vector<ValueMatrix> mtmp;
  for (int ig = 0; ig < ngroups; ig++)
  {
    const IndexType sid    = psi.getTWFGroupIndex(ig);
    const IndexType norbs  = psi.numOrbitals(sid);
    const IndexType first  = P.first(ig);
    const IndexType last   = P.last(ig);
    const IndexType nptcls = last - first;

    ValueMatrix zeromat;
    GradMatrix zerogradmat;

    zeromat.resize(nptcls, norbs);
    zerogradmat.resize(nptcls, norbs);

    mtmp.push_back(zeromat);
    M.push_back(zeromat);
    grad_M.push_back(zerogradmat);
    lapl_M.push_back(zeromat);
  }


  std::vector<std::vector<ValueMatrix>> dm, dlapl;
  std::vector<std::vector<GradMatrix>> dgmat;
  dm.push_back(mtmp);
  dm.push_back(mtmp);
  dm.push_back(mtmp);

  dlapl.push_back(mtmp);
  dlapl.push_back(mtmp);
  dlapl.push_back(mtmp);

  dgmat.push_back(grad_M);
  dgmat.push_back(grad_M);
  dgmat.push_back(grad_M);

  psi.getEGradELaplM(P, M, grad_M, lapl_M);
  psi.getIonGradIonGradELaplM(P, source, iat, dm, dgmat, dlapl);
  psi.evaluateJastrowVGL(P, G, L);
  psi.evaluateJastrowGradSource(P, source, iat, dG, dL);
  for (int idim = 0; idim < OHMMS_DIM; idim++)
    for (int ig = 0; ig < ngroups; ig++)
    {
      const IndexType sid    = psi.getTWFGroupIndex(ig);
      const IndexType norbs  = psi.numOrbitals(sid);
      const IndexType first  = P.first(ig);
      const IndexType last   = P.last(ig);
      const IndexType nptcls = last - first;

      for (int iel = first; iel < last; iel++)
      {
        for (int iorb = 0; iorb < norbs; iorb++)
        {
          Bforce[idim][sid][iel - first][iorb] = RealType(MinusOver2M[ig]) *
              (dlapl[idim][sid][iel - first][iorb] +
               RealType(2.0) *
                   (dot(GradType(G[iel]), dgmat[idim][sid][iel - first][iorb]) +
                    dot(GradType(dG[idim][iel]), grad_M[sid][iel - first][iorb])) +
               M[sid][iel - first][iorb] * ValueType(dL[idim][iel] + 2.0 * dot(dG[idim][iel], G[iel])) +
               ValueType(L[iel] + dot(G[iel], G[iel])) * dm[idim][sid][iel - first][iorb]);
        }
      }
    }
}

void BareKineticEnergy::createResource(ResourceCollection& collection) const
{
  auto new_res        = std::make_unique<MultiWalkerResource>();
  auto resource_index = collection.addResource(std::move(new_res));
}

void BareKineticEnergy::acquireResource(ResourceCollection& collection,
                                        const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<BareKineticEnergy>();
  auto res_ptr   = dynamic_cast<MultiWalkerResource*>(collection.lendResource().release());
  if (!res_ptr)
    throw std::runtime_error("BareKineticEnergy::acquireResource dynamic_cast failed");
  O_leader.mw_res_.reset(res_ptr);
}

void BareKineticEnergy::releaseResource(ResourceCollection& collection,
                                        const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<BareKineticEnergy>();
  collection.takebackResource(std::move(O_leader.mw_res_));
}

void BareKineticEnergy::mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                               const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                               const RefVectorWithLeader<ParticleSet>& p_list,
                                               const std::vector<ListenerVector<RealType>>& listeners,
                                               const std::vector<ListenerVector<RealType>>& ion_listeners) const
{
  auto& o_leader = o_list.getCastedLeader<BareKineticEnergy>();
  auto& p_leader = p_list.getLeader();
  assert(this == &o_list.getLeader());

  auto num_particles                        = p_leader.getTotalNum();
  auto& name                                = o_leader.name_;
  Vector<RealType>& t_samp                  = o_leader.mw_res_->t_samples;
  Vector<std::complex<RealType>>& tcmp_samp = o_leader.mw_res_->tcmp_samples;

  auto num_species = p_leader.getSpeciesSet().getTotalNum();
  t_samp.resize(num_particles);
  tcmp_samp.resize(num_particles);
  const RealType clambda(-OneOver2M);
  auto evaluate_walker_per_particle = [num_particles, name, &t_samp, &tcmp_samp,
                                       clambda](const int walker_index, const BareKineticEnergy& bke,
                                                const ParticleSet& pset,
                                                const std::vector<ListenerVector<RealType>>& listeners) {
    RealType value = 0;
    std::fill(t_samp.begin(), t_samp.end(), 0.0);
    std::fill(tcmp_samp.begin(), tcmp_samp.end(), 0.0);

    std::complex<RealType> t1 = 0.0;
    if (bke.SameMass)
      for (int i = 0; i < num_particles; i++)
      {
        t1           = pset.L[i] + dot(pset.G[i], pset.G[i]);
        t1           = clambda * t1;
        t_samp[i]    = real(t1);
        tcmp_samp[i] = t1;
        value += real(t1);
      }
    else
    {
      for (int s = 0; s < bke.MinusOver2M.size(); ++s)
      {
        FullPrecRealType mlambda = bke.MinusOver2M[s];
        for (int i = pset.first(s); i < pset.last(s); ++i)
        {
          //t1 = mlambda*( pset.L[i] + dot(pset.G[i],pset.G[i]) );
          t1 = pset.L[i] + dot(pset.G[i], pset.G[i]);
          t1 *= mlambda;
          t_samp[i]    = real(t1);
          tcmp_samp[i] = t1;
          value += real(t1);
        }
      }
    }
    for (auto& listener : listeners)
    {
      listener.report(walker_index, name, t_samp);
    }
    return value;
  };

  for (int iw = 0; iw < o_list.size(); iw++)
  {
    auto& bare_kinetic_energy  = o_list.getCastedElement<BareKineticEnergy>(iw);
    bare_kinetic_energy.value_ = evaluate_walker_per_particle(iw, bare_kinetic_energy, p_list[iw], listeners);
  }
}

void BareKineticEnergy::mw_evaluatePerParticleWithToperator(
    const RefVectorWithLeader<OperatorBase>& o_list,
    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
    const RefVectorWithLeader<ParticleSet>& p_list,
    const std::vector<ListenerVector<RealType>>& listeners,
    const std::vector<ListenerVector<RealType>>& ion_listeners) const
{
  mw_evaluatePerParticle(o_list, wf_list, p_list, listeners, ion_listeners);
}

#if !defined(REMOVE_TRACEMANAGER)
Return_t BareKineticEnergy::evaluate_sp(ParticleSet& P)
{
  Array<RealType, 1>& T_samp                    = *T_sample;
  Array<std::complex<RealType>, 1>& T_samp_comp = *T_sample_comp;
  Array<std::complex<RealType>, 2>& p_samp      = *p_sample;
  std::complex<RealType> t1                     = 0.0;
  const RealType clambda(-OneOver2M);
  value_ = 0.0;
  if (SameMass)
  {
    for (int i = 0; i < P.getTotalNum(); i++)
    {
      t1             = P.L[i] + dot(P.G[i], P.G[i]);
      t1             = clambda * t1;
      T_samp(i)      = real(t1);
      T_samp_comp(i) = t1;
      for (int d = 0; d < DIM; ++d)
        p_samp(i, d) = P.G[i][d];
      value_ += real(t1);
    }
  }
  else
  {
    for (int s = 0; s < MinusOver2M.size(); ++s)
    {
      FullPrecRealType mlambda = MinusOver2M[s];
      for (int i = P.first(s); i < P.last(s); ++i)
      {
        //t1 = mlambda*( P.L[i] + dot(P.G[i],P.G[i]) );
        t1 = P.L[i] + dot(P.G[i], P.G[i]);
        t1 *= mlambda;
        T_samp(i)      = real(t1);
        T_samp_comp(i) = t1;
        for (int d = 0; d < DIM; ++d)
          p_samp(i, d) = P.G[i][d];
        value_ += real(t1);
      }
    }
  }
#if defined(TRACE_CHECK)
  RealType Vnow = value_;
  RealType Vsum = T_samp.sum();
  RealType Vold = evaluate_orig(P);
  if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
  {
    app_log() << "accumtest: BareKineticEnergy::evaluate()" << std::endl;
    app_log() << "accumtest:   tot:" << Vnow << std::endl;
    app_log() << "accumtest:   sum:" << Vsum << std::endl;
    APP_ABORT("Trace check failed");
  }
  if (std::abs(Vold - Vnow) > TraceManager::trace_tol)
  {
    app_log() << "versiontest: BareKineticEnergy::evaluate()" << std::endl;
    app_log() << "versiontest:   orig:" << Vold << std::endl;
    app_log() << "versiontest:    mod:" << Vnow << std::endl;
    APP_ABORT("Trace check failed");
  }
#endif
  return value_;
}

#endif

Return_t BareKineticEnergy::evaluate_orig(ParticleSet& P)
{
  if (SameMass)
  {
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -OneOver2M;
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < MinusOver2M.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * MinusOver2M[i];
    }
  }
  return value_;
}

/** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */
bool BareKineticEnergy::put(xmlNodePtr) { return true; }

bool BareKineticEnergy::get(std::ostream& os) const
{
  os << "Kinetic energy";
  return true;
}

std::unique_ptr<OperatorBase> BareKineticEnergy::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<BareKineticEnergy>(*this);
}

#ifdef QMC_CUDA
////////////////////////////////
// Vectorized version for GPU //
////////////////////////////////
// Nothing is done on GPU here, just copy into vector
void BareKineticEnergy::addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
{
  auto& walkers = W.WalkerList;
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& w                                        = *(walkers[iw]);
    RealType KE                                        = -OneOver2M * (Dot(w.G, w.G) + Sum(w.L));
    w.getPropertyBase()[WP::NUMPROPERTIES + my_index_] = KE;
    LocalEnergy[iw] += KE;
  }
}
#endif
} // namespace qmcplusplus
