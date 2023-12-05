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
#include "type_traits/ConvertToReal.h"

namespace qmcplusplus
{
using WP       = WalkerProperties::Indexes;
using Return_t = BareKineticEnergy::Return_t;

struct BareKineticEnergy::MultiWalkerResource : public Resource
{
  MultiWalkerResource() : Resource("BareKineticEnergy") {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<MultiWalkerResource>(*this); }

  Vector<RealType> t_samples;
  Vector<std::complex<RealType>> tcmp_samples;
};

/** constructor with particleset
   * @param target particleset
   *
   * Store mass per species and use SameMass to choose the methods.
   * if SameMass, probably faster and easy to vectorize but no impact on the performance.
   */
BareKineticEnergy::BareKineticEnergy(ParticleSet& p, TrialWaveFunction& psi) : ps_(p), psi_(psi)
{
  setEnergyDomain(KINETIC);
  oneBodyQuantumDomain(p);
  update_mode_.set(OPTIMIZABLE, 1);
  SpeciesSet& tspecies(p.getSpeciesSet());
  minus_over_2m_.resize(tspecies.size());
  int massind    = tspecies.addAttribute("mass");
  same_mass_     = true;
  particle_mass_ = tspecies(massind, 0);
  one_over_2m_   = 0.5 / particle_mass_;
  for (int i = 0; i < tspecies.size(); ++i)
  {
    same_mass_ &= (std::abs(tspecies(massind, i) - particle_mass_) < 1e-6);
    minus_over_2m_[i] = -1.0 / (2.0 * tspecies(massind, i));
  }
}

///destructor
BareKineticEnergy::~BareKineticEnergy() = default;

bool BareKineticEnergy::dependsOnWaveFunction() const { return true; }

std::string BareKineticEnergy::getClassName() const { return "BareKineticEnergy"; }

void BareKineticEnergy::resetTargetParticleSet(ParticleSet& p) {}

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
    t_sample_      = tm.checkout_real<1>(name_, ps_);
    t_sample_comp_ = tm.checkout_complex<1>(name_ + "_complex", ps_);
    p_sample_      = tm.checkout_complex<2>("momentum", ps_, DIM);
  }
}

void BareKineticEnergy::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete t_sample_;
    delete t_sample_comp_;
    delete p_sample_;
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
      if (same_mass_)
  {
#ifdef QMC_COMPLEX
    value_ = std::real(CplxDot(P.G, P.G) + CplxSum(P.L));
    value_ *= -one_over_2m_;
#else
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -one_over_2m_;
#endif
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < minus_over_2m_.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * minus_over_2m_[i];
    }
  }
  return value_;
}

Return_t BareKineticEnergy::evaluateValueAndDerivatives(ParticleSet& P,
                                                        const opt_variables_type& optvars,
                                                        const Vector<ValueType>& dlogpsi,
                                                        Vector<ValueType>& dhpsioverpsi)
{
  // const_cast is needed because TWF::evaluateDerivatives calculates dlogpsi.
  // KineticEnergy must be the first element in the hamiltonian array.
  psi_.evaluateDerivatives(P, optvars, const_cast<Vector<ValueType>&>(dlogpsi), dhpsioverpsi);
  return evaluate(P);
}

void BareKineticEnergy::mw_evaluateWithParameterDerivatives(const RefVectorWithLeader<OperatorBase>& o_list,
                                                            const RefVectorWithLeader<ParticleSet>& p_list,
                                                            const opt_variables_type& optvars,
                                                            const RecordArray<ValueType>& dlogpsi,
                                                            RecordArray<ValueType>& dhpsioverpsi) const
{
  RefVectorWithLeader<TrialWaveFunction> wf_list(o_list.getCastedLeader<BareKineticEnergy>().psi_);
  for (int i = 0; i < o_list.size(); i++)
    wf_list.push_back(o_list.getCastedElement<BareKineticEnergy>(i).psi_);
  mw_evaluate(o_list, wf_list, p_list);
  // const_cast is needed because TWF::evaluateDerivatives calculates dlogpsi.
  // KineticEnergy must be the first element in the hamiltonian array.
  TrialWaveFunction::mw_evaluateParameterDerivatives(wf_list, p_list, optvars,
                                                     const_cast<RecordArray<ValueType>&>(dlogpsi), dhpsioverpsi);
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
    if (same_mass_)
    {
      for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
      {
        //These term[24]_ variables exist because I want to do complex math first, and then take the real part at the
        //end.  Sum() and Dot() perform the needed functions and spit out the real part at the end.
        term2_                 = P.L * iongradpsi_[iat][iondim];
        term4_                 = P.G * iongradpsi_[iat][iondim];
        pulaytmp_[iat][iondim] = -one_over_2m_ *
            (Sum(iongrad_lapl_[iondim]) + Sum(term2_) + 2.0 * Dot(iongrad_grad_[iondim], P.G) + Dot(P.G, term4_));
      }
    }
    else
    {
      for (int iondim = 0; iondim < OHMMS_DIM; iondim++)
      {
        for (int g = 0; g < minus_over_2m_.size(); g++)
        {
          for (int iel = P.first(g); iel < P.last(g); iel++)
          {
            pulaytmp_[iat][iondim] += minus_over_2m_[g] *
                (dlaplacian(P.G[iel], P.L[iel], iongrad_grad_[iondim][iel], iongrad_lapl_[iondim][iel],
                            iongradpsi_[iat][iondim]));
          }
        }
      }
    }
    convertToReal(pulaytmp_[iat], pulaytmpreal_[iat]);
  }

  if (same_mass_)
  {
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -one_over_2m_;
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < minus_over_2m_.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * minus_over_2m_[i];
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
        B[sid][iel - first][iorb] += RealType(minus_over_2m_[ig]) *
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
          Bforce[idim][sid][iel - first][iorb] = RealType(minus_over_2m_[ig]) *
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
  auto& O_leader   = o_list.getCastedLeader<BareKineticEnergy>();
  O_leader.mw_res_ = collection.lendResource<MultiWalkerResource>();
}

void BareKineticEnergy::releaseResource(ResourceCollection& collection,
                                        const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<BareKineticEnergy>();
  collection.takebackResource(O_leader.mw_res_);
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
  Vector<RealType>& t_samp                  = o_leader.mw_res_.getResource().t_samples;
  Vector<std::complex<RealType>>& tcmp_samp = o_leader.mw_res_.getResource().tcmp_samples;

  auto num_species = p_leader.getSpeciesSet().getTotalNum();
  t_samp.resize(num_particles);
  tcmp_samp.resize(num_particles);
  const RealType clambda(-one_over_2m_);
  auto evaluate_walker_per_particle = [num_particles, name, &t_samp, &tcmp_samp,
                                       clambda](const int walker_index, const BareKineticEnergy& bke,
                                                const ParticleSet& pset,
                                                const std::vector<ListenerVector<RealType>>& listeners) {
    RealType value = 0;
    std::fill(t_samp.begin(), t_samp.end(), 0.0);
    std::fill(tcmp_samp.begin(), tcmp_samp.end(), 0.0);

    std::complex<RealType> t1 = 0.0;
    if (bke.same_mass_)
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
      for (int s = 0; s < bke.minus_over_2m_.size(); ++s)
      {
        FullPrecRealType mlambda = bke.minus_over_2m_[s];
        for (int i = pset.first(s); i < pset.last(s); ++i)
        {
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
  Array<RealType, 1>& T_samp                    = *t_sample_;
  Array<std::complex<RealType>, 1>& T_samp_comp = *t_sample_comp_;
  Array<std::complex<RealType>, 2>& p_samp      = *p_sample_;
  std::complex<RealType> t1                     = 0.0;
  const RealType clambda(-one_over_2m_);
  value_ = 0.0;
  if (same_mass_)
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
    for (int s = 0; s < minus_over_2m_.size(); ++s)
    {
      FullPrecRealType mlambda = minus_over_2m_[s];
      for (int i = P.first(s); i < P.last(s); ++i)
      {
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
  if (same_mass_)
  {
    value_ = Dot(P.G, P.G) + Sum(P.L);
    value_ *= -one_over_2m_;
  }
  else
  {
    value_ = 0.0;
    for (int i = 0; i < minus_over_2m_.size(); ++i)
    {
      Return_t x = 0.0;
      for (int j = P.first(i); j < P.last(i); ++j)
        x += laplacian(P.G[j], P.L[j]);
      value_ += x * minus_over_2m_[i];
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
  return std::make_unique<BareKineticEnergy>(qp, psi);
}

} // namespace qmcplusplus
