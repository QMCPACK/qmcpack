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


#ifndef QMCPLUSPLUS_BAREKINETICENERGY_H
#define QMCPLUSPLUS_BAREKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCDrivers/WalkerProperties.h"
#include "QMCWaveFunctions/SPOSet.h"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif
#include "type_traits/scalar_traits.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/** compute real(laplacian)
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<T, D>& g, T l)
{
  return dot(g, g) + l;
}

/** specialization of laplacian with complex g & l
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<std::complex<T>, D>& g, const std::complex<T>& l)
{
  return l.real() + OTCDot<T, T, D>::apply(g, g);
}

/** Convenience function to compute \f$\Re( \nabla^2_i \partial \Psi_T/\Psi_T)\f$
 * @param g OHMMS_DIM dimensional vector for \f$\nabla_i \ln \Psi_T\f$ .  
 * @param l A number, representing \f$\nabla^2_i \ln \Psi_T\f$ .
 * @param gg OHMMS_DIM dimensional vector containing \f$\nabla_i \partial \ln \Psi_T\f$ . 
 * @param gl A number, representing \f$\nabla^2_i \partial \ln \Psi_T\f$
 * @param ideriv A number, representing \f$\partial \ln \Psi_T\f$
 *
 * @return A number corresponding to \f$\Re( \nabla^2_i \partial \Psi_T/\Psi_T)\f$
 */

template<typename T, unsigned D>
inline T dlaplacian(const TinyVector<T, D>& g, const T l, const TinyVector<T, D>& gg, const T gl, const T ideriv)
{
  return gl + l * ideriv + 2.0 * dot(g, gg) + dot(g, g) * ideriv;
}

template<typename T, unsigned D>
inline T dlaplacian(const TinyVector<std::complex<T>, D>& g,
                    const std::complex<T> l,
                    const TinyVector<std::complex<T>, D>& gg,
                    const std::complex<T> gl,
                    const std::complex<T> ideriv)
{
  std::complex<T> l_times_ideriv                = l * ideriv;
  TinyVector<std::complex<T>, D> g_times_ideriv = g * ideriv;

  return gl.real() + l_times_ideriv.real() + 2.0 * OTCDot<T, T, D>::apply(g, gg) +
      OTCDot<T, T, D>::apply(g, g_times_ideriv);
}


/** @ingroup hamiltonian
  @brief Evaluate the kinetic energy with a single mass

 *The unit of the mass is AU, i.e., the electron mass \f$ m_e = 1 \f$.
 * To evaluate the Bare Kinetic part of the local energy
 \f$E_L({\bf R}) = \Psi^{-1}({\bf R})\hat{H}\Psi({\bf R}),\f$
 it is useful to use the following trick
 \f{eqnarray*}
 \nabla^2\Psi({\bf R}) &=& \nabla^2(\exp(\ln \Psi({\bf R})))\\
 &=&\nabla\cdot(\nabla\exp(\ln \Psi({\bf R}))) \\
 &=&\nabla\cdot(\nabla\ln \Psi({\bf R}))\exp(\ln \Psi({\bf R}))\\
 -\frac{1}{2}\frac{\nabla^2\Psi({\bf R})}{\Psi({\bf R})} &=&
 -\frac{1}{2}\nabla^2\ln \Psi({\bf R})
 -\frac{1}{2}(\nabla\ln \Psi({\bf R}))^2
 \f}
 */


inline std::string int2string(const int& i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}


template<typename T>
struct BareKineticEnergy : public OperatorBase
{
  using ValueMatrix_t = SPOSet::ValueMatrix_t;
  using GradMatrix_t  = SPOSet::GradMatrix_t;
  ///true, if all the species have the same mass
  bool SameMass;
  ///mass of the particle
  T M;
  ///\f$ 1/(2 m^*) \f$
  T OneOver2M;
  ///MinusOver2M[i] = \f$ -1/2m[i]\f$ for the ith species
  std::vector<T> MinusOver2M;

  ParticleSet::ParticleGradient_t Gtmp;
  ParticleSet::ParticleLaplacian_t Ltmp;

  ///single particle trace samples
  bool streaming_kinetic;
  bool streaming_kinetic_comp;
  bool streaming_momentum;

#if !defined(REMOVE_TRACEMANAGER)
  Array<TraceReal, 1>* T_sample;
  Array<TraceComp, 1>* T_sample_comp;
  Array<TraceComp, 2>* p_sample;
#endif
  ParticleSet& Ps;

  /** constructor
   *
   * Kinetic operators need to be re-evaluated during optimization.
   */
  BareKineticEnergy(RealType m = 1.0) : SameMass(true), M(m), OneOver2M(0.5 / m)
  {
    set_energy_domain(kinetic);
    set_quantum_domain(quantum);
    streaming_kinetic      = false;
    streaming_kinetic_comp = false;
    streaming_momentum     = false;
    UpdateMode.set(OPTIMIZABLE, 1);
    is_nondiag = true;
  }

  /** constructor with particleset
   * @param target particleset
   *
   * Store mass per species and use SameMass to choose the methods.
   * if SameMass, probably faster and easy to vectorize but no impact on the performance.
   */
  BareKineticEnergy(ParticleSet& p) : Ps(p)
  {
    set_energy_domain(kinetic);
    one_body_quantum_domain(p);
    streaming_kinetic      = false;
    streaming_kinetic_comp = false;
    streaming_momentum     = false;
    UpdateMode.set(OPTIMIZABLE, 1);
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
    is_nondiag = true;
  }
  ///destructor
  ~BareKineticEnergy() override {}

  void resetTargetParticleSet(ParticleSet& P) override {}


#if !defined(REMOVE_TRACEMANAGER)
  void contribute_particle_quantities() override
  {
    request.contribute_array(myName);
    request.contribute_array(myName + "_complex");
    request.contribute_array("momentum");
  }

  void checkout_particle_quantities(TraceManager& tm) override
  {
    streaming_particles = request.streaming_array(myName) || request.streaming_array(myName + "_complex") ||
        request.streaming_array("momentum");
    if (streaming_particles)
    {
      T_sample      = tm.checkout_real<1>(myName, Ps);
      T_sample_comp = tm.checkout_complex<1>(myName + "_complex", Ps);
      p_sample      = tm.checkout_complex<2>("momentum", Ps, DIM);
    }
  }

  void delete_particle_quantities() override
  {
    if (streaming_particles)
    {
      delete T_sample;
      delete T_sample_comp;
      delete p_sample;
    }
  }
#endif


  inline Return_t evaluate(ParticleSet& P) override
  {
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles)
    {
      Value = evaluate_sp(P);
    }
    else
#endif
        if (SameMass)
    {
//app_log() << "Here" << std::endl;
#ifdef QMC_COMPLEX
      Value = std::real(CplxDot(P.G, P.G) + CplxSum(P.L));
      Value *= -OneOver2M;
#else
      Value = Dot(P.G, P.G) + Sum(P.L);
      Value *= -OneOver2M;
#endif
    }
    else
    {
      Value = 0.0;
      for (int i = 0; i < MinusOver2M.size(); ++i)
      {
        T x = 0.0;
        for (int j = P.first(i); j < P.last(i); ++j)
          x += laplacian(P.G[j], P.L[j]);
        Value += x * MinusOver2M[i];
      }
    }
    return Value;
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
  inline Return_t evaluateWithIonDerivs(ParticleSet& P,
                                        ParticleSet& ions,
                                        TrialWaveFunction& psi,
                                        ParticleSet::ParticlePos_t& hf_terms,
                                        ParticleSet::ParticlePos_t& pulay_terms) override
  {
    typedef ParticleSet::ParticlePos_t ParticlePos_t;
    typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
    typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;

    int Nions = ions.getTotalNum();
    int Nelec = P.getTotalNum();

    //These are intermediate arrays for potentially complex math.
    ParticleLaplacian_t term2_(Nelec);
    ParticleGradient_t term4_(Nelec);

    //Potentially complex temporary array for \partial \psi/\psi and \nabla^2 \partial \psi / \psi
    ParticleGradient_t iongradpsi_(Nions), pulaytmp_(Nions);
    //temporary arrays that will be explicitly real.
    ParticlePos_t pulaytmpreal_(Nions), iongradpsireal_(Nions);


    TinyVector<ParticleGradient_t, OHMMS_DIM> iongrad_grad_;
    TinyVector<ParticleLaplacian_t, OHMMS_DIM> iongrad_lapl_;

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
      convert(iongradpsi_[iat], iongradpsireal_[iat]);
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
      convert(pulaytmp_[iat], pulaytmpreal_[iat]);
    }

    if (SameMass)
    {
      Value = Dot(P.G, P.G) + Sum(P.L);
      Value *= -OneOver2M;
    }
    else
    {
      Value = 0.0;
      for (int i = 0; i < MinusOver2M.size(); ++i)
      {
        T x = 0.0;
        for (int j = P.first(i); j < P.last(i); ++j)
          x += laplacian(P.G[j], P.L[j]);
        Value += x * MinusOver2M[i];
      }
    }
    pulaytmpreal_ -= Value * iongradpsireal_;


    pulay_terms += pulaytmpreal_;
    return Value;
  }


#if !defined(REMOVE_TRACEMANAGER)
  inline Return_t evaluate_sp(ParticleSet& P)
  {
    Array<RealType, 1>& T_samp                    = *T_sample;
    Array<std::complex<RealType>, 1>& T_samp_comp = *T_sample_comp;
    Array<std::complex<RealType>, 2>& p_samp      = *p_sample;
    std::complex<RealType> t1                     = 0.0;
    const RealType clambda(-OneOver2M);
    Value = 0.0;
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
        Value += real(t1);
      }
    }
    else
    {
      for (int s = 0; s < MinusOver2M.size(); ++s)
      {
        T mlambda = MinusOver2M[s];
        for (int i = P.first(s); i < P.last(s); ++i)
        {
          //t1 = mlambda*( P.L[i] + dot(P.G[i],P.G[i]) );
          t1 = P.L[i] + dot(P.G[i], P.G[i]);
          t1 *= mlambda;
          T_samp(i)      = real(t1);
          T_samp_comp(i) = t1;
          for (int d = 0; d < DIM; ++d)
            p_samp(i, d) = P.G[i][d];
          Value += real(t1);
        }
      }
    }
#if defined(TRACE_CHECK)
    RealType Vnow = Value;
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
    return Value;
  }

#endif

  inline Return_t evaluate_orig(ParticleSet& P)
  {
    if (SameMass)
    {
      Value = Dot(P.G, P.G) + Sum(P.L);
      Value *= -OneOver2M;
    }
    else
    {
      Value = 0.0;
      for (int i = 0; i < MinusOver2M.size(); ++i)
      {
        T x = 0.0;
        for (int j = P.first(i); j < P.last(i); ++j)
          x += laplacian(P.G[j], P.L[j]);
        Value += x * MinusOver2M[i];
      }
    }
    return Value;
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */
  bool put(xmlNodePtr) override { return true; }

  bool get(std::ostream& os) const override
  {
    os << "Kinetic energy";
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final
  {
    return std::make_unique<BareKineticEnergy>(*this);
  }

#ifdef QMC_CUDA
  ////////////////////////////////
  // Vectorized version for GPU //
  ////////////////////////////////
  // Nothing is done on GPU here, just copy into vector
  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override
  {
    auto& walkers = W.WalkerList;
    for (int iw = 0; iw < walkers.size(); iw++)
    {
      Walker_t& w                                      = *(walkers[iw]);
      RealType KE                                      = -OneOver2M * (Dot(w.G, w.G) + Sum(w.L));
      w.getPropertyBase()[WP::NUMPROPERTIES + myIndex] = KE;
      LocalEnergy[iw] += KE;
    }
  }
#endif

  void evaluateOneBodyOpMatrix(ParticleSet& P, TWFPrototype& psi, std::vector<ValueMatrix_t>& B) override
  {
    IndexType ngroups = P.groups();
    assert(B.size() == ngroups);
    std::vector<ValueMatrix_t> M;
    std::vector<GradMatrix_t> grad_M;
    std::vector<ValueMatrix_t> lapl_M;
    for (int ig = 0; ig < ngroups; ig++)
    {
      IndexType norbs    = psi.num_orbitals(ig);
      IndexType numptcls = psi.num_particles(ig);

      ValueMatrix_t zeromat;
      GradMatrix_t zerogradmat;

      zeromat.resize(numptcls, norbs);
      zerogradmat.resize(numptcls, norbs);

      M.push_back(zeromat);
      grad_M.push_back(zerogradmat);
      lapl_M.push_back(zeromat);
    }

    psi.get_egrad_elapl_M(P, M, grad_M, lapl_M);
    for (int ig = 0; ig < ngroups; ig++)
    {
      lapl_M[ig] *= MinusOver2M[ig];
      B[ig] += lapl_M[ig];
    }
  }

  //Bforce.  First index is the x, y, or z component.  so Bforce[0] would be the std::vector<ValueMatrix_t> of d/dx B.
  void evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                         ParticleSet& source,
                                         TWFPrototype& psi,
                                         int iat,
                                         std::vector<std::vector<ValueMatrix_t>>& Bforce) override
  {
    ScopedTimer dBketimer(*timer_manager.createTimer("NEW::KE::dB"));
    IndexType ngroups = P.groups();
    assert(Bforce.size() == OHMMS_DIM);
    assert(Bforce[0].size() == ngroups);
    std::vector<ValueMatrix_t> mtmp;
    for (int ig = 0; ig < ngroups; ig++)
    {
      IndexType norbs    = psi.num_orbitals(ig);
      IndexType numptcls = psi.num_particles(ig);

      ValueMatrix_t zeromat;
      GradMatrix_t zerogradmat;

      zeromat.resize(numptcls, norbs);
      zerogradmat.resize(numptcls, norbs);

      mtmp.push_back(zeromat);
    }

    std::vector<std::vector<ValueMatrix_t>> dm, dlapl;
    dm.push_back(mtmp);
    dm.push_back(mtmp);
    dm.push_back(mtmp);

    dlapl.push_back(mtmp);
    dlapl.push_back(mtmp);
    dlapl.push_back(mtmp);

    {
      ScopedTimer gradorbtimer(*timer_manager.createTimer("NEW::KE::orb_grad"));
      psi.get_igrad_igradelapl_M(P, source, iat, dm, dlapl);
    }
    for (int idim = 0; idim < OHMMS_DIM; idim++)
      for (int ig = 0; ig < ngroups; ig++)
      {
        dlapl[idim][ig] *= MinusOver2M[ig];
        Bforce[idim][ig] += dlapl[idim][ig];
      }
  }
};
} // namespace qmcplusplus
#endif
