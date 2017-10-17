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
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif
#include "type_traits/scalar_traits.h"

namespace qmcplusplus
{

/** compute real(laplacian)
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<T,D>& g, T l)
{
  return dot(g,g)+l;
}

/** specialization of laplacian with complex g & l
 */
template<typename T, unsigned D>
inline T laplacian(const TinyVector<std::complex<T>,D>& g, const std::complex<T>& l)
{
  return  l.real()+OTCDot<T,T,D>::apply(g,g);
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
  ss<<i;
  return ss.str();
}


template<typename T>
struct BareKineticEnergy: public QMCHamiltonianBase
{

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
  Array<TraceReal,1>* T_sample;
  Array<TraceComp,1>* T_sample_comp;
  Array<TraceComp,2>* p_sample;
#endif
  ParticleSet& Ps;

  /** constructor
   *
   * Kinetic operators need to be re-evaluated during optimization.
   */
  BareKineticEnergy(RealType m=1.0): SameMass(true),M(m),OneOver2M(0.5/m)
  {
    set_energy_domain(kinetic);
    set_quantum_domain(quantum);
    streaming_kinetic      = false;
    streaming_kinetic_comp = false;
    streaming_momentum     = false;
    UpdateMode.set(OPTIMIZABLE,1);
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
    UpdateMode.set(OPTIMIZABLE,1);
    SpeciesSet& tspecies(p.getSpeciesSet());
    MinusOver2M.resize(tspecies.size());
    int massind=tspecies.addAttribute("mass");
    SameMass=true;
    M=tspecies(massind,0);
    OneOver2M=0.5/M;
    for(int i=0; i<tspecies.size(); ++i)
    {
      SameMass&=(std::abs(tspecies(massind,i)-M)<1e-6);
      MinusOver2M[i]=-1.0/(2.0*tspecies(massind,i));
    }
  }
  ///destructor
  ~BareKineticEnergy() { }

  void resetTargetParticleSet(ParticleSet& P) { }


#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities() 
  {
    request.contribute_array(myName);
    request.contribute_array(myName+"_complex");
    request.contribute_array("momentum");
  }

  virtual void checkout_particle_quantities(TraceManager& tm)
  {
    streaming_particles =  request.streaming_array(myName) 
                        || request.streaming_array(myName+"_complex")
                        || request.streaming_array("momentum");
    if( streaming_particles)
    {
      T_sample      = tm.checkout_real<1>(myName,Ps);
      T_sample_comp = tm.checkout_complex<1>(myName+"_complex",Ps);
      p_sample      = tm.checkout_complex<2>("momentum",Ps,DIM);
    }
  }

  virtual void delete_particle_quantities()
  {
    if( streaming_particles)
    {
      delete T_sample;
      delete T_sample_comp;
      delete p_sample;
    }
  }
#endif


  inline Return_t evaluate(ParticleSet& P)
  {
#if !defined(REMOVE_TRACEMANAGER)
    if( streaming_particles)
    {
      Value = evaluate_sp(P);
    }
    else
#endif
      if(SameMass)
      {
        Value = Dot(P.G,P.G) + Sum(P.L);
        Value*=-OneOver2M;
      }
      else
      {
        Value=0.0;
        for(int i=0; i<MinusOver2M.size(); ++i)
        {
          T x=0.0;
          for(int j=P.first(i); j<P.last(i); ++j)
            x += laplacian(P.G[j],P.L[j]);
          Value += x*MinusOver2M[i];
        }
      }
    return Value;
  }


  inline Return_t
  evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }


#if !defined(REMOVE_TRACEMANAGER)
  inline Return_t evaluate_sp(ParticleSet& P)
  {
    Array<RealType,1>& T_samp = *T_sample;
    Array<std::complex<RealType>,1>& T_samp_comp = *T_sample_comp;
    Array<std::complex<RealType>,2>& p_samp = *p_sample;
    std::complex<RealType> t1=0.0;
    const RealType clambda(-OneOver2M);
    Value = 0.0;
    if(SameMass)
    {
      for(int i=0; i<P.getTotalNum(); i++)
      {
        t1 = P.L[i] + dot(P.G[i],P.G[i]);
        t1 = clambda*t1;
        T_samp(i) = real(t1);
        T_samp_comp(i) = t1;
        for(int d=0; d<DIM; ++d)
          p_samp(i,d) = P.G[i][d];
        Value += real(t1);
      }
    }
    else
    {
      for(int s=0; s<MinusOver2M.size(); ++s)
      {
        T mlambda = MinusOver2M[s];
        for(int i=P.first(s); i<P.last(s); ++i)
        {
          //t1 = mlambda*( P.L[i] + dot(P.G[i],P.G[i]) );
          t1 =  P.L[i] + dot(P.G[i],P.G[i]);
          t1 *= mlambda;
          T_samp(i) = real(t1);
          T_samp_comp(i) = t1;
          for(int d=0; d<DIM; ++d)
            p_samp(i,d) = P.G[i][d];
          Value += real(t1);
        }
      }
    }
#if defined(TRACE_CHECK)
    RealType Vnow = Value;
    RealType Vsum = T_samp.sum();
    RealType Vold = evaluate_orig(P);
    if(std::abs(Vsum-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"accumtest: BareKineticEnergy::evaluate()"<< std::endl;
      app_log()<<"accumtest:   tot:"<< Vnow << std::endl;
      app_log()<<"accumtest:   sum:"<< Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    if(std::abs(Vold-Vnow)>TraceManager::trace_tol)
    {
      app_log()<<"versiontest: BareKineticEnergy::evaluate()"<< std::endl;
      app_log()<<"versiontest:   orig:"<< Vold << std::endl;
      app_log()<<"versiontest:    mod:"<< Vnow << std::endl;
      APP_ABORT("Trace check failed");
    }
#endif
    return Value;
  }

#endif

  inline Return_t evaluate_orig(ParticleSet& P)
  {
    if(SameMass)
    {
      Value = Dot(P.G,P.G) + Sum(P.L);
      Value*=-OneOver2M;
    }
    else
    {
      Value=0.0;
      for(int i=0; i<MinusOver2M.size(); ++i)
      {
        T x=0.0;
        for(int j=P.first(i); j<P.last(i); ++j)
          x += laplacian(P.G[j],P.L[j]);
        Value += x*MinusOver2M[i];
      }
    }
    return Value;
  }

  /** implements the virtual function.
   *
   * Nothing is done but should check the mass
   */
  bool put(xmlNodePtr)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "Kinetic energy";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new BareKineticEnergy(*this);
  }

#ifdef QMC_CUDA
  ////////////////////////////////
  // Vectorized version for GPU //
  ////////////////////////////////
  // Nothing is done on GPU here, just copy into vector
  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy)
  {
    std::vector<Walker_t*> &walkers = W.WalkerList;
    for (int iw=0; iw<walkers.size(); iw++)
    {
      Walker_t &w = *(walkers[iw]);
      RealType KE = - OneOver2M * (Dot(w.G,w.G)+Sum(w.L));
      w.getPropertyBase()[NUMPROPERTIES+myIndex] = KE;
      LocalEnergy[iw] += KE;
    }
  }
#endif

  //Not used anymore
  //void evaluate(WalkerSetRef& W, ValueVectorType& LE) {
  //  for(int iw=0; iw< W.walkers(); iw++) {
  //    RealType ke = 0.0;
  //    for(int iat=0; iat< W.particles(); iat++) {
  //      ke += dot(W.G(iw,iat),W.G(iw,iat)) + W.L(iw,iat);
  //    }
  //    LE[iw] -= M*ke;
  //  }
  //}
};
}
#endif


