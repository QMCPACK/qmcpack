//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#define QMCPLUSPLUS_QMCDRIFTOPERATORS_H
#include "type_traits/scalar_traits.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "ParticleBase/RandomSeqGenerator.h"
namespace qmcplusplus
{

template<class T, class TG, unsigned D>
inline T getDriftScale(T tau, const ParticleAttrib<TinyVector<TG,D> >& ga)
{
  T vsq=Dot(ga,ga);
  return (vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
}

/** evaluate a drift with a real force
 * @param tau timestep
 * @param qf quantum force
 * @param drift
 */
template<class Tt, class TG, class T, unsigned D>
inline void getScaledDrift(Tt tau, const TinyVector<TG,D>& qf, TinyVector<T,D>& drift)
{
  //We convert the complex gradient to real and temporarily store in drift.
  convert(qf,drift);
  T vsq=dot(drift,drift);
  T sc = (vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
  //Apply the umrigar scaled drift.
  drift*=sc;
}

/** scale drift
 * @param tau_au timestep au
 * @param qf quantum forces
 * @param drift scaled quantum forces
 * @param return correction term
 *
 * Assume, mass=1
 */
template<class T, class T1, unsigned D>
inline T setScaledDriftPbyPandNodeCorr(T tau,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  T norm=0.0, norm_scaled=0.0, tau2=tau*tau;
  for(int iat=0; iat<qf.size(); ++iat)
  {
    convert(qf[iat],drift[iat]);
    T vsq=dot(drift[iat],drift[iat]);
    T sc=(vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
    norm_scaled+=vsq*sc*sc;
    norm+=vsq*tau2;
    drift[iat]*=sc;
  }
  return std::sqrt(norm_scaled/norm);
}

/** scale drift
 * @param tau_au timestep au
 * @param massinv 1/m per particle
 * @param qf quantum forces, i.e. grad_psi_over_psi
 * @param drift scaled importance-sampling drift
 * @param return correction term
 *
 * Fill the drift vector one particle at a time (pbyp).
 *
 * The naive drift is tau/mass*|grad_psi_over_psi|,
 *  grad_psi_over_psi the log derivative of the guiding wavefunction;
 *  tau is timestep; mass is particle mass
 * The naive drift diverges at a node and is large near a nucleus.
 * Both cases may cause persistent configurations, because the
 *  reverse proposal probability is virtually zero.
 *
 * The norm of the drift vector should be limited in two ways:
 *  1. Umrigar: suppress drift divergence with a sclaing factor
 *  2. Ceperley: limit max drift rate to diffusion rate -> set Umrigar "a" parameter
 * The choice of drift vector does not affect VMC correctness
 *  so long as the proposal probabilities are correctly calculated.
 * The choice of drift vector changes the DMC Green's function. BE CAREFUL!
 *
 * T should be either float or double
 * T1 may be real (float or double) or complex<float or double>
 * D should be the number of spatial dimensions (int)
 */
template<class T, class T1, unsigned D>
inline T setScaledDriftPbyPandNodeCorr(T tau_au, const std::vector<T>& massinv,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  // Umrigar eq. (34) "a" parameter is set to 1.0
  T norm=0.0, norm_scaled=0.0; // variables to be accumulated
  for(int iat=0; iat<massinv.size(); ++iat)
  {
    // !!!! assume timestep is scaled by mass
    T tau_over_mass = tau_au*massinv[iat]; 
    // save real part of wf log derivative in drift
    convert(qf[iat],drift[iat]);
    T vsq = dot(drift[iat],drift[iat]);
    // calculate drift scalar "sc" of Umrigar, JCP 99, 2865 (1993); eq. (34) * tau
    // use naive drift if vsq may cause numerical instability in the denominator
    T sc  = (vsq < std::numeric_limits<T>::epsilon()) ? tau_over_mass : (-1.0+std::sqrt(1.0+2.0*tau_over_mass*vsq))/vsq;
    drift[iat] *= sc;

    norm_scaled+=vsq*sc*sc;
    norm+=vsq*tau_over_mass*tau_over_mass;
  }
  return std::sqrt(norm_scaled/norm);
}

//NOTE: While poorly named, setScaledDrift is the all-electron analogue of
//      getScaledDrift.

/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf real quantum forces
 * @param drift drift
 */
template<class T, class TG, unsigned D>
inline void setScaledDrift(T tau,
                           const ParticleAttrib<TinyVector<TG,D> >& qf,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  T s = getDriftScale(tau,qf);
  PAOps<T,D,TG>::scale(s,qf,drift);
}

/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf real quantum forces
 * @param drift drift
 */
template<class T, unsigned D>
inline void setScaledDrift(T tau,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  T s = getDriftScale(tau,drift);
  drift *=s;
}


/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf complex quantum forces
 * @param drift drift
 */
template<class T, class TG, unsigned D>
inline void setScaledDrift(T tau,
                           const ParticleAttrib<TinyVector<std::complex<TG>,D> >& qf,
                           ParticleAttrib<TinyVector<T,D> >& drift)
{
  for(int iat=0; iat<qf.size(); ++iat)
    convert(qf[iat],drift[iat]);

  T s = getDriftScale(tau,drift);
  drift*=s;
}

template<class T, class TG, unsigned D>
inline void assignDrift(T s,
                        const ParticleAttrib<TinyVector<TG,D> >& ga,
                        ParticleAttrib<TinyVector<T,D> >& da)
{
  PAOps<T,D,TG>::scale(s,ga,da);
}

template<class T, class TG, unsigned D>
inline void assignDrift(T s,
                        const ParticleAttrib<TinyVector<std::complex<TG>,D> >& ga,
                        ParticleAttrib<TinyVector<T,D> >& da)
{
  //This operation does s*ga, and takes the real part.
  PAOps<T,D,TG>::scale(s,ga,da);
}

//Assign drift does pbyp calculation of the scaled drift on all drift components.
template<class T, class T1, unsigned D>
inline void assignDrift(T tau_au, const std::vector<T>& massinv,
                                       const ParticleAttrib<TinyVector<T1,D> >& qf,
                                       ParticleAttrib<TinyVector<T,D> >& drift)
{
  for(int iat=0; iat<massinv.size(); ++iat)
  {
    T tau_over_mass=tau_au*massinv[iat];
    // naive drift "tau/mass*qf" can diverge
    getScaledDrift(tau_over_mass,qf[iat],drift[iat]);
  }
}

}
#endif
