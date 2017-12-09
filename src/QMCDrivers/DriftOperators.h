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

template<class T, class TG, unsigned D>
inline T getDriftScale(T tau, const ParticleAttrib<TinyVector<std::complex<TG>,D> >& ga)
{
  APP_ABORT("getDriftScale: Scaled drift computable only with real wavefunction forces.\n")
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

/** evaluate \f$\gamma\f$ for \f$ \bar V= \gamma V\f$
 *
 * Using eq. 34 of JCP 99, 2865 (1993)
 * \f$ \bar v(i)=\gamma_i \frac{-1+\sqrt{1+2*a*v^2*\tau}}{av^2\tau} v(i)\f$
 */
template<class T, unsigned D>
inline T getNodeCorrectionP(T tau, const ParticleAttrib<TinyVector<T,D> >& ga, T a=1)
{
  APP_ABORT("getNodeCorrectionP: unexpected deprecated code path.  Contact a developer.\n");
  T norm=0.0, norm_scaled=0.0;
  for(int i=0; i<ga.size(); ++i)
  {
    T vsq=dot(ga[i],ga[i]);
    T x=a*vsq*tau;
    T scale= (vsq<std::numeric_limits<T>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
    norm_scaled+=vsq*scale*scale;
    norm+=vsq;
  }
  return std::sqrt(norm_scaled/norm);
}

/** evaluate \f$\gamma\f$ for \f$ \bar V= \gamma V\f$
 *
 * Similar to getNodeCorrectionP but scale all the gradients with the same factor
 */
template<class T, unsigned D>
inline T getNodeCorrectionW(T tau, const ParticleAttrib<TinyVector<T,D> >& ga)
{
  APP_ABORT("getNodeCorrectionW: unexpected deprecated code path.  Contact a developer.\n");
  T vsq=Dot(ga,ga);
  T x=tau*vsq;
  return (vsq<std::numeric_limits<T>::epsilon())? 1.0:((-1.0+std::sqrt(1.0+2.0*x))/x);
}


/** evaluate drift using the scaling function by JCP93
 * @param tau timestep
 * @param qf quantum force
 * @param drift drift
 */
template<class T, unsigned D>
inline void setScaledDriftPbyP(T tau,
                               const ParticleAttrib<TinyVector<T,D> >& qf,
                               ParticleAttrib<TinyVector<T,D> >& drift)
{
  APP_ABORT("setScaledDriftPbyP: unexpected deprecated code path.  Contact a developer.\n");
  for(int iat=0; iat<qf.size(); ++iat)
  {
    T vsq=dot(qf[iat],qf[iat]);
    T sc=(vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
    drift[iat]=sc*qf[iat];
  }
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
    //T vsq=dot(qf[iat],qf[iat]);
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
 * @param qf quantum forces
 * @param drift scaled quantum forces
 * @param return correction term
 *
 * Fill the drift vector one particle at a time (pbyp).
 *
 * The naive drift is tau/mass*grad_psi_over_psi,
 *  namely the log derivative of the guiding wavefunction; tau is timestep; mass is particle mass
 * The naive drift diverges at a node, causing persistent configurations.
 * The norm of the drift vector should be limited in two ways:
 *  1. Umrigar: suppress drift divergence to mimic wf divergence
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
  T norm=0.0, norm_scaled=0.0; // variables to be accumulated
  for(int iat=0; iat<massinv.size(); ++iat)
  {
    // effectively set the Umrigar "a" parameter to particle mass
    T tau_over_mass = tau_au*massinv[iat]; // !!!! assume timestep is scaled by mass
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

template<class T, unsigned D>
inline T setLargestScaledDriftPbyP(T tau,
                                   const ParticleAttrib<TinyVector<T,D> >& qf,
                                   ParticleAttrib<TinyVector<T,D> >& drift)
{
  APP_ABORT("setLargestScaledDriftPbyP: unexpected deprecated code path.  Contact a developer.\n");
  T maxSC=tau;
  for(int iat=0; iat<qf.size(); ++iat)
  {
    T vsq=dot(qf[iat],qf[iat]);
    T sc=(vsq<std::numeric_limits<T>::epsilon())? tau:((-1.0+std::sqrt(1.0+2.0*tau*vsq))/vsq);
    maxSC=(sc<maxSC? sc:maxSC);
  }
  for(int iat=0; iat<qf.size(); ++iat)
  {
    drift[iat]=maxSC*qf[iat];
  }
  return maxSC;
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

/** da = scaled(tau)*ga
 * @param tau time step
 * @param qf complex quantum forces
 * @param drift drift
 */
template<class T, unsigned D>
inline void setScaledDrift(T tau,
                           const ParticleAttrib<TinyVector<std::complex<T>,D> >& qf,
                           ParticleAttrib<TinyVector<std::complex<T>,D> >& drift)
{
  APP_ABORT("setScaledDrift: untested code pathway.  Please contact a developer");
  T s = getDriftScale(tau,qf);
  ///INCOMPLETE implementation
  //PAOps<T,D>::scale(s,qf,drift);
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
