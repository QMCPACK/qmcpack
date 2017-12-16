//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HEPRESETHF_H
#define QMCPLUSPLUS_HEPRESETHF_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/DummyBasisSet.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/**
 *@brief A specialized class that implments Helium Clementi-Roetti
 *type orbitals.
 *
 *The single-particle orbital \f$ \psi \f$ is represented by
 \f[
 \psi ({\bf r}_j) = \sum_i C_{i} \chi_{i}({\bf r}_j-{\bf R}),
 \f]
 *where \f$ \chi_{i} \f$ is the ith basis function and
 *\f$ {\bf R} \f$ is the position of the ion.
 *
 *The Clementi-Rotti basis has the form
 \f[ \chi_{i}(r) = e^{-Z_{i}r} \f]
 *
 *For reference see: E Clementi & C Roetti, Atomic Data &
 *Nuclear Data Tables, 14 177 (1974).
 *
 *@note Used to debug the package.
 */
struct HePresetHF: public QMCTraits
{

  typedef DummyBasisSet                     BasisSet_t;
  enum {N=5};
  ///array containing the coefficients
  TinyVector<RealType,N> C;
  ///array containing the coefficients
  TinyVector<RealType,N> Z, ZZ;
  ///constructor
  HePresetHF()
  {
    C[0]=0.76838;
    C[1]=0.22346;
    C[2]=0.04082;
    C[3]=-0.00994;
    C[4]=0.00230;
    Z[0]=1.41714;
    Z[1]=2.37682;
    Z[2]=4.39628;
    Z[3]=6.52699;
    Z[4]=7.94252;
    const RealType fourpi = 4.0*(4.0*std::atan(1.0));
    for(int i=0; i<N; i++)
    {
      C[i] *= std::sqrt(std::pow(2.0*Z[i],3.0)/(2.0*fourpi));
      //C[i] *= sqrt(pow(2.0*Z[i],3.0)/2)/fourpi;
    }
    for(int i=0; i<N; i++)
      ZZ[i] = Z[i]*Z[i];
  }

  inline void reset() { }

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(d_table->origin(),P);
  }

  template<class VV>
  inline
  void
  evaluate(const ParticleSet& P, int iat, VV& phi)
  {
    RealType r(d_table->Temp[0].r1);
    RealType chi(0.0);
    for(int i=0; i<C.size(); i++)
    {
      chi += C[i]*std::exp(-Z[i]*r);
    }
    phi[0]=chi;
  }

  /** @ingroup particlebyparticle
   *@brief evalaute the single-particle-orbital values
   *
   *HeSTOClementiRottie class should not be used with particle-by-particle update
   */
  template<class VV, class GV>
  inline
  void
  evaluate(const ParticleSet& P, int iat, VV& phi, GV& dphi, VV& d2phi )
  {
    RealType r = d_table->Temp[0].r1;
    RealType rinv = d_table->Temp[0].rinv1;
    PosType dr = d_table->Temp[0].dr1;
    RealType chi = 0.0, d2chi = 0.0;
    PosType dchi;
    for(int i=0; i<C.size(); i++)
    {
      RealType u = C[i]*std::exp(-Z[i]*r);
      RealType du = -u*Z[i]*rinv; // 1/r du/dr
      RealType d2u = u*ZZ[i]; // d2u/dr2
      chi += u;
      dchi += du*dr;
      d2chi += (d2u+2.0*du);
    }
    phi[0]=chi;
    dphi[0] = dchi;
    d2phi[0]=d2chi;
  }

  /**
   *@param P input configuration containing N particles
   *@param first index of the first particle
   *@param last index of the last particle
   *@param logdet matrix \f$ logdet[0,0] =
   \sum_i C_{i} \chi_{i}({\bf r}_0-{\bf R}) \f$
   *@param dlogdet vector matrix \f$ dlogdet[0,0] =
   \sum_i C_{i} \nabla_j \phi_{i}({\bf r}_0-{\bf R}) \f$
   *@param d2logdet matrix \f$ d2logdet[0,0] =
   \sum_i C_{i} \nabla^2_j \phi_{i}({\bf r}_0-{\bf R}) \f$
   *@brief Evaluate the single-particle orbital and the derivatives.
   *
   \f[ \chi_{i}(r) = e^{-Z_{i}r} \f]
   \f[
   \nabla \chi_{i}(r) = \frac{d\chi_{i}}{dr}\hat{r} =
   -Z_{i}e^{-Z_{i}r}\hat{r}
   \f]  \f[
   \nabla^2 \chi_{i}(r) = \frac{1}{r^2}\frac{d}{dr}\left( r^2
   \frac{d\chi_{i}}{dr} \right) = Z_{i}^2 e^{-Z_{i}r}
   - \frac{2Z_{i}}{r} e^{-Z_{i}r}
   \f]
  */

  template<class VM, class GM>
  inline void
  evaluate(const ParticleSet& P, int first, int last,
           VM& logdet, GM& dlogdet, VM& d2logdet)
  {
    RealType r = d_table->r(first);
    RealType rinv = d_table->rinv(first);
    PosType dr = d_table->dr(first);
    RealType rinv2 = rinv*rinv;
    RealType chi = 0.0, d2chi = 0.0;
    PosType dchi;
    for(int i=0; i<C.size(); i++)
    {
      RealType u = C[i]*std::exp(-Z[i]*r);
      RealType du = -u*Z[i]*rinv; // 1/r du/dr
      RealType d2u = u*ZZ[i]; // d2u/dr2
      chi += u;
      dchi += du*dr;
      d2chi += (d2u+2.0*du);
    }
    logdet(0,0) =chi;
    dlogdet(0,0) = dchi;
    d2logdet(0,0) = d2chi;
  }

  template<class VM, class GM>
  inline void
  evaluate(const WalkerSetRef& W, int first, int last,
           std::vector<VM>& logdet, std::vector<GM>& dlogdet, std::vector<VM>& d2logdet)
  {
    int nptcl = last-first;
    for(int iw=0; iw<W.walkers(); iw++)
    {
      int nn = first;///first pair of the particle subset
      RealType r = d_table->r(iw,first);
      RealType rinv = d_table->rinv(iw,first);
      PosType dr = d_table->dr(iw,first);
      RealType rinv2 = rinv*rinv;
      RealType chi = 0.0, d2chi = 0.0;
      PosType dchi;
      for(int i=0; i<C.size(); i++)
      {
        RealType u = C[i]*std::exp(-Z[i]*r);
        RealType du = -u*Z[i]*rinv; // 1/r du/dr
        RealType d2u = u*ZZ[i]; // d2u/dr2
        chi += u;
        dchi += du*dr;
        d2chi += (d2u+2.0*du);
      }
      logdet[iw](0,0) =chi;
      dlogdet[iw](0,0) = dchi;
      d2logdet[iw](0,0) = d2chi;
    }
  }
///set the distance table for the single particle orbital
  void setTable(DistanceTableData* dtable)
  {
    d_table = dtable;
  }
  DistanceTableData* d_table;
};


struct HePresetHFBuilder: public OrbitalBuilderBase
{

  HePresetHFBuilder(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions);
  bool put(xmlNodePtr cur);

};

}
#endif
