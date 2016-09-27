//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: D. Das, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_QMC_WOSPARTICLES_H
#define OHMMS_QMC_WOSPARTICLES_H
#include <vector>
#include "Numerics/Spline3D/Config.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

struct WOSParticles
{

  int m_ions;
  int m_elcs;
  int m_qpts;
  double qwt;

  /// position of common centre
  posvec_t R0;

  /// position of all particles \vec{q}_{i}
  posarray_t R;

  /// the charge Qi of all particles
  std::vector<double> Q;

  std::vector<double> wt;

  ///Constructor
  WOSParticles(ParticleSet& ions,
               const ParticleSet& elcs)
  {
    int iz = ions.Species.addAttribute("charge");
    m_ions = ions.getTotalNum();
    m_elcs = elcs.getTotalNum();
    m_qpts = m_ions + m_elcs;
    R.resize(m_qpts);
    wt.resize(m_qpts);
    Q.resize(m_qpts,-1.0);
    R0 = 0.0;
    /// assign coordinates and charges for the ions
    for(int i = 0; i < m_ions; i++)
    {
      R[i] =  ions.R[i];
      Q[i] = ions.Species(iz,ions.GroupID[i]);
    }
  }

  inline void setP(const ParticleSet& P)
  {
    for(int i = m_ions; i < m_qpts; i++)
      R[i] = P.R[i-m_ions];
  }

  inline void calcwt()
  {
    qwt = 0.0;
    for(int i = 0; i < m_qpts; i++)
      qwt += Q[i] * wt[i];
    qwt *= 0.5;
    return;
  }


};
}
#endif
