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
    
    
#ifndef OHMMS_QMC_WOSPOTENTIAL_H
#define OHMMS_QMC_WOSPOTENTIAL_H
#include <algorithm>
#include <vector>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Numerics/Spline3D/Config.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/WOS/Domain.h"
#include "QMCHamiltonians/WOS/Device.h"
#include "QMCHamiltonians/WOS/WOSParticles.h"

namespace qmcplusplus
{

struct WOSPotential: public QMCHamiltonianBase
{

  //    typedef double RealType;
  int Mode;
  int m_runs;
  int dmc_runs;
  double m_norm;
  Device* device;
  WOSParticles* WP;

  /// constructor
  WOSPotential(int mode,
               int mruns,
               int Druns,
               Device* adevice,
               ParticleSet& ions,
               ParticleSet& elcs)
  {
    WP = new WOSParticles(ions,elcs);
    device = adevice;
    set_mrun(mruns);
    Mode = mode;
    dmc_runs = Druns;
  }

  ~WOSPotential() {}

  inline void set_mrun(int mruns)
  {
    m_runs = mruns;
    m_norm = 1.0/double(mruns);
  }

  /// evaluate potential by random walk
  inline ValueType evaluate(ParticleSet& P)
  {
    switch(Mode)
    {
    case 0:
      return method0(P);
      break;
    case 1:
      return method1(P);
      break;
    case 2:
      return method2(P);
      break;
    case 3:
      return method3(P);
      break;
    case 4:
      return method4(P);
      break;
    case 5:
      return method5(P);
      break;
    case 6:
      return method6(P);
      break;
    default:
      return method0(P);
    }
  }

  inline ValueType evaluate(ParticleSet& P, RealType& x)
  {
    return x = evaluate(P);
  }

  ValueType method0(ParticleSet&);
  ValueType method1(ParticleSet&);
  ValueType method2(ParticleSet&);
  ValueType method3(ParticleSet&);
  ValueType method4(ParticleSet&);
  ValueType method5(ParticleSet&);
  ValueType method6(ParticleSet&);


#ifdef USE_FASTWALKER
  inline void evaluate(WalkerSetRef& W, ValueVectorType& LE) {}
#else
  inline void evaluate(WalkerSetRef& W, ValueVectorType& LE) {}
#endif


};
}
#endif
