//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <omp.h>
#include <miniapps/FakeWaveFunction.h>
#include <miniapps/input.hpp>

namespace qmcplusplus
{
  SoAWaveFunction::SoAWaveFunction(ParticleSet& ions, ParticleSet& els)
  {
    FirstTime=true;

    ions.RSoA=ions.R;
    els.RSoA=els.R;

    d_ee=DistanceTable::add(els,DT_SOA);
    d_ie=DistanceTable::add(ions,els,DT_SOA);

    double r2_cut=std::min(6.4,double(els.Lattice.WignerSeitzRadius));

    int ip=omp_get_thread_num();
    J2=new J2OrbType(els,ip);
    buildJ2(*J2,r2_cut);
  }

  SoAWaveFunction::~SoAWaveFunction()
  {
    delete J2;
  }

  void SoAWaveFunction::evaluateLog(ParticleSet& P)
  {
    constexpr valT czero(0);
    P.G=czero;
    P.L=czero;
    LogValue=J2->evaluateLog(P,P.G,P.L);
    FirstTime=false;
  }

  FakeWaveFunctionBase::posT SoAWaveFunction::evalGrad(ParticleSet& P, int iat)
  {
    return J2->evalGrad(P,iat);
  }

  FakeWaveFunctionBase::valT SoAWaveFunction::ratioGrad(ParticleSet& P, int iat, posT& grad)
  {
    return J2->ratioGrad(P,iat,grad);
  }

  FakeWaveFunctionBase::valT SoAWaveFunction::ratio(ParticleSet& P, int iat)
  {
    return J2->ratio(P,iat);
  }
  void SoAWaveFunction::acceptMove(ParticleSet& P, int iat)
  {
    J2->acceptMove(P,iat);
  }
  void SoAWaveFunction::restore(int iat) {}

  void SoAWaveFunction::evaluateGL(ParticleSet& P)
  {
    constexpr valT czero(0);
    P.G=czero;
    P.L=czero;
    J2->evaluateGL(P,P.G,P.L);
  }
}
