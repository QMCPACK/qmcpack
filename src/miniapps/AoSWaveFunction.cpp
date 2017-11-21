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
#include <iostream>
using namespace std;

namespace qmcplusplus
{
  AoSWaveFunction::AoSWaveFunction(ParticleSet& ions, ParticleSet& els)
  {
    FirstTime=true;

    d_ee=DistanceTable::add(els,DT_AOS);
    d_ie=DistanceTable::add(ions,els,DT_AOS);

    int ip=omp_get_thread_num();
    double r2_cut=std::min(6.4,double(els.Lattice.WignerSeitzRadius));
    J2=new J2OrbType(els,ip);
    buildJ2(*J2,r2_cut);
  }

  AoSWaveFunction::~AoSWaveFunction()
  {
    delete J2;
  }

  void AoSWaveFunction::evaluateLog(ParticleSet& P)
  {
    constexpr valT czero(0);
    P.G=czero;
    P.L=czero;
    LogValue=J2->evaluateLog(P,P.G,P.L);
    FirstTime=false;
  }

  FakeWaveFunctionBase::posT AoSWaveFunction::evalGrad(ParticleSet& P, int iat)
  {
    return J2->evalGrad(P,iat);
  }

  FakeWaveFunctionBase::valT AoSWaveFunction::ratioGrad(ParticleSet& P, int iat, posT& grad)
  {
    return J2->ratioGrad(P,iat,grad);
  }

  FakeWaveFunctionBase::valT AoSWaveFunction::ratio(ParticleSet& P, int iat)
  {
    return J2->ratio(P,iat);
  }
  void AoSWaveFunction::acceptMove(ParticleSet& P, int iat)
  {
    J2->acceptMove(P,iat);
  }

  void AoSWaveFunction::restore(int iat) 
  {
    J2->restore(iat);
  }

  void AoSWaveFunction::evaluateGL(ParticleSet& P)
  {
    constexpr valT czero(0);
    P.G=czero;
    P.L=czero;
    J2->evaluateGL(P);
  }
}
