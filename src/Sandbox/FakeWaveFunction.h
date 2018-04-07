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
#ifndef QMCPLUSPLUS_FAKEWAVEFUNCTIONS_H
#define QMCPLUSPLUS_FAKEWAVEFUNCTIONS_H
#include <Configuration.h>
#include <QMCWaveFunctions/Jastrow/BsplineFunctor.h>
#include <QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h>
#include <QMCWaveFunctions/Jastrow/J2OrbitalSoA.h>

namespace qmcplusplus
{
  /** A minimal TrialWavefunction
   */
  struct FakeWaveFunctionBase
  {
    using valT=OHMMS_PRECISION;
    using posT=TinyVector<OHMMS_PRECISION,OHMMS_DIM>;

    valT LogValue;
    DistanceTableData* d_ee;
    DistanceTableData* d_ie;

    inline void setRmax(valT x)
    {
      d_ie->setRmax(x);
    }

    virtual ~FakeWaveFunctionBase(){}
    virtual void evaluateLog(ParticleSet& P)=0;
    virtual posT evalGrad(ParticleSet& P, int iat)=0;
    virtual valT ratioGrad(ParticleSet& P, int iat, posT& grad)=0;
    virtual valT ratio(ParticleSet& P, int iat)=0;
    virtual void acceptMove(ParticleSet& P, int iat)=0;
    virtual void restore(int iat)=0;
    virtual void evaluateGL(ParticleSet& P)=0;
  };

  struct AoSWaveFunction: public FakeWaveFunctionBase
  {

    using J2OrbType=TwoBodyJastrowOrbital<BsplineFunctor<valT> >;
    bool FirstTime;
    J2OrbType* J2;

    AoSWaveFunction(ParticleSet& ions, ParticleSet& els);
    ~AoSWaveFunction();
    void evaluateLog(ParticleSet& P);
    posT evalGrad(ParticleSet& P, int iat);
    valT ratioGrad(ParticleSet& P, int iat, posT& grad);
    valT ratio(ParticleSet& P, int iat);
    void acceptMove(ParticleSet& P, int iat);
    void restore(int iat);
    void evaluateGL(ParticleSet& P);
  };

  struct SoAWaveFunction: public FakeWaveFunctionBase
  {
    using J2OrbType=J2OrbitalSoA<BsplineFunctor<valT> >;

    bool FirstTime;
    J2OrbType* J2;

    SoAWaveFunction(ParticleSet& ions, ParticleSet& els);
    ~SoAWaveFunction();
    void evaluateLog(ParticleSet& P);
    posT evalGrad(ParticleSet& P, int iat);
    valT ratioGrad(ParticleSet& P, int iat, posT& grad);
    valT ratio(ParticleSet& P, int iat);
    void acceptMove(ParticleSet& P, int iat);
    void restore(int iat);
    void evaluateGL(ParticleSet& P);
  };
}
#endif

