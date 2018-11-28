//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file CuspCorrection.h
  * @brief Corrections to electron-nucleus cusp for all-electron molecular calculations.
  */

#ifndef QMCPLUSPLUS_CUSPCORRECTION_H
#define QMCPLUSPLUS_CUSPCORRECTION_H

#include <cmath>
#include <iostream>
#include "Configuration.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/ParticleSet.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaAtomicBasisSet.h"
#include "QMCWaveFunctions/lcao/MultiQuinticSpline1D.h"
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"
#include "QMCWaveFunctions/lcao/SoaSphericalTensor.h"


namespace qmcplusplus
{
/**
  * @brief Cusp correction parameters
  *
  * From "Scheme for adding electron-nuclear cusps to Gaussian orbitals"  Ma, Towler, Drummond, and Needs
  *  JCP 122, 224322 (2005)
  *
  * Equations 7 and 8 in the paper define the correction.  These are the parameters in those equations.
  */

struct CuspCorrectionParameters
{
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::RealType RealType;

  /// The cutoff radius
  RealType Rc;

  /// A shift to keep correction to a single sign
  RealType C;

  /// The sign of the wavefunction at the nucleus
  RealType sg;

  /// The coefficients of the polynomial \f$p(r)\f$ in Eq 8
  TinyVector<ValueType, 5> alpha;

  /// Flag to indicate the correction should be recalculated
  int redo;

  CuspCorrectionParameters() : Rc(0.0), C(0.0), sg(1.0), redo(0), alpha(0.0) {}
};

class OneMolecularOrbital
{
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::ValueType ValueType;
  typedef QMCTraits::GradType GradType;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t GradVector_t;
  typedef SPOSet* SPOSetPtr;
public:
  RealType phi(RealType r)
  {
    TinyVector<RealType, 3> dr = 0;
    dr[0]                      = r;

    targetPtcl->R[0]             = sourcePtcl->R[curCenter];
    TinyVector<RealType, 3> ddr2 = targetPtcl->makeMove(0, dr);
    Psi1->evaluate(*targetPtcl, 0, val1);

    return val1[curOrb];
  }

  void phi_vgl(RealType r, RealType &val, GradType &grad, RealType &lap)
  {
    TinyVector<RealType, 3> dr = 0;
    dr[0]                      = r;

    targetPtcl->R[0]             = sourcePtcl->R[curCenter];
    TinyVector<RealType, 3> ddr2 = targetPtcl->makeMove(0, dr);
    Psi1->evaluate(*targetPtcl, 0, val1, grad1, lap1);

    val = val1[curOrb];
    grad = grad1[curOrb];
    lap = lap1[curOrb];
  }

  OneMolecularOrbital(ParticleSet* targetP, ParticleSet* sourceP, SPOSetPtr Phi) : targetPtcl(targetP), sourcePtcl(sourceP), curOrb(0), curCenter(0) {
    Psi1     = Phi;
    int norb = Psi1->OrbitalSetSize;
    val1.resize(norb);
    grad1.resize(norb);
    lap1.resize(norb);
  }

  void changeOrbital(int centerIdx, int orbIdx) {
    curCenter = centerIdx;
    curOrb = orbIdx;
    
  }

private:
  /// Index of atomic center
  int curCenter;

  /// Index of orbital
  int curOrb;

  /// Temporary storage for real wavefunction values
  ValueVector_t val1;
  GradVector_t grad1;
  ValueVector_t lap1;

  /// target ParticleSet
  ParticleSet* targetPtcl;
  /// source ParticleSet
  ParticleSet* sourcePtcl;

  SPOSetPtr Psi1;
};

/// Formulas for applying the cusp correction

class CuspCorrection
{
  typedef QMCTraits::RealType RealType;

public:
  inline RealType phiBar(RealType r, OneMolecularOrbital &phiMO)
  {
    if (r <= cparam.Rc)
      return cparam.C + Rr(r);
    else
      return phiMO.phi(r);
  }

  inline RealType Rr(RealType r) { return cparam.sg * std::exp(pr(r)); }

  inline RealType pr(RealType r)
  {
    auto& alpha = cparam.alpha;
    return alpha[0] + alpha[1] * r + alpha[2] * r * r + alpha[3] * r * r * r + alpha[4] * r * r * r * r;
  }

  inline RealType dpr(RealType r)
  {
    auto& alpha = cparam.alpha;
    return alpha[1] + 2.0*alpha[2] * r + 3.0*alpha[3]* r * r + 4.0*alpha[4] * r * r * r;
  }

  inline RealType d2pr(RealType r)
  {
    auto& alpha = cparam.alpha;
    return 2.0*alpha[2] + 6.0*alpha[3] * r + 12.0*alpha[4] * r * r;
  }

  CuspCorrection(const CuspCorrectionParameters &param) : cparam(param) {}

  CuspCorrectionParameters cparam;
};

/// Read cusp correction parameters from XML file
bool readCuspInfo(const std::string& cuspInfoFile,
                  const std::string& objectName,
                  int OrbitalSetSize,
                  Matrix<CuspCorrectionParameters>& info);

/// Divide molecular orbital into atomic S-orbitals on this center (phi), and everything else (eta).
void splitPhiEta(int center, const std::vector<bool>& corrCenter, LCAOrbitalSet& phi, LCAOrbitalSet& eta);

/// Remove S atomic orbitals from all molecular orbitals on all centers.
void removeSTypeOrbitals(const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi);

/// Compute the radial part of the corrected wavefunction
void computeRadialPhiBar(ParticleSet* targetP,
                         ParticleSet* sourceP,
                         int curOrb_,
                         int curCenter_,
                         SPOSet* Phi,
                         Vector<QMCTraits::RealType>& xgrid,
                         Vector<QMCTraits::RealType>& rad_orb,
                         const CuspCorrectionParameters& data);

typedef QMCTraits::RealType RealType;
typedef QMCTraits::ValueType ValueType;
typedef QMCTraits::GradType GradType;
typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;

RealType getOneIdealLocalEnergy(RealType r, RealType Z, RealType beta0);

void getIdealLocalEnergy(const ValueVector_t& pos, RealType Z, RealType Rc, RealType phiAtRc, ValueVector_t& ELideal);

void evalX(RealType valRc, GradType gradRc, ValueType lapRc, RealType Rc, RealType Z, RealType C,
           RealType valAtZero, RealType eta0, TinyVector<ValueType, 5> &X);

void X2alpha(const TinyVector<ValueType, 5> &X, RealType Rc, TinyVector<ValueType, 5> &alpha);

RealType getZeff(RealType Z, RealType etaAtZero, RealType phiBarAtZero);

void getCurrentLocalEnergy(const ValueVector_t& pos, RealType Zeff, RealType Rc, RealType originalELatRc, CuspCorrection &cusp, OneMolecularOrbital& phiMO, ValueVector_t& ELcurr);

RealType getOriginalLocalEnergy(const ValueVector_t& pos, RealType Zeff, RealType Rc, OneMolecularOrbital &phiMO, ValueVector_t& Elorig);

RealType getELchi2(const ValueVector_t& ELcurr, const ValueVector_t& ELideal);


RealType minimizeForPhiAtZero(CuspCorrection &cusp, OneMolecularOrbital &phiMO, RealType Z, RealType eta0, ValueVector_t &pos, ValueVector_t &ELcurr, ValueVector_t& ELideal);


void minimizeForRc(CuspCorrection &cusp, OneMolecularOrbital &phiMO, RealType Z, RealType Rc_init, RealType Rc_max, RealType eta0, ValueVector_t &pos,
ValueVector_t &ELcurr, ValueVector_t& ELideal);


} // namespace qmcplusplus

#endif
