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
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "LCAOrbitalSet.h"

class Communicate;
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
  using ValueType = QMCTraits::ValueType;
  using RealType  = QMCTraits::RealType;

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

  CuspCorrectionParameters() : Rc(0.0), C(0.0), sg(1.0), alpha(0.0), redo(0) {}
};


/// Broadcast cusp correction parameters
void broadcastCuspInfo(CuspCorrectionParameters& param, Communicate& Comm, int root);

class OneMolecularOrbital
{
  using RealType    = QMCTraits::RealType;
  using ValueType   = QMCTraits::ValueType;
  using GradType    = QMCTraits::GradType;
  using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;
  using GradVector  = OrbitalSetTraits<ValueType>::GradVector;
  using SPOSetPtr   = SPOSet*;

public:
  RealType phi(RealType r)
  {
    TinyVector<RealType, 3> dr = 0;
    dr[0]                      = r;

    targetPtcl->R[0] = sourcePtcl->R[curCenter];
    targetPtcl->makeMove(0, dr);
    Psi1->evaluateValue(*targetPtcl, 0, val1);

    return val1[curOrb];
  }

  void phi_vgl(RealType r, RealType& val, GradType& grad, RealType& lap)
  {
    TinyVector<RealType, 3> dr = 0;
    dr[0]                      = r;

    targetPtcl->R[0] = sourcePtcl->R[curCenter];
    targetPtcl->makeMove(0, dr);
    Psi1->evaluateVGL(*targetPtcl, 0, val1, grad1, lap1);

    val  = val1[curOrb];
    grad = grad1[curOrb];
    lap  = lap1[curOrb];
  }

  OneMolecularOrbital(ParticleSet* targetP, ParticleSet* sourceP, SPOSetPtr Phi)
      : targetPtcl(targetP), sourcePtcl(sourceP), curOrb(0), curCenter(0)
  {
    Psi1     = Phi;
    int norb = Psi1->getOrbitalSetSize();
    val1.resize(norb);
    grad1.resize(norb);
    lap1.resize(norb);
  }

  void changeOrbital(int centerIdx, int orbIdx)
  {
    curCenter = centerIdx;
    curOrb    = orbIdx;
  }

private:
  /// Temporary storage for real wavefunction values
  ValueVector val1;
  GradVector grad1;
  ValueVector lap1;

  /// target ParticleSet
  ParticleSet* targetPtcl;
  /// source ParticleSet
  ParticleSet* sourcePtcl;

  /// Index of orbital
  int curOrb;

  /// Index of atomic center
  int curCenter;

  SPOSetPtr Psi1;
};

/// Formulas for applying the cusp correction

class CuspCorrection
{
  using RealType = QMCTraits::RealType;

public:
  inline RealType phiBar(RealType r, OneMolecularOrbital& phiMO)
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
    return alpha[1] + 2.0 * alpha[2] * r + 3.0 * alpha[3] * r * r + 4.0 * alpha[4] * r * r * r;
  }

  inline RealType d2pr(RealType r)
  {
    auto& alpha = cparam.alpha;
    return 2.0 * alpha[2] + 6.0 * alpha[3] * r + 12.0 * alpha[4] * r * r;
  }

  CuspCorrection(const CuspCorrectionParameters& param) : cparam(param) {}

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

using RealType    = QMCTraits::RealType;
using ValueType   = QMCTraits::ValueType;
using GradType    = QMCTraits::GradType;
using ValueVector = OrbitalSetTraits<ValueType>::ValueVector;

/** Ideal local energy at one point
 * @param r  input radial distance
 * @param Z  nuclear charge
 * @param beta0  adjustable parameter to make energy continuous at Rc
 */
RealType getOneIdealLocalEnergy(RealType r, RealType Z, RealType beta0);

/** Ideal local energy at a vector of points
 * @param pos input vector of radial distances
 * @param Z nuclear charge
 * @param Rc cutoff radius where the correction meets the actual orbital
 * @param ELorigAtRc local energy at Rc.  beta0 is adjusted to make energy continuous at Rc
 * @param ELideal - output the ideal local energy at pos values
 */
void getIdealLocalEnergy(const ValueVector& pos, RealType Z, RealType Rc, RealType ELorigAtRc, ValueVector& ELideal);

/** Evaluate various orbital quantities that enter as constraints on the correction
 * @param valRc  orbital value at Rc
 * @param gradRc  orbital gradient at Rc
 * @param lapRc  orbital laplacian at Rc
 * @param Rc cutoff radius
 * @param Z nuclear charge
 * @param C offset to keep correction to a single sign
 * @param valAtZero orbital value at zero
 * @param eta0 value of non-corrected pieces of the orbital at zero
 * @param X output
 */
void evalX(RealType valRc,
           GradType gradRc,
           ValueType lapRc,
           RealType Rc,
           RealType Z,
           RealType C,
           RealType valAtZero,
           RealType eta0,
           TinyVector<ValueType, 5>& X);

/** Convert constraints to polynomial parameters
 * @param X input from evalX
 * @param Rc cutoff radius
 * @param alpha output the polynomial parameters for the correction
 */
void X2alpha(const TinyVector<ValueType, 5>& X, RealType Rc, TinyVector<ValueType, 5>& alpha);

/** Effective nuclear charge to keep effective local energy finite at zero
 * @param Z nuclear charge
 * @param etaAtZero value of non-S orbitals at this center
 * @param phiBarAtZero value of corrected orbital at zero
 */
RealType getZeff(RealType Z, RealType etaAtZero, RealType phiBarAtZero);

/**  Compute effective local energy at vector of points
 * @param pos input vector of radial distances
 * @param Zeff effective charge from getZeff
 * @param Rc cutoff radius
 * @param originalELatRc  Local energy at the center from the uncorrected orbital
 * @param cusp cusp correction parameters
 * @param phiMO uncorrected orbital (S-orbitals on this center only)
 * @param ELcurr output local energy at each distance in pos
 */
void getCurrentLocalEnergy(const ValueVector& pos,
                           RealType Zeff,
                           RealType Rc,
                           RealType originalELatRc,
                           CuspCorrection& cusp,
                           OneMolecularOrbital& phiMO,
                           ValueVector& ELcurr);

/** Local energy from uncorrected orbital
 * @param pos input vector of radial distances
 * @param Zeff nuclear charge
 * @param Rc cutoff radius
 * @param phiMO uncorrected orbital (S-orbitals on this center only)
 * @param ELorig output local energy at each distance in pos
 *
 * Return is value of local energy at zero.  This is the value needed for subsequent computations.
 * The routine can be called with an empty vector of positions to get just this value.
 */
RealType getOriginalLocalEnergy(const ValueVector& pos,
                                RealType Zeff,
                                RealType Rc,
                                OneMolecularOrbital& phiMO,
                                ValueVector& Elorig);

/** Sum of squares difference between the current and ideal local energies
 * This is the objective function to be minimized.
 * @param Elcurr  current local energy
 * @param Elideal  ideal local energy
 */
RealType getELchi2(const ValueVector& ELcurr, const ValueVector& ELideal);


/** Minimize chi2 with respect to phi at zero for a fixed Rc
 * @param cusp correction parameters
 * @param phiMO uncorrected orbital (S-orbitals on this center only)
 * @param Z nuclear charge
 * @param eta0 value at zero for parts of the orbital that don't require correction - the non-S-orbitals on this center and all orbitals on other centers
 * @param pos vector of radial positions
 * @param Elcurr storage for current local energy
 * @param Elideal storage for ideal local energy
 */
RealType minimizeForPhiAtZero(CuspCorrection& cusp,
                              OneMolecularOrbital& phiMO,
                              RealType Z,
                              RealType eta0,
                              ValueVector& pos,
                              ValueVector& ELcurr,
                              ValueVector& ELideal,
                              RealType start_phi0);


/** Minimize chi2 with respect to Rc and phi at zero.
 * @param cusp correction parameters
 * @param phiMO uncorrected orbital (S-orbitals on this center only)
 * @param Z nuclear charge
 * @param Rc_init initial value for Rc
 * @param Rc_max maximum value for Rc
 * @param eta0 value at zero for parts of the orbital that don't require correction - the non-S-orbitals on this center and all orbitals on other centers
 * @param pos vector of radial positions
 * @param Elcurr storage for current local energy
 * @param Elideal storage for ideal local energy
 *
 * Output is parameter values in cusp.cparam
 */
void minimizeForRc(CuspCorrection& cusp,
                   OneMolecularOrbital& phiMO,
                   RealType Z,
                   RealType Rc_init,
                   RealType Rc_max,
                   RealType eta0,
                   ValueVector& pos,
                   ValueVector& ELcurr,
                   ValueVector& ELideal);


} // namespace qmcplusplus

#endif
