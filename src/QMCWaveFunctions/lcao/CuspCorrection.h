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

  CuspCorrectionParameters() : Rc(0.0), C(0.0), sg(0.0), redo(0), alpha(0.0) {}
};


/// Formulas for applying the cusp correction

class CuspCorrection
{
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::ValueType ValueType;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef SPOSet* SPOSetPtr;

public:
  inline RealType phiBar(RealType r)
  {
    if (r <= cparam.Rc)
      return cparam.C + Rr(r);
    else
      return phi(r);
  }

  RealType phi(RealType r)
  {
    TinyVector<RealType, 3> dr = 0;
    dr[0]                      = r;

    targetPtcl->R[0]             = sourcePtcl->R[curCenter];
    TinyVector<RealType, 3> ddr2 = targetPtcl->makeMove(0, dr);
    Psi1->evaluate(*targetPtcl, 0, val1);

    return val1[curOrb];
  }

  inline RealType Rr(RealType r) { return cparam.sg * std::exp(pr(r)); }

  inline RealType pr(RealType r)
  {
    auto& alpha = cparam.alpha;
    return alpha[0] + alpha[1] * r + alpha[2] * r * r + alpha[3] * r * r * r + alpha[4] * r * r * r * r;
  }

  CuspCorrection(ParticleSet* targetP, ParticleSet* sourceP) : targetPtcl(targetP), sourcePtcl(sourceP) {}

  void setPsi(SPOSetPtr Phi)
  {
    Psi1     = Phi;
    int norb = Psi1->OrbitalSetSize;
    val1.resize(norb);
  }

  CuspCorrectionParameters cparam;

  /// Index of orbital
  int curOrb;

  /// Index of atomic center
  int curCenter;

  /// Temporary storage for real wavefunction values
  ValueVector_t val1;

  /// target ParticleSet
  ParticleSet* targetPtcl;
  /// source ParticleSet
  ParticleSet* sourcePtcl;

  SPOSetPtr Psi1;
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


} // namespace qmcplusplus

#endif
