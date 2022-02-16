//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "Configuration.h"
#include "CuspCorrection.h"
#include "SoaLocalizedBasisSet.h"
#include "SoaAtomicBasisSet.h"
#include "MultiQuinticSpline1D.h"
#include "Numerics/MinimizeOneDim.h"


namespace qmcplusplus
{
bool readCuspInfo(const std::string& cuspInfoFile,
                  const std::string& objectName,
                  int OrbitalSetSize,
                  Matrix<CuspCorrectionParameters>& info)
{
  bool success = true;
  int ncenter  = info.rows();
  app_log() << "Reading cusp info from : " << cuspInfoFile << std::endl;
  Libxml2Document adoc;
  if (!adoc.parse(cuspInfoFile))
  {
    app_log() << "Could not find precomputed cusp data for spo set: " << objectName << std::endl;
    app_log() << "Recalculating data.\n";
    return false;
  }
  xmlNodePtr head = adoc.getRoot();
  head            = head->children;
  xmlNodePtr cur  = NULL, ctr;
  while (head != NULL)
  {
    std::string cname(getNodeName(head));
    if (cname == "sposet")
    {
      std::string name;
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add(name, "name");
      spoAttrib.put(head);
      if (name == objectName)
      {
        cur = head;
        break;
      }
    }
    head = head->next;
  }
  if (cur == NULL)
  {
    app_log() << "Could not find precomputed cusp data for spo set: " << objectName << std::endl;
    app_log() << "Recalculating data.\n";
    return false;
  }
  else
  {
    app_log() << "Found precomputed cusp data for spo set: " << objectName << std::endl;
  }
  cur = cur->children;
  while (cur != NULL)
  {
    std::string cname(getNodeName(cur));
    if (cname == "center")
    {
      int num = -1;
      OhmmsAttributeSet Attrib;
      Attrib.add(num, "num");
      Attrib.put(cur);
      if (num < 0 || num >= ncenter)
      {
        APP_ABORT("Error with cusp info xml block. incorrect center number. \n");
      }
      ctr = cur->children;
      while (ctr != NULL)
      {
        std::string cname(getNodeName(ctr));
        if (cname == "orbital")
        {
          int orb = -1;
          OhmmsAttributeSet orbAttrib;
          QMCTraits::RealType a1(0.0), a2, a3, a4, a5, a6, a7, a8, a9;
          orbAttrib.add(orb, "num");
          orbAttrib.add(a1, "redo");
          orbAttrib.add(a2, "C");
          orbAttrib.add(a3, "sg");
          orbAttrib.add(a4, "rc");
          orbAttrib.add(a5, "a1");
          orbAttrib.add(a6, "a2");
          orbAttrib.add(a7, "a3");
          orbAttrib.add(a8, "a4");
          orbAttrib.add(a9, "a5");
          orbAttrib.put(ctr);
          if (orb < OrbitalSetSize)
          {
            info(num, orb).redo     = a1;
            info(num, orb).C        = a2;
            info(num, orb).sg       = a3;
            info(num, orb).Rc       = a4;
            info(num, orb).alpha[0] = a5;
            info(num, orb).alpha[1] = a6;
            info(num, orb).alpha[2] = a7;
            info(num, orb).alpha[3] = a8;
            info(num, orb).alpha[4] = a9;
          }
        }
        ctr = ctr->next;
      }
    }
    cur = cur->next;
  }
  return success;
}

void broadcastCuspInfo(CuspCorrectionParameters& param, Communicate& Comm, int root)
{
#ifdef HAVE_MPI
  std::vector<double> buffer(9);
  buffer[0] = param.Rc;
  buffer[1] = param.C;
  buffer[2] = param.sg;
  buffer[3] = param.alpha[0];
  buffer[4] = param.alpha[1];
  buffer[5] = param.alpha[2];
  buffer[6] = param.alpha[3];
  buffer[7] = param.alpha[4];
  buffer[8] = param.redo;

  Comm.comm.broadcast(buffer.begin(), buffer.end(), root);

  param.Rc       = buffer[0];
  param.C        = buffer[1];
  param.sg       = buffer[2];
  param.alpha[0] = buffer[3];
  param.alpha[1] = buffer[4];
  param.alpha[2] = buffer[5];
  param.alpha[3] = buffer[6];
  param.alpha[4] = buffer[7];
  param.redo     = buffer[8] == 0.0 ? 0 : 1;
#endif
}

void splitPhiEta(int center, const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi, LCAOrbitalSet& Eta)
{
  using RealType = QMCTraits::RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);
  std::vector<bool> correct_this_center(corrCenter.size(), false);
  correct_this_center[center] = corrCenter[center];

  Phi.myBasisSet->queryOrbitalsForSType(correct_this_center, is_s_orbital);

  int nOrbs = Phi.getOrbitalSetSize();
  int bss   = Phi.getBasisSetSize();

  for (int i = 0; i < bss; i++)
  {
    if (is_s_orbital[i])
    {
      auto& cref(*(Eta.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0; //Eta->C(k,i) = 0.0;
    }
    else
    {
      auto& cref(*(Phi.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0; //Phi->C(k,i) = 0.0;
    }
  }
}

void removeSTypeOrbitals(const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi)
{
  using RealType = QMCTraits::RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);

  Phi.myBasisSet->queryOrbitalsForSType(corrCenter, is_s_orbital);

  int nOrbs = Phi.getOrbitalSetSize();
  int bss   = Phi.getBasisSetSize();

  for (int i = 0; i < bss; i++)
  {
    if (is_s_orbital[i])
    {
      auto& cref(*(Phi.C));
      for (int k = 0; k < nOrbs; k++)
        cref(k, i) = 0.0;
    }
  }
}


// Will be the corrected value for r < rc and the original wavefunction for r > rc
void computeRadialPhiBar(ParticleSet* targetP,
                         ParticleSet* sourceP,
                         int curOrb_,
                         int curCenter_,
                         SPOSet* Phi,
                         Vector<QMCTraits::RealType>& xgrid,
                         Vector<QMCTraits::RealType>& rad_orb,
                         const CuspCorrectionParameters& data)
{
  OneMolecularOrbital phiMO(targetP, sourceP, Phi);
  phiMO.changeOrbital(curCenter_, curOrb_);
  CuspCorrection cusp(data);

  for (int i = 0; i < xgrid.size(); i++)
  {
    rad_orb[i] = cusp.phiBar(xgrid[i], phiMO);
  }
}

using RealType = QMCTraits::RealType;

// Get the ideal local energy at one point
// Eq. 17 in the paper.  Coefficients are taken from the paper.
RealType getOneIdealLocalEnergy(RealType r, RealType Z, RealType beta0)
{
  RealType beta[7] = {3.25819, -15.0126, 33.7308, -42.8705, 31.2276, -12.1316, 1.94692};
  RealType idealEL = beta0;
  RealType r1      = r * r;
  for (int i = 0; i < 7; i++)
  {
    idealEL += beta[i] * r1;
    r1 *= r;
  }
  return idealEL * Z * Z;
}

// Get the ideal local energy for a vector of positions
void getIdealLocalEnergy(const ValueVector& pos, RealType Z, RealType Rc, RealType ELorigAtRc, ValueVector& ELideal)
{
  // assert(pos.size() == ELideal.size()
  RealType beta0 = 0.0;
  RealType tmp   = getOneIdealLocalEnergy(Rc, Z, beta0);
  beta0          = (ELorigAtRc - tmp) / (Z * Z);
  for (int i = 0; i < pos.size(); i++)
  {
    ELideal[i] = getOneIdealLocalEnergy(pos[i], Z, beta0);
  }
}

// Evaluate constraints. Equations 9-13 in the paper.
void evalX(RealType valRc,
           GradType gradRc,
           ValueType lapRc,
           RealType Rc,
           RealType Z,
           RealType C,
           RealType valAtZero,
           RealType eta0,
           TinyVector<ValueType, 5>& X)
{
  X[0] = std::log(std::abs(valRc - C));
  X[1] = gradRc[0] / (valRc - C);
  X[2] = (lapRc - 2.0 * gradRc[0] / Rc) / (valRc - C);
  X[3] = -Z * (valAtZero + eta0) / (valAtZero - C);
  X[4] = std::log(std::abs(valAtZero - C));
}

// Compute polynomial coefficients from constraints.  Eq. 14 in the paper.
void X2alpha(const TinyVector<ValueType, 5>& X, RealType Rc, TinyVector<ValueType, 5>& alpha)
{
  RealType RcInv = 1.0 / Rc, RcInv2 = RcInv * RcInv;
  alpha[0] = X[4];
  alpha[1] = X[3];
  alpha[2] = 6.0 * X[0] * RcInv2 - 3.0 * X[1] * RcInv + X[2] * 0.5 - 3.0 * X[3] * RcInv - 6.0 * X[4] * RcInv2 -
      0.5 * X[1] * X[1];
  alpha[3] = -8.0 * X[0] * RcInv2 * RcInv + 5.0 * X[1] * RcInv2 - X[2] * RcInv + 3.0 * X[3] * RcInv2 +
      8.0 * X[4] * RcInv2 * RcInv + X[1] * X[1] * RcInv;
  alpha[4] = 3.0 * X[0] * RcInv2 * RcInv2 - 2.0 * X[1] * RcInv2 * RcInv + 0.5 * X[2] * RcInv2 - X[3] * RcInv2 * RcInv -
      3.0 * X[4] * RcInv2 * RcInv2 - 0.5 * X[1] * X[1] * RcInv2;
}

// Eq. 16 in the paper.
RealType getZeff(RealType Z, RealType etaAtZero, RealType phiBarAtZero) { return Z * (1.0 + etaAtZero / phiBarAtZero); }

// Compute the effective one-electron local energy at a vector of points.
// Eq. 15 in the paper for r < Rc.  Normal local energy for R > Rc.
void getCurrentLocalEnergy(const ValueVector& pos,
                           RealType Zeff,
                           RealType Rc,
                           RealType originalELatRc,
                           CuspCorrection& cusp,
                           OneMolecularOrbital& phiMO,
                           ValueVector& ELcurr)
{
  // assert(pos.size() == ELcurr.size());
  ValueType val;
  GradType grad;
  ValueType lap;
  phiMO.phi_vgl(Rc, val, grad, lap);
  RealType dE = originalELatRc - (-0.5 * lap / val - Zeff / Rc);
  for (int i = 0; i < pos.size(); i++)
  {
    RealType r = pos[i];
    // prevent NaN's if phiBar is zero
    RealType offset = 1e-12;
    if (r <= Rc)
    {
      RealType dp = cusp.dpr(r);
      ELcurr[i]   = -0.5 * cusp.Rr(r) * (2.0 * dp / r + cusp.d2pr(r) + dp * dp) / (offset + cusp.phiBar(r, phiMO)) -
          Zeff / r + dE;
    }
    else
    {
      phiMO.phi_vgl(pos[i], val, grad, lap);
      ELcurr[i] = -0.5 * lap / val - Zeff / r + dE;
    }
  }
}

// Return value is local energy at Rc
RealType getOriginalLocalEnergy(const ValueVector& pos,
                                RealType Zeff,
                                RealType Rc,
                                OneMolecularOrbital& phiMO,
                                ValueVector& ELorig)
{
  // assert(pos.size() == ELorig.size());

  ValueType val;
  GradType grad;
  ValueType lap;
  for (int i = 0; i < pos.size(); i++)
  {
    RealType r = pos[i];
    phiMO.phi_vgl(r, val, grad, lap);
    ELorig[i] = -0.5 * lap / val - Zeff / r;
  }

  phiMO.phi_vgl(Rc, val, grad, lap);
  return -0.5 * lap / val - Zeff / Rc;
}

// Sum of squares difference between the current local energy and the ideal local energy.
//  This is the objective function to minimize.
RealType getELchi2(const ValueVector& ELcurr, const ValueVector& ELideal)
{
  assert(ELcurr.size() == ELideal.size());

  RealType chi2 = 0.0;
  for (int i = 0; i < ELcurr.size(); i++)
  {
    RealType diff = ELcurr[i] - ELideal[i];
    chi2 += diff * diff;
  }
  return chi2;
}

struct ValGradLap
{
  ValueType val;
  GradType grad;
  ValueType lap;
};

//  Compute the chi squared distance given a value for phi at zero.
RealType evaluateForPhi0Body(RealType phi0,
                             ValueVector& pos,
                             ValueVector& ELcurr,
                             ValueVector& ELideal,
                             CuspCorrection& cusp,
                             OneMolecularOrbital& phiMO,
                             ValGradLap phiAtRc,
                             RealType etaAtZero,
                             RealType ELorigAtRc,
                             RealType Z)
{
  cusp.cparam.sg = phi0 > 0.0 ? 1.0 : -1.0;
  cusp.cparam.C  = (phiAtRc.val * phi0 < 0.0) ? 1.5 * phiAtRc.val : 0.0;
  TinyVector<ValueType, 5> X;
  evalX(phiAtRc.val, phiAtRc.grad, phiAtRc.lap, cusp.cparam.Rc, Z, cusp.cparam.C, phi0, etaAtZero, X);
  X2alpha(X, cusp.cparam.Rc, cusp.cparam.alpha);
  RealType Zeff = getZeff(Z, etaAtZero, cusp.phiBar(0.0, phiMO));
  getCurrentLocalEnergy(pos, Zeff, cusp.cparam.Rc, ELorigAtRc, cusp, phiMO, ELcurr);
  RealType chi2 = getELchi2(ELcurr, ELideal);
  return chi2;
}

// Optimize free parameter (value of phi at zero) to minimize distance to ideal local energy.
// Output is return value and parameter values are in cusp.cparam
RealType minimizeForPhiAtZero(CuspCorrection& cusp,
                              OneMolecularOrbital& phiMO,
                              RealType Z,
                              RealType eta0,
                              ValueVector& pos,
                              ValueVector& ELcurr,
                              ValueVector& ELideal,
                              RealType start_phi0)
{
  ValGradLap vglAtRc;
  ValueVector tmp_pos(0);
  ValueVector ELorig(0);
  RealType Zeff = getZeff(Z, eta0, cusp.phiBar(0.0, phiMO));

  RealType ELorigAtRc = getOriginalLocalEnergy(tmp_pos, Zeff, cusp.cparam.Rc, phiMO, ELorig);
  getIdealLocalEnergy(pos, Z, cusp.cparam.Rc, ELorigAtRc, ELideal);
  phiMO.phi_vgl(cusp.cparam.Rc, vglAtRc.val, vglAtRc.grad, vglAtRc.lap);

  Bracket_min_t<RealType> bracket(start_phi0, 0.0, 0.0, false);
  try
  {
    bracket = bracket_minimum(
        [&](RealType x) -> RealType {
          return evaluateForPhi0Body(x, pos, ELcurr, ELideal, cusp, phiMO, vglAtRc, eta0, ELorigAtRc, Z);
        },
        start_phi0);
  }
  catch (const std::runtime_error& e)
  {
    APP_ABORT("Bracketing minimum failed for finding phi0. \n");
  }

  auto min_res = find_minimum(
      [&](RealType x) -> RealType {
        return evaluateForPhi0Body(x, pos, ELcurr, ELideal, cusp, phiMO, vglAtRc, eta0, ELorigAtRc, Z);
      },
      bracket);

  start_phi0 = min_res.first;

  return min_res.second;
}


// Optimize the cutoff radius.  There is an inner loop optimizing for phi0 for each value of Rc.
// Elcurr and ELideal are expected to have the correct size on input (same size as pos)
// Output is parameter values in cusp.cparam
void minimizeForRc(CuspCorrection& cusp,
                   OneMolecularOrbital& phiMO,
                   RealType Z,
                   RealType Rc_init,
                   RealType Rc_max,
                   RealType eta0,
                   ValueVector& pos,
                   ValueVector& ELcurr,
                   ValueVector& ELideal)
{
  Bracket_min_t<RealType> bracket(Rc_init, 0.0, 0.0, false);
  RealType start_phi0 = phiMO.phi(0.0);
  try
  {
    bracket = bracket_minimum(
        [&](RealType x) -> RealType {
          cusp.cparam.Rc = x;
          return minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal, start_phi0);
        },
        Rc_init, Rc_max);
  }
  catch (const std::runtime_error& e)
  {
    APP_ABORT("Bracketing minimum failed for finding rc. \n");
  }


  if (bracket.success)
  {
    auto min_res = find_minimum(
        [&](RealType x) -> RealType {
          cusp.cparam.Rc = x;
          return minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal, start_phi0);
        },
        bracket);
  }
  else
  {
    cusp.cparam.Rc = bracket.a;
    minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal, start_phi0);
  }
}


}; // namespace qmcplusplus
