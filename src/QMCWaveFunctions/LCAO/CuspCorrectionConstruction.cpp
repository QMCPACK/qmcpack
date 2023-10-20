//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "CuspCorrectionConstruction.h"
#include "Message/Communicate.h"
#include "SoaCuspCorrectionBasisSet.h"
#include "Utilities/FairDivide.h"
#include "SoaLocalizedBasisSet.h"
#include "SoaAtomicBasisSet.h"
#include "MultiQuinticSpline1D.h"
#include "Numerics/MinimizeOneDim.h"
#include "OhmmsData/AttributeSet.h"


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

void saveCusp(const std::string& filename, const Matrix<CuspCorrectionParameters>& info, const std::string& id)
{
  const int num_centers      = info.rows();
  const int orbital_set_size = info.cols();
  xmlDocPtr doc              = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr cuspRoot        = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  xmlNodePtr spo             = xmlNewNode(NULL, (const xmlChar*)"sposet");
  xmlNewProp(spo, (const xmlChar*)"name", (const xmlChar*)id.c_str());
  xmlAddChild(cuspRoot, spo);
  xmlDocSetRootElement(doc, cuspRoot);

  for (int center_idx = 0; center_idx < num_centers; center_idx++)
  {
    xmlNodePtr ctr = xmlNewNode(NULL, (const xmlChar*)"center");
    std::ostringstream num;
    num << center_idx;
    xmlNewProp(ctr, (const xmlChar*)"num", (const xmlChar*)num.str().c_str());

    for (int mo_idx = 0; mo_idx < orbital_set_size; mo_idx++)
    {
      std::ostringstream num0, C, sg, rc, a1, a2, a3, a4, a5;
      xmlNodePtr orb = xmlNewNode(NULL, (const xmlChar*)"orbital");
      num0 << mo_idx;
      xmlNewProp(orb, (const xmlChar*)"num", (const xmlChar*)num0.str().c_str());


      C.setf(std::ios::scientific, std::ios::floatfield);
      C.precision(14);
      C << info(center_idx, mo_idx).C;
      sg.setf(std::ios::scientific, std::ios::floatfield);
      sg.precision(14);
      sg << info(center_idx, mo_idx).sg;
      rc.setf(std::ios::scientific, std::ios::floatfield);
      rc.precision(14);
      rc << info(center_idx, mo_idx).Rc;
      a1.setf(std::ios::scientific, std::ios::floatfield);
      a1.precision(14);
      a1 << info(center_idx, mo_idx).alpha[0];
      a2.setf(std::ios::scientific, std::ios::floatfield);
      a2.precision(14);
      a2 << info(center_idx, mo_idx).alpha[1];
      a3.setf(std::ios::scientific, std::ios::floatfield);
      a3.precision(14);
      a3 << info(center_idx, mo_idx).alpha[2];
      a4.setf(std::ios::scientific, std::ios::floatfield);
      a4.precision(14);
      a4 << info(center_idx, mo_idx).alpha[3];
      a5.setf(std::ios::scientific, std::ios::floatfield);
      a5.precision(14);
      a5 << info(center_idx, mo_idx).alpha[4];
      xmlNewProp(orb, (const xmlChar*)"C", (const xmlChar*)C.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"sg", (const xmlChar*)sg.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"rc", (const xmlChar*)rc.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"a1", (const xmlChar*)a1.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"a2", (const xmlChar*)a2.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"a3", (const xmlChar*)a3.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"a4", (const xmlChar*)a4.str().c_str());
      xmlNewProp(orb, (const xmlChar*)"a5", (const xmlChar*)a5.str().c_str());
      xmlAddChild(ctr, orb);
    }
    xmlAddChild(spo, ctr);
  }

  app_log() << "Saving resulting cusp Info xml block to: " << filename << std::endl;
  xmlSaveFormatFile(filename.c_str(), doc, 1);
  xmlFreeDoc(doc);
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
    rad_orb[i] = phiBar(cusp, xgrid[i], phiMO);
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

RealType phiBar(const CuspCorrection& cusp, RealType r, OneMolecularOrbital& phiMO)
{
  if (r <= cusp.cparam.Rc)
    return cusp.cparam.C + cusp.Rr(r);
  else
    return phiMO.phi(r);
}

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
      ELcurr[i]   = -0.5 * cusp.Rr(r) * (2.0 * dp / r + cusp.d2pr(r) + dp * dp) / (offset + phiBar(cusp, r, phiMO)) -
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
  RealType Zeff = getZeff(Z, etaAtZero, phiBar(cusp, 0.0, phiMO));
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
  RealType Zeff = getZeff(Z, eta0, phiBar(cusp, 0.0, phiMO));

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

// Modifies orbital set lcwc
void applyCuspCorrection(const Matrix<CuspCorrectionParameters>& info,
                         ParticleSet& targetPtcl,
                         ParticleSet& sourcePtcl,
                         LCAOrbitalSet& lcao,
                         SoaCuspCorrection& cusp,
                         const std::string& id)
{
  const int num_centers      = info.rows();
  const int orbital_set_size = info.cols();
  using RealType             = QMCTraits::RealType;

  NewTimer& cuspApplyTimer = createGlobalTimer("CuspCorrectionConstruction::applyCuspCorrection", timer_level_medium);

  ScopedTimer cuspApplyTimerWrapper(cuspApplyTimer);

  LCAOrbitalSet phi("phi", std::unique_ptr<LCAOrbitalSet::basis_type>(lcao.myBasisSet->makeClone()));
  phi.setOrbitalSetSize(lcao.getOrbitalSetSize());

  LCAOrbitalSet eta("eta", std::unique_ptr<LCAOrbitalSet::basis_type>(lcao.myBasisSet->makeClone()));
  eta.setOrbitalSetSize(lcao.getOrbitalSetSize());

  std::vector<bool> corrCenter(num_centers, "true");

  //What's this grid's lifespan?  Why on the heap?
  auto radial_grid = std::make_unique<LogGrid<RealType>>();
  radial_grid->set(0.000001, 100.0, 1001);


  Vector<RealType> xgrid;
  Vector<RealType> rad_orb;
  xgrid.resize(radial_grid->size());
  rad_orb.resize(radial_grid->size());
  for (int ig = 0; ig < radial_grid->size(); ig++)
  {
    xgrid[ig] = radial_grid->r(ig);
  }

  for (int ic = 0; ic < num_centers; ic++)
  {
    *eta.C = *lcao.C;
    *phi.C = *lcao.C;

    splitPhiEta(ic, corrCenter, phi, eta);

    // loop over MO index - cot must be an array (of len MO size)
    //   the loop is inside cot - in the multiqunitic
    auto cot = std::make_unique<CuspCorrectionAtomicBasis<RealType>>();
    cot->initializeRadialSet(*radial_grid, orbital_set_size);
    //How is this useful?
    // cot->ID.resize(orbital_set_size);
    // for (int mo_idx = 0; mo_idx < orbital_set_size; mo_idx++) {
    //   cot->ID[mo_idx] = mo_idx;
    // }

    for (int mo_idx = 0; mo_idx < orbital_set_size; mo_idx++)
    {
      computeRadialPhiBar(&targetPtcl, &sourcePtcl, mo_idx, ic, &phi, xgrid, rad_orb, info(ic, mo_idx));
      RealType yprime_i = (rad_orb[1] - rad_orb[0]) / (radial_grid->r(1) - radial_grid->r(0));
      OneDimQuinticSpline<RealType> radial_spline(radial_grid->makeClone(), rad_orb);
      radial_spline.spline(0, yprime_i, rad_orb.size() - 1, 0.0);
      cot->addSpline(mo_idx, radial_spline);

      if (outputManager.isDebugActive())
      {
        // For testing against AoS output
        // Output phiBar to soaOrbs.downdet.C0.MO0
        int nElms   = 500;
        RealType dx = info(ic, mo_idx).Rc * 1.2 / nElms;
        Vector<RealType> pos;
        Vector<RealType> output_orb;
        pos.resize(nElms);
        output_orb.resize(nElms);
        for (int i = 0; i < nElms; i++)
        {
          pos[i] = (i + 1.0) * dx;
        }
        computeRadialPhiBar(&targetPtcl, &sourcePtcl, mo_idx, ic, &phi, pos, output_orb, info(ic, mo_idx));
        std::string filename = "soaOrbs." + id + ".C" + std::to_string(ic) + ".MO" + std::to_string(mo_idx);
        std::cout << "Writing to " << filename << std::endl;
        std::ofstream out(filename.c_str());
        out << "# r phiBar(r)" << std::endl;
        for (int i = 0; i < nElms; i++)
        {
          out << pos[i] << "  " << output_orb[i] << std::endl;
        }
        out.close();
      }
    }
    cusp.add(ic, std::move(cot));
  }
  removeSTypeOrbitals(corrCenter, lcao);
}

void generateCuspInfo(Matrix<CuspCorrectionParameters>& info,
                      const ParticleSet& targetPtcl,
                      const ParticleSet& sourcePtcl,
                      const LCAOrbitalSet& lcao,
                      const std::string& id,
                      Communicate& Comm)
{
  const int num_centers      = info.rows();
  const int orbital_set_size = info.cols();
  using RealType             = QMCTraits::RealType;

  NewTimer& cuspCreateTimer = createGlobalTimer("CuspCorrectionConstruction::createCuspParameters", timer_level_medium);
  NewTimer& splitPhiEtaTimer = createGlobalTimer("CuspCorrectionConstruction::splitPhiEta", timer_level_fine);
  NewTimer& computeTimer     = createGlobalTimer("CuspCorrectionConstruction::computeCorrection", timer_level_fine);

  ScopedTimer createCuspTimerWrapper(cuspCreateTimer);

  LCAOrbitalSet phi("phi", std::unique_ptr<LCAOrbitalSet::basis_type>(lcao.myBasisSet->makeClone()));
  phi.setOrbitalSetSize(lcao.getOrbitalSetSize());

  LCAOrbitalSet eta("eta", std::unique_ptr<LCAOrbitalSet::basis_type>(lcao.myBasisSet->makeClone()));
  eta.setOrbitalSetSize(lcao.getOrbitalSetSize());


  std::vector<bool> corrCenter(num_centers, "true");

  using GridType = OneDimGridBase<RealType>;
  int npts       = 500;

  // Parallelize correction of MO's across MPI ranks
  std::vector<int> offset;
  FairDivideLow(orbital_set_size, Comm.size(), offset);

  int start_mo = offset[Comm.rank()];
  int end_mo   = offset[Comm.rank() + 1];
  app_log() << "  Number of molecular orbitals to compute correction on this rank: " << end_mo - start_mo << std::endl;

// Specify dynamic scheduling explicitly for load balancing.   Each iteration should take enough
// time that scheduling overhead is not an issue.
#pragma omp parallel for schedule(dynamic) collapse(2)
  for (int center_idx = 0; center_idx < num_centers; center_idx++)
  {
    for (int mo_idx = start_mo; mo_idx < end_mo; mo_idx++)
    {
      ParticleSet localTargetPtcl(targetPtcl);
      ParticleSet localSourcePtcl(sourcePtcl);

      LCAOrbitalSet local_phi("local_phi", std::unique_ptr<LCAOrbitalSet::basis_type>(phi.myBasisSet->makeClone()));
      local_phi.setOrbitalSetSize(phi.getOrbitalSetSize());

      LCAOrbitalSet local_eta("local_eta", std::unique_ptr<LCAOrbitalSet::basis_type>(eta.myBasisSet->makeClone()));
      local_eta.setOrbitalSetSize(eta.getOrbitalSetSize());

#pragma omp critical
      app_log() << "   Working on MO: " << mo_idx << " Center: " << center_idx << std::endl;

      {
        ScopedTimer local_timer(splitPhiEtaTimer);

        *local_eta.C = *lcao.C;
        *local_phi.C = *lcao.C;
        splitPhiEta(center_idx, corrCenter, local_phi, local_eta);
      }

      bool corrO = false;
      auto& cref(*(local_phi.C));
      for (int ip = 0; ip < cref.cols(); ip++)
      {
        if (std::abs(cref(mo_idx, ip)) > 0)
        {
          corrO = true;
          break;
        }
      }

      if (corrO)
      {
        OneMolecularOrbital etaMO(&localTargetPtcl, &localSourcePtcl, &local_eta);
        etaMO.changeOrbital(center_idx, mo_idx);

        OneMolecularOrbital phiMO(&localTargetPtcl, &localSourcePtcl, &local_phi);
        phiMO.changeOrbital(center_idx, mo_idx);

        SpeciesSet& tspecies(localSourcePtcl.getSpeciesSet());
        int iz     = tspecies.addAttribute("charge");
        RealType Z = tspecies(iz, localSourcePtcl.GroupID[center_idx]);

        RealType Rc_max = 0.2;
        RealType rc     = 0.1;

        RealType dx = rc * 1.2 / npts;
        ValueVector pos(npts);
        ValueVector ELideal(npts);
        ValueVector ELcurr(npts);
        for (int i = 0; i < npts; i++)
        {
          pos[i] = (i + 1.0) * dx;
        }

        RealType eta0 = etaMO.phi(0.0);
        ValueVector ELorig(npts);
        CuspCorrection cusp(info(center_idx, mo_idx));
        {
          ScopedTimer local_timer(computeTimer);
          minimizeForRc(cusp, phiMO, Z, rc, Rc_max, eta0, pos, ELcurr, ELideal);
        }
        // Update shared object.  Each iteration accesses a different element and
        // this is an array (no bookkeeping data to update), so no synchronization
        // is necessary.
        info(center_idx, mo_idx) = cusp.cparam;
      }
    }
  }

  for (int root = 0; root < Comm.size(); root++)
  {
    int start_mo = offset[root];
    int end_mo   = offset[root + 1];
    for (int mo_idx = start_mo; mo_idx < end_mo; mo_idx++)
    {
      for (int center_idx = 0; center_idx < num_centers; center_idx++)
      {
        broadcastCuspInfo(info(center_idx, mo_idx), Comm, root);
      }
    }
  }
}

} // namespace qmcplusplus
