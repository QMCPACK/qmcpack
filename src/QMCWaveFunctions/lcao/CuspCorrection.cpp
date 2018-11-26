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
#include "QMCWaveFunctions/lcao/SoaLocalizedBasisSet.h"
#include "QMCWaveFunctions/lcao/SoaAtomicBasisSet.h"
#include "QMCWaveFunctions/lcao/MultiQuinticSpline1D.h"
#include "QMCWaveFunctions/lcao/SoaCartesianTensor.h"
#include "Numerics/MinimizeOneDim.h"


namespace qmcplusplus
{
bool readCuspInfo(const std::string& cuspInfoFile,
                  const std::string& objectName,
                  int OrbitalSetSize,
                  Matrix<CuspCorrectionParameters>& info)
{
  bool success = true;
  std::string cname;
  int ncenter = info.rows();
  int nOrbs   = info.cols();
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
    getNodeName(cname, head);
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
    getNodeName(cname, cur);
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
        getNodeName(cname, ctr);
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

void splitPhiEta(int center, const std::vector<bool>& corrCenter, LCAOrbitalSet& Phi, LCAOrbitalSet& Eta)
{
  typedef QMCTraits::RealType RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);
  std::vector<bool> correct_this_center(corrCenter.size(), false);
  correct_this_center[center] = corrCenter[center];

  Phi.myBasisSet->queryOrbitalsForSType(correct_this_center, is_s_orbital);

  int nOrbs = Phi.OrbitalSetSize;
  int bss   = Phi.BasisSetSize;

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
  typedef QMCTraits::RealType RealType;

  std::vector<bool> is_s_orbital(Phi.myBasisSet->BasisSetSize, false);

  Phi.myBasisSet->queryOrbitalsForSType(corrCenter, is_s_orbital);

  int nOrbs = Phi.OrbitalSetSize;
  int bss   = Phi.BasisSetSize;

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
  //CuspCorrection cusp(targetP, sourceP);
  //cusp.setPsi(Phi);
  //cusp.cparam    = data;
  //cusp.curOrb    = curOrb_;
  //cusp.curCenter = curCenter_;

  for (int i = 0; i < xgrid.size(); i++)
  {
    rad_orb[i] = cusp.phiBar(xgrid[i], phiMO);
  }
}

typedef QMCTraits::RealType RealType;
RealType getOneIdealLocalEnergy(RealType r, RealType Z, RealType beta0)
{
  RealType beta[7] = {3.25819, -15.0126, 33.7308, -42.8705, 31.2276, -12.1316, 1.94692};
  RealType idealEL = beta0;
  RealType r1 = r*r;
  for (int i = 0; i < 7; i++) {
    idealEL += beta[i]*r1;
    r1 *= r;
  }
  return idealEL*Z*Z;
}

void
getIdealLocalEnergy(const ValueVector_t& pos, RealType Z, RealType Rc, RealType ELorigAtRc, ValueVector_t& ELideal)
{
  // assert(pos.size() == ELideal.size()
  RealType beta0 = 0.0;
  RealType tmp = getOneIdealLocalEnergy(Rc, Z, beta0);
  beta0 = (ELorigAtRc - tmp)/(Z*Z);
  for (int i = 0; i < pos.size(); i++) {
      ELideal[i] = getOneIdealLocalEnergy(pos[i], Z, beta0);
  }
}

void evalX(RealType valRc, GradType gradRc, ValueType lapRc, RealType Rc, RealType Z, RealType C,
           RealType valAtZero, RealType eta0, TinyVector<ValueType, 5> &X)
{
  X[0] = std::log(std::abs(valRc - C));
  X[1] = gradRc[0]/(valRc - C);
  X[2] = (lapRc - 2.0*gradRc[0]/Rc)/(valRc - C);
  X[3] = -Z*(valAtZero + eta0)/(valAtZero - C);
  X[4] = std::log(std::abs(valAtZero-C));
}

void X2alpha(const TinyVector<ValueType, 5> &X, RealType Rc, TinyVector<ValueType, 5> &alpha)
{
   RealType RcInv=1.0/Rc, RcInv2=RcInv*RcInv;
    alpha[0] = X[4];
    alpha[1] = X[3];
    alpha[2] = 6.0*X[0]*RcInv2 - 3.0*X[1]*RcInv + X[2]*0.5
               - 3.0*X[3]*RcInv - 6.0*X[4]*RcInv2 - 0.5*X[1]*X[1];
    alpha[3] = -8.0*X[0]*RcInv2*RcInv + 5.0*X[1]*RcInv2 - X[2]*RcInv
               + 3.0*X[3]*RcInv2 + 8.0*X[4]*RcInv2*RcInv + X[1]*X[1]*RcInv;
    alpha[4] = 3.0*X[0]*RcInv2*RcInv2 - 2.0*X[1]*RcInv2*RcInv
               + 0.5*X[2]*RcInv2 - X[3]*RcInv2*RcInv - 3.0*X[4]*RcInv2*RcInv2
               - 0.5*X[1]*X[1]*RcInv2;

}

RealType getZeff(RealType Z, RealType etaAtZero, RealType phiBarAtZero)
{
  return Z*(1.0 + etaAtZero/phiBarAtZero);
}

void getCurrentLocalEnergy(const ValueVector_t& pos, RealType Zeff, RealType Rc, RealType originalELatRc, CuspCorrection &cusp, OneMolecularOrbital &phiMO, ValueVector_t& ELcurr)
{
  // assert(pos.size() == ELcurr.size());
  ValueType val;
  GradType grad;
  ValueType lap;
  phiMO.phi_vgl(Rc, val, grad, lap);
  RealType dE = originalELatRc - (-0.5*lap/val - Zeff/Rc);
  //std::cout << "dE = " << dE << std::endl;
  for (int i = 0; i  < pos.size(); i++) {
    RealType r = pos[i];
    if (r <= Rc) {
      RealType dp = cusp.dpr(r);
      ELcurr[i] = -0.5*cusp.Rr(r)*(2.0*dp/r + cusp.d2pr(r) + dp*dp)/cusp.phiBar(r, phiMO) - Zeff/r + dE;
    } else {
      phiMO.phi_vgl(pos[i], val, grad, lap);
      ELcurr[i] = -0.5*lap/val - Zeff/r + dE;
    }
  }
}

// Returns value is local energy at Rc
RealType getOriginalLocalEnergy(const ValueVector_t& pos, RealType Zeff, RealType Rc, OneMolecularOrbital &phiMO, ValueVector_t& ELorig)
{
  // assert(pos.size() == ELorig.size());

  ValueType val;
  GradType grad;
  ValueType lap;
  for (int i = 0; i  < pos.size(); i++) {
    RealType r = pos[i];
    phiMO.phi_vgl(r, val, grad, lap);
    ELorig[i] = -0.5*lap/val - Zeff/r;
  }

  phiMO.phi_vgl(Rc, val, grad, lap);
  return -0.5*lap/val - Zeff/Rc;
  
}

RealType getELchi2(const ValueVector_t& ELcurr, const ValueVector_t& ELideal)
{
   assert(ELcurr.size() == ELideal.size());

  RealType chi2 = 0.0;
  for (int i = 0; i  < ELcurr.size(); i++) {
    RealType diff = ELcurr[i] - ELideal[i];
    chi2 += diff*diff;
  }
  return chi2;
}

class MinimizePhiAtZero
{
public:
  
  CuspCorrection& cusp;

  OneMolecularOrbital& phiMO;

  ValueType etaAtZero; // phiEta.phi(Rc);

  ValueType valAtRc; // phiMO.phi(Rc);
  GradType gradAtRc;
  ValueType lapAtRc;

  RealType Rc;
  RealType Z;

  RealType ELorigAtRc;

  ValueVector_t& pos;
  ValueVector_t& ELcurr;
  ValueVector_t& ELideal;

  MinimizePhiAtZero(ValueVector_t &pos_, ValueVector_t& ELcurr_, ValueVector_t &ELideal_, CuspCorrection &cusp_, OneMolecularOrbital &phiMO_) : pos(pos_), ELcurr(ELcurr_), ELideal(ELideal_), cusp(cusp_), phiMO(phiMO_) {}


  RealType operator()(RealType phi0) const {
    return evaluateForPhi0(phi0);
  }

  RealType evaluateForPhi0(RealType phi0) const
  {
      //std::cout << "Start cycle, phi0 =  " << phi0 << std::endl;
      cusp.cparam.sg = phi0 > 0.0 ? 1.0:-1.0;
      cusp.cparam.C = (valAtRc*phi0 < 0.0) ? 1.5*valAtRc:0.0;
      TinyVector<ValueType, 5> X;
      evalX(valAtRc, gradAtRc, lapAtRc, Rc, Z, cusp.cparam.C, phi0, etaAtZero, X);
      X2alpha(X, Rc, cusp.cparam.alpha);
      RealType Zeff = getZeff(Z, etaAtZero, cusp.phiBar(0.0, phiMO));
      getCurrentLocalEnergy(pos, Zeff, Rc, ELorigAtRc, cusp, phiMO, ELcurr);
      //getCurrentLocalEnergy(pos, Zeff, Rc, 0.0, cusp, ELcurr);
      //std::cout << "  ELideal = " << ELideal[0] << std::endl;
      //std::cout << "  ELcurr = " << ELcurr[0] << std::endl;
      RealType chi2 = getELchi2(ELcurr, ELideal);
      std::cout << "  phi0 = " << phi0 << " chi2 =  " << chi2 << std::endl;
      return chi2;
  }
};
struct ValGradLap
{
  ValueType val; // phiMO.phi(Rc);
  GradType grad;
  ValueType lap;
};

  RealType evaluateForPhi0Body(RealType phi0, ValueVector_t &pos, ValueVector_t &ELcurr, ValueVector_t &ELideal, CuspCorrection &cusp, OneMolecularOrbital &phiMO, ValGradLap phiAtRc, RealType etaAtZero, RealType ELorigAtRc, RealType Z)
  {
      //std::cout << "Start cycle, phi0 =  " << phi0 << std::endl;
      cusp.cparam.sg = phi0 > 0.0 ? 1.0:-1.0;
      cusp.cparam.C = (phiAtRc.val*phi0 < 0.0) ? 1.5*phiAtRc.val:0.0;
      TinyVector<ValueType, 5> X;
      ValueType  phiBarAtRc; // phiMO.phi(Rc);
      evalX(phiAtRc.val, phiAtRc.grad, phiAtRc.lap, cusp.cparam.Rc, Z, cusp.cparam.C, phi0, etaAtZero, X);
      X2alpha(X, cusp.cparam.Rc, cusp.cparam.alpha);
      RealType Zeff = getZeff(Z, etaAtZero, cusp.phiBar(0.0, phiMO));
      getCurrentLocalEnergy(pos, Zeff, cusp.cparam.Rc, ELorigAtRc, cusp, phiMO, ELcurr);
      //getCurrentLocalEnergy(pos, Zeff, Rc, 0.0, cusp, ELcurr);
      //std::cout << "  ELideal = " << ELideal[0] << std::endl;
      //std::cout << "  ELcurr = " << ELcurr[0] << std::endl;
      RealType chi2 = getELchi2(ELcurr, ELideal);
//      std::cout << "  phi0 = " << phi0 << " chi2 =  " << chi2 << std::endl;
      return chi2;
  }

// output is return value and parameter values in cusp.cparam
RealType minimizeForPhiAtZero(CuspCorrection &cusp, OneMolecularOrbital &phiMO, RealType Z, RealType eta0, ValueVector_t &pos, ValueVector_t &ELcurr, ValueVector_t& ELideal)
{
  //MinimizePhiAtZero minPhi0(pos, ELcurr, ELideal, cusp);

  ValGradLap vglAtRc;
  ValueVector_t tmp_pos(0);
  ValueVector_t ELorig(0);
  RealType Zeff = getZeff(Z, eta0, cusp.phiBar(0.0, phiMO) );

  RealType ELorigAtRc = getOriginalLocalEnergy(tmp_pos, Zeff, cusp.cparam.Rc, phiMO, ELorig);
  getIdealLocalEnergy(pos, Z, cusp.cparam.Rc, ELorigAtRc, ELideal);
#if 0
  minPhi0.ELorigAtRc =  getOriginalLocalEnergy(tmp_pos, Zeff, Rc, cusp.phiMO, ELorig);
  minPhi0.Rc = Rc;
  minPhi0.Z = Z;
  minPhi0.etaAtZero = eta0;
  cusp.phiMO.phi_vgl(Rc, minPhi0.valAtRc, minPhi0.gradAtRc, minPhi0.lapAtRc);
#endif
  phiMO.phi_vgl(cusp.cparam.Rc, vglAtRc.val, vglAtRc.grad, vglAtRc.lap);

  RealType start_phi0 = phiMO.phi(0.0);
  //std::cout << "start phi0 = " << start_phi0 << std::endl;
  Bracket_min_t<RealType> bracket = bracket_minimum([&](RealType x)->RealType{return evaluateForPhi0Body(x, pos, ELcurr, ELideal, cusp, phiMO, vglAtRc, eta0, ELorigAtRc, Z);}, start_phi0);
  //Bracket_min_t<RealType> bracket = bracket_minimum([&minPhi0](RealType x)->RealType{return minPhi0.oneCycle(x);}, start_phi0);
  //Bracket_min_t<RealType> bracket = bracket_minimum(minPhi0, start_phi0);
  //std::cout << "bracket okay = " << bracket.success << std::endl;
  //std::cout << "bracket = " << bracket.a << " " << bracket.b << " " << bracket.c << std::endl;

  //auto min_res = find_mininum([&minPhi0](RealType x)->RealType{return minPhi0.oneCycle(x);}, bracket);
  auto min_res = find_minimum([&](RealType x)->RealType{return evaluateForPhi0Body(x, pos, ELcurr, ELideal, cusp, phiMO, vglAtRc, eta0, ELorigAtRc, Z) ;}, bracket);
  //auto min_res = find_minimum(minPhi0, bracket);
  //std::cout << "phi0 min = " << min_res.first << " " << min_res.second  << std::endl;
  return min_res.second;
}


void minimizeForRc(CuspCorrection &cusp, OneMolecularOrbital &phiMO, RealType Z, RealType Rc_init, RealType Rc_max, RealType eta0, ValueVector_t &pos,
ValueVector_t &ELcurr, ValueVector_t& ELideal)
{
  RealType Rc = Rc_init;
  Bracket_min_t<RealType> bracket = bracket_minimum([&](RealType x)->RealType{cusp.cparam.Rc = x; return minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal);}, Rc_init, Rc_max);

//  std::cout << "rc bracket okay = " << bracket.success << std::endl;
//  std::cout << "rc bracket = " << bracket.a << " " << bracket.b << " " << bracket.c << std::endl;
  if (bracket.success) {
    auto min_res  = find_minimum([&](RealType x)->RealType{cusp.cparam.Rc = x; return minimizeForPhiAtZero(cusp, phiMO, Z , eta0, pos, ELcurr, ELideal);}, bracket);
//    std::cout << "rc min = " << min_res.first << " " << min_res.second  << std::endl;
  } else {
//    std::cout << "best rc = " << bracket.a << std::endl;
    cusp.cparam.Rc = bracket.a;
    minimizeForPhiAtZero(cusp, phiMO, Z, eta0, pos, ELcurr, ELideal);
  }
//  std::cout << "  alpha = " << cusp.cparam.alpha[0] << " " << cusp.cparam.alpha[1] << " " << cusp.cparam.alpha[2];
//  std::cout << " " << cusp.cparam.alpha[3] << " " << cusp.cparam.alpha[4] << std::endl;

}



}; // namespace qmcplusplus
