//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file RPAJastrow.cpp
 * @brief Definitions of RPAJastrow
 */

#include "RPAJastrow.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
RPAJastrow::RPAJastrow(ParticleSet& target) : WaveFunctionComponent("RPAJastrow"), targetPtcl(target)
{
  Optimizable = true;
}

RPAJastrow::~RPAJastrow() = default;

bool RPAJastrow::put(xmlNodePtr cur)
{
  ReportEngine PRE("RPAJastrow", "put");
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  app_log() << "!!!  WARNING:  RPAJastrow is not fully tested for production !!!\n";
  app_log() << "!!!      level calculations.  Use at your own risk!          !!!\n";
  app_log() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
  //capture attribute jastrow/@name
  MyName           = "RPA_Jee";
  std::string useL = "yes";
  std::string useS = "yes";
  rpafunc          = "rpa";
  OhmmsAttributeSet a;
  a.add(MyName, "name");
  a.add(useL, "longrange");
  a.add(useS, "shortrange");
  a.add(rpafunc, "function");
  a.put(cur);
  Rs                = -1.0;
  Kc                = -1.0;
  std::string ID_Rs = "RPA_rs";
  ParameterSet params;
  params.add(Rs, "rs");
  params.add(Kc, "kc");
  params.put(cur);
  buildOrbital(MyName, useL, useS, rpafunc, Rs, Kc);
  return true;
}

void RPAJastrow::buildOrbital(const std::string& name,
                              const std::string& UL,
                              const std::string& US,
                              const std::string& RF,
                              RealType R,
                              RealType K)
{
  std::string ID_Rs = "RPA_rs";
  MyName            = name;
  std::string useL  = UL;
  std::string useS  = US;
  rpafunc           = RF;
  Rs                = R;
  Kc                = K;
  app_log() << std::endl << "   LongRangeForm is " << rpafunc << std::endl;
  DropLongRange  = (useL == "no");
  DropShortRange = (useS == "no");
  RealType tlen =
      std::pow(3.0 / 4.0 / M_PI * targetPtcl.Lattice.Volume / static_cast<RealType>(targetPtcl.getTotalNum()),
               1.0 / 3.0);
  if (Rs < 0)
  {
    if (targetPtcl.Lattice.SuperCellEnum)
    {
      Rs = tlen;
    }
    else
    {
      std::cout << "  Error finding rs. Is this an open system?!" << std::endl;
      Rs = 100.0;
    }
  }
  int indx      = targetPtcl.SK->getKLists().ksq.size() - 1;
  double Kc_max = std::pow(targetPtcl.SK->getKLists().ksq[indx], 0.5);
  if (Kc < 0)
  {
    Kc = 2.0 * std::pow(2.25 * M_PI, 1.0 / 3.0) / tlen;
  }
  if (Kc > Kc_max)
  {
    Kc = Kc_max;
    app_log() << "    Kc set too high. Resetting to the maximum value" << std::endl;
  }
  app_log() << "    RPAJastrowBuilder::addTwoBodyPart Rs = " << Rs << "  Kc= " << Kc << std::endl;
  if (rpafunc == "yukawa" || rpafunc == "breakup")
    myHandler = std::make_unique<LRHandlerTemp<YukawaBreakup<RealType>, LPQHIBasis>>(targetPtcl, Kc);
  else if (rpafunc == "rpa")
    myHandler = std::make_unique<LRRPAHandlerTemp<RPABreakup<RealType>, LPQHIBasis>>(targetPtcl, Kc);
  else if (rpafunc == "dyukawa")
    myHandler = std::make_unique<LRHandlerTemp<DerivYukawaBreakup<RealType>, LPQHIBasis>>(targetPtcl, Kc);
  else if (rpafunc == "drpa")
    myHandler = std::make_unique<LRRPAHandlerTemp<DerivRPABreakup<RealType>, LPQHIBasis>>(targetPtcl, Kc);
  else
    throw std::invalid_argument("RPAJastrowBuilder::buildOrbital:  Unrecognized rpa function type.\n");
  myHandler->Breakup(targetPtcl, Rs);
  app_log() << "  Maximum K shell " << myHandler->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << myHandler->Fk.size() << std::endl;
  if (!DropLongRange)
  {
    makeLongRange();
    app_log() << "  Using LongRange part" << std::endl;
  }
  if (!DropShortRange)
  {
    makeShortRange();
    app_log() << "  Using ShortRange part" << std::endl;
  }
}

void RPAJastrow::makeLongRange()
{
  // create two-body kSpaceJastrow
  kSpaceJastrow::SymmetryType oneBodySymm, twoBodySymm;
  bool oneBodySpin, twoBodySpin;
  oneBodySymm  = kSpaceJastrow::ISOTROPIC;
  twoBodySymm  = kSpaceJastrow::ISOTROPIC;
  oneBodySpin  = false;
  twoBodySpin  = false;
  auto LongRangeRPA_uptr =
      std::make_unique<kSpaceJastrow>(targetPtcl, targetPtcl, oneBodySymm, -1, "cG1", oneBodySpin, // no one-body part
                                      twoBodySymm, Kc, "cG2", twoBodySpin);
  LongRangeRPA = LongRangeRPA_uptr.get();
  // fill in CG2 coefficients
  std::vector<RealType> oneBodyCoefs, twoBodyCoefs;
  twoBodyCoefs.resize(myHandler->MaxKshell);
  //  need to cancel prefactor in kSpaceJastrow
  RealType prefactorInv = -targetPtcl.Lattice.Volume;
  for (size_t is = 0; is < myHandler->MaxKshell; is++)
  {
    twoBodyCoefs[is] = prefactorInv * myHandler->Fk_symm[is];
  }
  LongRangeRPA->setCoefficients(oneBodyCoefs, twoBodyCoefs);
  Psi.push_back(std::move(LongRangeRPA_uptr));
}

void RPAJastrow::makeShortRange()
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!  WARNING: tiny, nparam, npts should be input parameters  !!!
  // !!!                    fix in a future PR!                   !!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  app_log() << "  Adding Short Range part of RPA function" << std::endl;
  //short-range uses realHandler
  RealType tiny = 1e-6;
  Rcut          = myHandler->get_rc() - tiny;
  //create numerical functor of type BsplineFunctor<RealType>.
  auto nfunc_uptr = std::make_unique<FuncType>();
  nfunc           = nfunc_uptr.get();
  ShortRangePartAdapter<RealType> SRA(myHandler.get());
  SRA.setRmax(Rcut);
  auto j2 = std::make_unique<J2OrbitalSoA<BsplineFunctor<RealType>>>("RPA", targetPtcl);
  size_t nparam                              = 12;  // number of Bspline parameters
  size_t npts                                = 100; // number of 1D grid points for basis functions
  RealType cusp                              = SRA.df(0);
  RealType delta                             = Rcut / static_cast<double>(npts);
  std::vector<RealType> X(npts + 1), Y(npts + 1);
  for (size_t i = 0; i < npts; ++i)
  {
    X[i] = i * delta;
    Y[i] = SRA.evaluate(X[i]);
  }
  X[npts]              = npts * delta;
  Y[npts]              = 0.0;
  std::string functype = "rpa";
  std::string useit    = "no";
  nfunc->initialize(nparam, X, Y, cusp, Rcut + tiny, functype, useit);
  for (size_t i = 0; i < npts; ++i)
  {
    X[i] = i * delta;
    Y[i] = SRA.evaluate(X[i]);
  }
  j2->addFunc(0, 0, std::move(nfunc_uptr));
  ShortRangeRPA = j2.get();
  Psi.push_back(std::move(j2));
}

void RPAJastrow::resetParameters(const opt_variables_type& active)
{
  //This code was removed in 6 April 2017.  To reimplement, please consult a revision
  //earlier than this.
}

void RPAJastrow::checkOutVariables(const opt_variables_type& active) {}

void RPAJastrow::checkInVariables(opt_variables_type& active) {}

void RPAJastrow::reportStatus(std::ostream& os)
{
  for (int i = 0; i < Psi.size(); i++)
    Psi[i]->reportStatus(os);
}

RPAJastrow::LogValueType RPAJastrow::evaluateLog(const ParticleSet& P,
                                                 ParticleSet::ParticleGradient_t& G,
                                                 ParticleSet::ParticleLaplacian_t& L)
{
  log_value_ = 0.0;
  for (int i = 0; i < Psi.size(); i++)
    log_value_ += Psi[i]->evaluateLog(P, G, L);
  return log_value_;
}

RPAJastrow::PsiValueType RPAJastrow::ratio(ParticleSet& P, int iat)
{
  ValueType r(1.0);
  for (int i = 0; i < Psi.size(); i++)
    r *= Psi[i]->ratio(P, iat);
  return static_cast<PsiValueType>(r);
}

RPAJastrow::GradType RPAJastrow::evalGrad(ParticleSet& P, int iat)
{
  GradType grad(0);
  for (int i = 0; i < Psi.size(); i++)
    grad += Psi[i]->evalGrad(P, iat);
  return grad;
}

RPAJastrow::PsiValueType RPAJastrow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  ValueType r(1);
  for (int i = 0; i < Psi.size(); i++)
  {
    r *= Psi[i]->ratioGrad(P, iat, grad_iat);
  }
  return static_cast<PsiValueType>(r);
}


void RPAJastrow::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  for (int i = 0; i < Psi.size(); i++)
    Psi[i]->acceptMove(P, iat, safe_to_delay);
}

void RPAJastrow::restore(int iat)
{
  for (int i = 0; i < Psi.size(); i++)
    Psi[i]->restore(iat);
}

void RPAJastrow::registerData(ParticleSet& P, WFBufferType& buf)
{
  for (int i = 0; i < Psi.size(); i++)
    Psi[i]->registerData(P, buf);
}

RPAJastrow::LogValueType RPAJastrow::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{
  log_value_ = 0.0;
  for (int i = 0; i < Psi.size(); i++)
    log_value_ += Psi[i]->updateBuffer(P, buf, fromscratch);
  return log_value_;
}

void RPAJastrow::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  for (int i = 0; i < Psi.size(); i++)
    Psi[i]->copyFromBuffer(P, buf);
}

/** this is a great deal of logic for make clone I'm wondering what is going on
 */
std::unique_ptr<WaveFunctionComponent> RPAJastrow::makeClone(ParticleSet& tpq) const
{
  std::unique_ptr<HandlerType> tempHandler;
  if (rpafunc == "yukawa" || rpafunc == "breakup")
  {
    tempHandler =
        std::make_unique<LRHandlerTemp<YukawaBreakup<RealType>, LPQHIBasis>>(dynamic_cast<const LRHandlerTemp<
                                                                                 YukawaBreakup<RealType>, LPQHIBasis>&>(
                                                                                 *myHandler),
                                                                             tpq);
  }
  else if (rpafunc == "rpa")
  {
    tempHandler =
        std::make_unique<LRRPAHandlerTemp<RPABreakup<RealType>, LPQHIBasis>>(dynamic_cast<const LRRPAHandlerTemp<
                                                                                 RPABreakup<RealType>, LPQHIBasis>&>(
                                                                                 *myHandler),
                                                                             tpq);
  }
  else if (rpafunc == "dyukawa")
  {
    tempHandler = std::make_unique<LRHandlerTemp<
        DerivYukawaBreakup<RealType>,
        LPQHIBasis>>(dynamic_cast<const LRHandlerTemp<DerivYukawaBreakup<RealType>, LPQHIBasis>&>(*myHandler), tpq);
  }
  else if (rpafunc == "drpa")
  {
    tempHandler = std::make_unique<LRRPAHandlerTemp<
        DerivRPABreakup<RealType>,
        LPQHIBasis>>(dynamic_cast<const LRRPAHandlerTemp<DerivRPABreakup<RealType>, LPQHIBasis>&>(*myHandler), tpq);
  }

  auto myClone  = std::make_unique<RPAJastrow>(tpq);
  myClone->Rcut = Rcut;
  myClone->Kc   = Kc;
  myClone->setHandler(std::move(tempHandler));
  if (!DropLongRange)
    myClone->makeLongRange();
  if (!DropShortRange)
    myClone->makeShortRange();
  return myClone;
}

}; // namespace qmcplusplus
