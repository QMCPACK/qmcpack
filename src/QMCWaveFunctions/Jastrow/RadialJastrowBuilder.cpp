//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "RadialJastrowBuilder.h"

#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"
#else
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#endif

#if defined(QMC_CUDA)
#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBspline.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBspline.h"
#else
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBsplineAoS.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBsplineAoS.h"
#endif
#endif

#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodySpinJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodySpinJastrowOrbital.h"

#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/PadeFunctors.h"
#include "QMCWaveFunctions/Jastrow/ShortRangeCuspFunctor.h"
#include "QMCWaveFunctions/Jastrow/UserFunctor.h"
#include <iostream>


namespace qmcplusplus
{
RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, ParticleSet& source)
    : WaveFunctionComponentBuilder(target, psi), SourcePtcl(&source)
{
  ClassName    = "RadialJastrowBuilder";
  NameOpt      = "0";
  TypeOpt      = "unknown";
  Jastfunction = "unknown";
  SpinOpt      = "no";
}

RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi)
    : WaveFunctionComponentBuilder(target, psi), SourcePtcl(NULL)
{
  ClassName    = "RadialJastrowBuilder";
  NameOpt      = "0";
  TypeOpt      = "unknown";
  Jastfunction = "unknown";
  SpinOpt      = "no";
}

// helper method for dealing with functor incompatible with Open Boundaries
void RadialJastrowBuilder::guardAgainstOBC()
{
  if (targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
  {
    app_error() << Jastfunction << " relies on the total density for its form\n";
    app_error() << "but open boundary conditions are requested.  Please choose other forms of Jastrow\n";
  }
}

// helper method for dealing with functor incompatible with PBC
void RadialJastrowBuilder::guardAgainstPBC()
{
  if (targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN)
  {
    app_error() << Jastfunction << " does not support a cutoff, but is requested with\n";
    app_error() << "periodic boundary conditions, please choose other forms of Jastrow\n";
  }
}

// quick helper class to allow use of RPA
class RPAFunctor
{};

// helper class to simplify and localize ugly ifdef stuff for types
template<class RadFuncType>
class JastrowTypeHelper
{
public:
#if defined(ENABLE_SOA)
  using J1OrbitalType = J1OrbitalSoA<RadFuncType>;
  using J2OrbitalType = J2OrbitalSoA<RadFuncType>;
#else
  using J1OrbitalType = OneBodyJastrowOrbital<RadFuncType>;
  using J2OrbitalType = TwoBodyJastrowOrbital<RadFuncType>;
#endif
  using DiffJ1OrbitalType = DiffOneBodyJastrowOrbital<RadFuncType>;
  using DiffJ2OrbitalType = DiffTwoBodyJastrowOrbital<RadFuncType>;
  // spin polarized J1
  using J1OrbitalTypeSpin     = OneBodySpinJastrowOrbital<RadFuncType>;
  using DiffJ1OrbitalTypeSpin = DiffOneBodyJastrowOrbital<RadFuncType>;
};

template<>
class JastrowTypeHelper<BsplineFunctor<RadialJastrowBuilder::RealType>>
{
public:
  using RadFuncType = BsplineFunctor<RadialJastrowBuilder::RealType>;
#if defined(QMC_CUDA) && defined(ENABLE_SOA)
  using J1OrbitalType = OneBodyJastrowOrbitalBspline<RadFuncType>;
  using J2OrbitalType = TwoBodyJastrowOrbitalBspline<RadFuncType>;
#endif
#if defined(QMC_CUDA) && !defined(ENABLE_SOA)
  using J1OrbitalType = OneBodyJastrowOrbitalBsplineAoS;
  using J2OrbitalType = TwoBodyJastrowOrbitalBsplineAoS;
#endif
#if !defined(QMC_CUDA) && defined(ENABLE_SOA)
  using J1OrbitalType = J1OrbitalSoA<RadFuncType>;
  using J2OrbitalType = J2OrbitalSoA<RadFuncType>;
#endif
#if !defined(QMC_CUDA) && !defined(ENABLE_SOA)
  using J1OrbitalType = OneBodyJastrowOrbital<RadFuncType>;
  using J2OrbitalType = TwoBodyJastrowOrbital<RadFuncType>;
#endif
  using DiffJ1OrbitalType = DiffOneBodyJastrowOrbital<RadFuncType>;
  using DiffJ2OrbitalType = DiffTwoBodyJastrowOrbital<RadFuncType>;
  // spin polarized J1
  using J1OrbitalTypeSpin     = OneBodySpinJastrowOrbital<RadFuncType>;
  using DiffJ1OrbitalTypeSpin = DiffOneBodyJastrowOrbital<RadFuncType>;
};

template<class RadFuncType>
void RadialJastrowBuilder::initTwoBodyFunctor(RadFuncType& functor, double fac)
{}

template<>
void RadialJastrowBuilder::initTwoBodyFunctor(BsplineFunctor<RealType>& bfunc, double fac)
{
  if (targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN) // for open systems, do nothing
  {
    return;
  }
  app_log() << "  Initializing Two-Body with RPA Jastrow " << std::endl;
  std::vector<RealType> rpaValues;
  int npts = bfunc.NumParams;
  if (rpaValues.empty())
  {
    rpaValues.resize(npts);
    LRRPAHandlerTemp<RPABreakup<RealType>, LPQHIBasis> rpa(targetPtcl, -1.0);
    rpa.Breakup(targetPtcl, -1.0);
    RealType dr = bfunc.cutoff_radius / static_cast<RealType>(npts);
    RealType r  = 0;
    for (int i = 0; i < npts; i++)
    {
      rpaValues[i] = rpa.evaluate(r, 1.0 / r); //y[i]=fac*rpa.evaluate(r,1.0/r);
      r += dr;
    }
  }
  RealType last = rpaValues[npts - 1];

  for (int i = 0; i < npts; i++)
    bfunc.Parameters[i] = fac * (rpaValues[i] - last);
  bfunc.reset();
}


template<class RadFuncType>
bool RadialJastrowBuilder::createJ2(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ2(xmlNodePtr)");
  using RT                = typename RadFuncType::real_type;
  using J2OrbitalType     = typename JastrowTypeHelper<RadFuncType>::J2OrbitalType;
  using DiffJ2OrbitalType = typename JastrowTypeHelper<RadFuncType>::DiffJ2OrbitalType;

  std::string j2name = "J2_" + Jastfunction;
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int taskid = (targetPsi.is_manager()) ? targetPsi.getGroupID() : -1;
  auto* J2   = new J2OrbitalType(targetPtcl, taskid);
  auto* dJ2  = new DiffJ2OrbitalType(targetPtcl);

  std::string init_mode("0");
  {
    OhmmsAttributeSet hAttrib;
    hAttrib.add(init_mode, "init");
    hAttrib.put(cur);
  }

  cur = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "correlation")
    {
      OhmmsAttributeSet rAttrib;
      RealType cusp = -1e10;
      std::string spA(species.speciesName[0]);
      std::string spB(species.speciesName[0]);
      std::string pairType("0");
      rAttrib.add(spA, "speciesA");
      rAttrib.add(spB, "speciesB");
      rAttrib.add(pairType, "pairType");
      rAttrib.add(cusp, "cusp");

      rAttrib.put(cur);
      if (pairType[0] == '0')
      {
        pairType = spA + spB;
      }
      else
      {
        PRE.warning("pairType is deprecated. Use speciesA/speciesB");
        //overwrite the species
        spA = pairType[0];
        spB = pairType[1];
      }

      int ia        = species.findSpecies(spA);
      int ib        = species.findSpecies(spB);
      int chargeInd = species.addAttribute("charge");
      std::string illegal_species;
      if (ia == species.size())
        illegal_species = spA;
      if (ib == species.size())
      {
        if (illegal_species.size())
          illegal_species += " and ";
        illegal_species += spB;
      }
      if (illegal_species.size())
        PRE.error("species " + illegal_species + " requested for Jastrow " + j2name +
                      " does not exist in ParticleSet " + targetPtcl.getName(),
                  true);
      if (ia == ib && (targetPtcl.last(ia) - targetPtcl.first(ia) == 1))
        PRE.error("Failed to add " + spA + spB + " correlation for only 1 " + spA +
                      " particle. Please remove it from two-body Jastrow.",
                  true);
      if (cusp < -1e6)
      {
        RealType qq = species(chargeInd, ia) * species(chargeInd, ib);
        cusp        = (ia == ib) ? -0.25 * qq : -0.5 * qq;
      }
      app_summary() << "    Radial function for species: " << spA << " - " << spB << std::endl;
      app_debug() << "    RadialJastrowBuilder adds a functor with cusp = " << cusp << std::endl;

      auto* functor = new RadFuncType();
      functor->setCusp(cusp);
      functor->setPeriodic(targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN);
      functor->cutoff_radius   = targetPtcl.Lattice.WignerSeitzRadius;
      bool functor_initialized = functor->put(cur);
      if (!functor_initialized && init_mode == "rpa")
      {
        initTwoBodyFunctor(*functor, -cusp / 0.5);
      }

      app_summary() << std::endl;

      J2->addFunc(ia, ib, functor);
      dJ2->addFunc(ia, ib, functor);

      if (qmc_common.io_node)
      {
        char fname[32];
        if (qmc_common.mpi_groups > 1)
          sprintf(fname, "J2.%s.g%03d.dat", pairType.c_str(), taskid);
        else
          sprintf(fname, "J2.%s.dat", pairType.c_str());
        std::ofstream os(fname);
        print(*functor, os);
      }
    }
    cur = cur->next;
  }
  J2->dPsi = dJ2;
  targetPsi.addComponent(J2, j2name.c_str());
  J2->setOptimizable(true);
  return true;
}

// specialiation for J2 RPA jastrow.
template<>
bool RadialJastrowBuilder::createJ2<RPAFunctor>(xmlNodePtr cur)
{
  RPAJastrow* rpajastrow = new RPAJastrow(targetPtcl, targetPsi.is_manager());
  rpajastrow->put(cur);
  targetPsi.addComponent(rpajastrow, NameOpt);
  return true;
}

template<class RadFuncType>
bool RadialJastrowBuilder::createJ1(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ1(xmlNodePtr)");
  using RT                = typename RadFuncType::real_type;
  using J1OrbitalType     = typename JastrowTypeHelper<RadFuncType>::J1OrbitalType;
  using DiffJ1OrbitalType = typename JastrowTypeHelper<RadFuncType>::DiffJ1OrbitalType;

  int taskid             = targetPsi.getGroupID();
  J1OrbitalType* J1      = new J1OrbitalType(*SourcePtcl, targetPtcl);
  DiffJ1OrbitalType* dJ1 = new DiffJ1OrbitalType(*SourcePtcl, targetPtcl);

  xmlNodePtr kids = cur->xmlChildrenNode;

  // Find the number of the source species
  SpeciesSet& sSet = SourcePtcl->getSpeciesSet();
  SpeciesSet& tSet = targetPtcl.getSpeciesSet();
  bool success     = false;
  bool Opt(true);
  std::string jname = "J1_" + Jastfunction;
  while (kids != NULL)
  {
    std::string kidsname = (char*)kids->name;
    tolower(kidsname);
    if (kidsname == "correlation")
    {
      std::string speciesA;
      std::string speciesB;
      RealType cusp(0);
      OhmmsAttributeSet rAttrib;
      rAttrib.add(speciesA, "elementType");
      rAttrib.add(speciesA, "speciesA");
      rAttrib.add(speciesB, "speciesB");
      rAttrib.add(cusp, "cusp");
      rAttrib.put(kids);
      auto* functor = new RadFuncType();
      functor->setPeriodic(SourcePtcl->Lattice.SuperCellEnum != SUPERCELL_OPEN);
      functor->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
      functor->setCusp(cusp);
      const int ig = sSet.findSpecies(speciesA);
      const int jg = speciesB.size() ? tSet.findSpecies(speciesB) : -1;
      if (ig == sSet.getTotalNum())
      {
        PRE.error("species " + speciesA + " requested for Jastrow " + jname + " does not exist in ParticleSet " +
                      SourcePtcl->getName(),
                  true);
      }
      if (jg == tSet.getTotalNum())
      {
        PRE.error("species " + speciesB + " requested for Jastrow " + jname + " does not exist in ParticleSet " +
                      targetPtcl.getName(),
                  true);
      }
      app_summary() << "    Radial function for element: " << speciesA << std::endl;
      functor->put(kids);
      app_summary() << std::endl;
      J1->addFunc(ig, functor, jg);
      dJ1->addFunc(ig, functor, jg);
      success = true;
      if (qmc_common.io_node)
      {
        char fname[128];
        if (qmc_common.mpi_groups > 1)
        {
          if (speciesB.size())
            sprintf(fname, "%s.%s%s.g%03d.dat", jname.c_str(), speciesA.c_str(), speciesB.c_str(), taskid);
          else
            sprintf(fname, "%s.%s.g%03d.dat", jname.c_str(), speciesA.c_str(), taskid);
        }
        else
        {
          if (speciesB.size())
            sprintf(fname, "%s.%s%s.dat", jname.c_str(), speciesA.c_str(), speciesB.c_str());
          else
            sprintf(fname, "%s.%s.dat", jname.c_str(), speciesA.c_str());
        }
        std::ofstream os(fname);
        print(*functor, os);
      }
    }
    kids = kids->next;
  }
  if (success)
  {
    J1->dPsi = dJ1;
    targetPsi.addComponent(J1, jname.c_str());
    J1->setOptimizable(Opt);
    return true;
  }
  else
  {
    PRE.warning("BsplineJastrowBuilder failed to add an One-Body Jastrow.");
    delete J1;
    delete dJ1;
    return false;
  }
}

// specialiation for J1 RPA jastrow.  Note that the long range part is not implemented
template<>
bool RadialJastrowBuilder::createJ1<RPAFunctor>(xmlNodePtr cur)
{
  using RT               = RealType;
  using SplineEngineType = CubicBspline<RT, LINEAR_1DGRID, FIRSTDERIV_CONSTRAINTS>;
  using RadFunctorType   = CubicSplineSingle<RT, SplineEngineType>;
  using GridType         = LinearGrid<RT>;
  using HandlerType      = LRHandlerBase;
#if defined(ENABLE_SOA)
  using J1OrbitalType = J1OrbitalSoA<RadFunctorType>;
#endif
#if !defined(ENABLE_SOA)
  using J1OrbitalType = OneBodyJastrowOrbital<RadFunctorType>;
#endif
  using DiffJ1OrbitalType = DiffOneBodyJastrowOrbital<RadFunctorType>;

  /*
  using J1OrbitalType = JastrowTypeHelper<RT, RadFunctorType>::J1OrbitalType;
  using DiffJ1OrbitalType = JastrowTypeHelper<RT, RadFunctorType>::DiffJ1OrbitalType;
  */

  std::string MyName  = "Jep";
  std::string rpafunc = "RPA";
  OhmmsAttributeSet a;
  a.add(MyName, "name");
  a.add(rpafunc, "function");
  a.put(cur);
  ParameterSet params;
  RealType Rs(-1.0);
  RealType Kc(-1.0);
  params.add(Rs, "rs", "double");
  params.add(Kc, "kc", "double");
  params.put(cur);
  bool Opt(true);

  HandlerType* myHandler;
  if (Rs < 0)
  {
    Rs = std::pow(3.0 / 4.0 / M_PI * targetPtcl.Lattice.Volume / static_cast<RealType>(targetPtcl.getTotalNum()),
                  1.0 / 3.0);
  }
  if (Kc < 0)
  {
    Kc = 1e-6;
  }
  if (rpafunc == "RPA")
  {
    myHandler = new LRRPAHandlerTemp<EPRPABreakup<RealType>, LPQHIBasis>(targetPtcl, Kc);
    app_log() << "  using e-p RPA" << std::endl;
  }
  else if (rpafunc == "dRPA")
  {
    myHandler = new LRRPAHandlerTemp<derivEPRPABreakup<RealType>, LPQHIBasis>(targetPtcl, Kc);
    app_log() << "  using e-p derivRPA" << std::endl;
  }
  myHandler->Breakup(targetPtcl, Rs);

  RT Rcut          = myHandler->get_rc() - 0.1;
  GridType* myGrid = new GridType;
  int npts         = static_cast<int>(Rcut / 0.01) + 1;
  myGrid->set(0, Rcut, npts);

  //create the numerical functor
  RadFunctorType* nfunc          = new RadFunctorType;
  ShortRangePartAdapter<RT>* SRA = new ShortRangePartAdapter<RT>(myHandler);
  SRA->setRmax(Rcut);
  nfunc->initialize(SRA, myGrid);

  J1OrbitalType* J1      = new J1OrbitalType(*SourcePtcl, targetPtcl);
  DiffJ1OrbitalType* dJ1 = new DiffJ1OrbitalType(*SourcePtcl, targetPtcl);

  SpeciesSet& sSet = SourcePtcl->getSpeciesSet();
  for (int ig = 0; ig < sSet.getTotalNum(); ig++)
  {
    J1->addFunc(ig, nfunc);
    dJ1->addFunc(ig, nfunc);
  }

  J1->dPsi          = dJ1;
  std::string jname = "J1_" + Jastfunction;
  targetPsi.addComponent(J1, jname.c_str());
  J1->setOptimizable(Opt);
  return true;
}


bool RadialJastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "put(xmlNodePtr)");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(NameOpt, "name");
  aAttrib.add(TypeOpt, "type");
  aAttrib.add(Jastfunction, "function");
  aAttrib.add(SpinOpt, "spin");
  aAttrib.put(cur);
  tolower(NameOpt);
  tolower(TypeOpt);
  tolower(Jastfunction);
  tolower(SpinOpt);

  bool success = false;

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd = species.addAttribute("charge");

  if (TypeOpt.find("one") < TypeOpt.size())
  {
    // it's a one body jastrow factor
    if (Jastfunction == "bspline")
    {
      success = createJ1<BsplineFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "pade")
    {
      guardAgainstPBC();
      success = createJ1<PadeFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "shortrangecusp")
    {
      //guardAgainstPBC(); // is this needed?
      success = createJ1<ShortRangeCuspFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "user")
    {
      success = createJ1<UserFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "rpa")
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
#if defined(ENABLE_SOA)
      app_error() << "one body RPA jastrow is not compatible with SOA at the moment\n";
      success = false;
#else
      success = createJ1<RPAFunctor>(cur);
#endif
    }
    else
    {
      app_error() << "Unknown one jastrow body function: " << Jastfunction << ".\n";
    }
  }
  else if (TypeOpt.find("two") < TypeOpt.size())
  {
    // it's a two body jastrow factor
    if (Jastfunction == "bspline")
    {
      success = createJ2<BsplineFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "pade")
    {
      guardAgainstPBC();
      success = createJ2<PadeFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "user")
    {
      success = createJ2<UserFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "rpa" || Jastfunction == "yukawa")
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
      success = createJ2<RPAFunctor>(cur);
    }
    else
    {
      app_error() << "Unknown two body jastrow function: " << Jastfunction << ".\n";
    }
  }
  return success;
}

} // namespace qmcplusplus
