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
#include <type_traits>
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/J1Spin.h"
#include "QMCWaveFunctions/Jastrow/J2OrbitalSoA.h"

#if defined(ENABLE_OFFLOAD)
#include "QMCWaveFunctions/Jastrow/J2OMPTarget.h"
#endif

#if defined(QMC_CUDA)
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbitalBspline.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbitalBspline.h"
#endif


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
// quick helper class to allow use of RPA
class RPAFunctor
{};

// helper class to simplify and localize ugly ifdef stuff for types
template<class RadFuncType, unsigned Implementation = RadialJastrowBuilder::detail::CPU>
class JastrowTypeHelper
{
public:
  using J1Type     = J1OrbitalSoA<RadFuncType>;
  using J1SpinType = J1Spin<RadFuncType>;
  using J2Type     = J2OrbitalSoA<RadFuncType>;
};

#if defined(QMC_CUDA)
template<>
class JastrowTypeHelper<BsplineFunctor<RadialJastrowBuilder::RealType>, RadialJastrowBuilder::detail::CUDA_LEGACY>
{
public:
  using RadFuncType = BsplineFunctor<RadialJastrowBuilder::RealType>;
  using J1Type      = OneBodyJastrowOrbitalBspline<RadFuncType>;
  using J1SpinType  = void;
  using J2Type      = TwoBodyJastrowOrbitalBspline<RadFuncType>;
};
#endif

#if defined(ENABLE_OFFLOAD)
template<>
class JastrowTypeHelper<BsplineFunctor<RadialJastrowBuilder::RealType>, RadialJastrowBuilder::detail::OMPTARGET>
{
public:
  using RadFuncType = BsplineFunctor<RadialJastrowBuilder::RealType>;
  using J2Type      = J2OMPTarget<RadFuncType>;
};
#endif

RadialJastrowBuilder::RadialJastrowBuilder(Communicate* comm, ParticleSet& target, ParticleSet& source)
    : WaveFunctionComponentBuilder(comm, target), SourcePtcl(&source)
{
  ClassName    = "RadialJastrowBuilder";
  NameOpt      = "0";
  TypeOpt      = "unknown";
  Jastfunction = "unknown";
  SpinOpt      = "no";
}

RadialJastrowBuilder::RadialJastrowBuilder(Communicate* comm, ParticleSet& target)
    : WaveFunctionComponentBuilder(comm, target), SourcePtcl(NULL)
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
  if (targetPtcl.getLattice().SuperCellEnum == SUPERCELL_OPEN)
  {
    app_error() << Jastfunction << " relies on the total density for its form\n";
    app_error() << "but open boundary conditions are requested.  Please choose other forms of Jastrow\n";
  }
}

// helper method for dealing with functor incompatible with PBC
void RadialJastrowBuilder::guardAgainstPBC()
{
  if (targetPtcl.getLattice().SuperCellEnum != SUPERCELL_OPEN)
  {
    app_error() << Jastfunction << " does not support a cutoff, but is requested with\n";
    app_error() << "periodic boundary conditions, please choose other forms of Jastrow\n";
  }
}

template<class RadFuncType>
void RadialJastrowBuilder::initTwoBodyFunctor(RadFuncType& functor, double fac)
{}

template<>
void RadialJastrowBuilder::initTwoBodyFunctor(BsplineFunctor<RealType>& bfunc, double fac)
{
  if (targetPtcl.getLattice().SuperCellEnum == SUPERCELL_OPEN) // for open systems, do nothing
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


template<class RadFuncType, unsigned Implementation>
std::unique_ptr<WaveFunctionComponent> RadialJastrowBuilder::createJ2(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ2(xmlNodePtr)");
  using Real       = typename RadFuncType::real_type;
  using J2Type     = typename JastrowTypeHelper<RadFuncType, Implementation>::J2Type;

  std::string input_name(getXMLAttributeValue(cur, "name"));
  std::string j2name = input_name.empty() ? "J2_" + Jastfunction : input_name;
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  auto J2  = std::make_unique<J2Type>(j2name, targetPtcl);

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
      int massInd   = species.addAttribute("mass");
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
        RealType qq       = species(chargeInd, ia) * species(chargeInd, ib);
        RealType red_mass = species(massInd, ia) * species(massInd, ib) / (species(massInd, ia) + species(massInd, ib));
#if OHMMS_DIM == 1
        RealType dim_factor = 1.0 / (OHMMS_DIM + 1);
#else
        RealType dim_factor = (ia == ib) ? 1.0 / (OHMMS_DIM + 1) : 1.0 / (OHMMS_DIM - 1);
#endif
        cusp = -2 * qq * red_mass * dim_factor;
      }
      app_summary() << "    Radial function for species: " << spA << " - " << spB << std::endl;
      app_debug() << "    RadialJastrowBuilder adds a functor with cusp = " << cusp << std::endl;

      auto functor = std::make_unique<RadFuncType>();
      functor->setCusp(cusp);
      functor->setPeriodic(targetPtcl.getLattice().SuperCellEnum != SUPERCELL_OPEN);
      functor->cutoff_radius   = targetPtcl.getLattice().WignerSeitzRadius;
      bool functor_initialized = functor->put(cur);
      if (!functor_initialized && init_mode == "rpa")
      {
        initTwoBodyFunctor(*functor, -cusp / 0.5);
      }

      app_summary() << std::endl;

      if (is_manager())
      {
        char fname[32];
        sprintf(fname, "J2.%s.%s.g%03d.dat", NameOpt.c_str(), pairType.c_str(), getGroupID());
        std::ofstream os(fname);
        print(*functor, os);
      }

      J2->addFunc(ia, ib, std::move(functor));
    }
    cur = cur->next;
  }
  J2->setOptimizable(true);

  // compute Chiesa Correction based on the current J2 parameters
  J2->ChiesaKEcorrection();

  // Ye: actually don't know what uk.dat is used for
  if (targetPtcl.getLattice().SuperCellEnum)
    computeJ2uk(J2->getPairFunctions());

  return J2;
}


template<class RadFuncType>
void RadialJastrowBuilder::computeJ2uk(const std::vector<RadFuncType*>& functors)
{
  const int numPoints = 1000;
  RealType vol        = targetPtcl.getLattice().Volume;
  int nsp             = targetPtcl.groups();
  FILE* fout          = 0;
  if (is_manager())
  {
    char fname[16];
    sprintf(fname, "uk.%s.g%03d.dat", NameOpt.c_str(), getGroupID());
    fout = fopen(fname, "w");
  }
  for (int iG = 0; iG < targetPtcl.getSimulationCell().getKLists().ksq.size(); iG++)
  {
    RealType Gmag = std::sqrt(targetPtcl.getSimulationCell().getKLists().ksq[iG]);
    RealType sum  = 0.0;
    RealType uk   = 0.0;
    for (int i = 0; i < targetPtcl.groups(); i++)
    {
      int Ni          = targetPtcl.last(i) - targetPtcl.first(i);
      RealType aparam = 0.0;
      for (int j = 0; j < targetPtcl.groups(); j++)
      {
        int Nj = targetPtcl.last(j) - targetPtcl.first(j);
        if (functors[i * nsp + j])
        {
          auto& ufunc     = *functors[i * nsp + j];
          RealType radius = ufunc.cutoff_radius;
          RealType k      = Gmag;
          RealType dr     = radius / (RealType)(numPoints - 1);
          for (int ir = 0; ir < numPoints; ir++)
          {
            RealType r = dr * (RealType)ir;
            RealType u = ufunc.evaluate(r);
            aparam += (1.0 / 4.0) * k * k * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr;
            uk += 0.5 * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr * (RealType)Nj / (RealType)(Ni + Nj);
          }
        }
      }
      //app_log() << "A = " << aparam << std::endl;
      sum += Ni * aparam / vol;
    }
    if (fout)
      fprintf(fout, "%1.8f %1.12e %1.12e\n", Gmag, uk, sum);
  }
  if (fout)
    fclose(fout);
}

// specialiation for J2 RPA jastrow.
template<>
std::unique_ptr<WaveFunctionComponent> RadialJastrowBuilder::createJ2<RPAFunctor>(xmlNodePtr cur)
{
  auto rpajastrow = std::make_unique<RPAJastrow>(targetPtcl);
  rpajastrow->put(cur);
  return rpajastrow;
}

template<class RadFuncType, bool SPIN, unsigned Implementation>
std::unique_ptr<WaveFunctionComponent> RadialJastrowBuilder::createJ1(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ1(xmlNodePtr)");
  using Real   = typename RadFuncType::real_type;
  using TH     = JastrowTypeHelper<RadFuncType, Implementation>;
  using J1Type = typename std::conditional<SPIN, typename TH::J1SpinType, typename TH::J1Type>::type;

  std::string input_name(getXMLAttributeValue(cur, "name"));
  std::string jname = input_name.empty() ? Jastfunction : input_name;

  auto J1 = std::make_unique<J1Type>(jname, *SourcePtcl, targetPtcl);

  xmlNodePtr kids = cur->xmlChildrenNode;

  // Find the number of the source species
  SpeciesSet& sSet = SourcePtcl->getSpeciesSet();
  SpeciesSet& tSet = targetPtcl.getSpeciesSet();
  bool success     = false;
  bool Opt(true);
  while (kids != NULL)
  {
    std::string kidsname(lowerCase(castXMLCharToChar(kids->name)));
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
      auto functor = std::make_unique<RadFuncType>();
      functor->setPeriodic(SourcePtcl->getLattice().SuperCellEnum != SUPERCELL_OPEN);
      functor->cutoff_radius = targetPtcl.getLattice().WignerSeitzRadius;
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
      app_summary() << "    Radial function for element: " << speciesA << " - "
                    << (speciesB.empty() ? targetPtcl.getName() : speciesB) << std::endl;
      functor->put(kids);
      app_summary() << std::endl;
      if (is_manager())
      {
        char fname[128];
        if (speciesB.size())
          sprintf(fname, "%s.%s.%s%s.g%03d.dat", jname.c_str(), NameOpt.c_str(), speciesA.c_str(), speciesB.c_str(),
                  getGroupID());
        else
          sprintf(fname, "%s.%s.%s.g%03d.dat", jname.c_str(), NameOpt.c_str(), speciesA.c_str(), getGroupID());
        std::ofstream os(fname);
        if (std::is_same<RadFuncType, PadeFunctor<RealType>>::value ||
            std::is_same<RadFuncType, Pade2ndOrderFunctor<RealType>>::value)
        {
          double plotextent = 10.0;
          print(*functor.get(), os, plotextent);
        }
        else
        {
          print(*functor.get(), os);
        }
      }
      J1->addFunc(ig, std::move(functor), jg);
      success = true;
    }
    kids = kids->next;
  }
  if (success)
  {
    J1->setOptimizable(Opt);
    return J1;
  }
  else
  {
    PRE.error("BsplineJastrowBuilder failed to add an One-Body Jastrow.");
    return std::unique_ptr<WaveFunctionComponent>();
  }
}

// specialiation for J1 RPA jastrow.  Note that the long range part is not implemented
template<>
std::unique_ptr<WaveFunctionComponent> RadialJastrowBuilder::createJ1<RPAFunctor>(xmlNodePtr cur)
{
  using Real             = RealType;
  using SplineEngineType = CubicBspline<Real, LINEAR_1DGRID, FIRSTDERIV_CONSTRAINTS>;
  using RadFunctorType   = CubicSplineSingle<Real, SplineEngineType>;
  using GridType         = LinearGrid<Real>;
  using HandlerType      = LRHandlerBase;
  using J1Type           = J1OrbitalSoA<RadFunctorType>;

  std::string input_name;
  std::string rpafunc = "RPA";
  OhmmsAttributeSet a;
  a.add(input_name, "name");
  a.add(rpafunc, "function");
  a.put(cur);
  ParameterSet params;
  RealType Rs(-1.0);
  RealType Kc(-1.0);
  params.add(Rs, "rs");
  params.add(Kc, "kc");
  params.put(cur);
  bool Opt(true);

  std::string jname = input_name.empty() ? Jastfunction : input_name;

  HandlerType* myHandler = nullptr;
  if (Rs < 0)
  {
    Rs = std::pow(3.0 / 4.0 / M_PI * targetPtcl.getLattice().Volume / static_cast<RealType>(targetPtcl.getTotalNum()),
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

  Real Rcut        = myHandler->get_rc() - 0.1;
  GridType* myGrid = new GridType;
  int npts         = static_cast<int>(Rcut / 0.01) + 1;
  myGrid->set(0, Rcut, npts);

  //create the numerical functor
  auto nfunc                       = std::make_unique<RadFunctorType>();
  ShortRangePartAdapter<Real>* SRA = new ShortRangePartAdapter<Real>(myHandler);
  SRA->setRmax(Rcut);
  nfunc->initialize(SRA, myGrid);

  auto J1 = std::make_unique<J1Type>(jname, *SourcePtcl, targetPtcl);

  SpeciesSet& sSet = SourcePtcl->getSpeciesSet();
  for (int ig = 0; ig < sSet.getTotalNum(); ig++)
  {
    J1->addFunc(ig, std::move(nfunc));
  }

  J1->setOptimizable(Opt);
  return J1;
}


std::unique_ptr<WaveFunctionComponent> RadialJastrowBuilder::buildComponent(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "put(xmlNodePtr)");
  std::string useGPU;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(NameOpt, "name");
  aAttrib.add(TypeOpt, "type");
  aAttrib.add(Jastfunction, "function");
  aAttrib.add(SpinOpt, "spin", {"no", "yes"});
#if defined(ENABLE_OFFLOAD)
  aAttrib.add(useGPU, "gpu", {"yes", "no"});
#endif
  aAttrib.put(cur);
  NameOpt = lowerCase(NameOpt);
  TypeOpt = lowerCase(TypeOpt);
  Jastfunction = lowerCase(Jastfunction);
  SpinOpt = lowerCase(SpinOpt);

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd = species.addAttribute("charge");

  if (TypeOpt.find("one") < TypeOpt.size())
  {
    // it's a one body jastrow factor
    if (Jastfunction == "bspline")
    {
      if (SpinOpt == "yes")
      {
#if defined(QMC_CUDA)
        myComm->barrier_and_abort("RadialJastrowBuilder::buildComponent spin resolved bspline Jastrow is not supported in legacy CUDA build.");
#else
        return createJ1<BsplineFunctor<RealType>, true>(cur);
#endif
      }
      else
      {
#if defined(QMC_CUDA)
        return createJ1<BsplineFunctor<RealType>, false, detail::CUDA_LEGACY>(cur);
#else
        return createJ1<BsplineFunctor<RealType>>(cur);
#endif
      }
    }
    else if (Jastfunction == "pade")
    {
      guardAgainstPBC();
      if (SpinOpt == "yes")
        return createJ1<PadeFunctor<RealType>, true>(cur);
      else
        return createJ1<PadeFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "pade2")
    {
      guardAgainstPBC();
      return createJ1<Pade2ndOrderFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "shortrangecusp")
    {
      //guardAgainstPBC(); // is this needed?
      if (SpinOpt == "yes")
        return createJ1<ShortRangeCuspFunctor<RealType>, true>(cur);
      else
        return createJ1<ShortRangeCuspFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "user")
    {
      if (SpinOpt == "yes")
        return createJ1<UserFunctor<RealType>, true>(cur);
      else
        return createJ1<UserFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "rpa")
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
      app_error() << "one body RPA jastrow is not supported at the moment\n";
      //return createJ1<RPAFunctor>(cur);
    }
    else
      app_error() << "Unknown one jastrow body function: " << Jastfunction << ".\n";
  }
  else if (TypeOpt.find("two") < TypeOpt.size())
  {
    // it's a two body jastrow factor
    if (Jastfunction == "bspline")
    {
#if defined(QMC_CUDA)
      return createJ2<BsplineFunctor<RealType>, detail::CUDA_LEGACY>(cur);
#else
#if defined(ENABLE_OFFLOAD)
      if (useGPU == "yes")
      {
        static_assert(std::is_same<JastrowTypeHelper<BsplineFunctor<RealType>, OMPTARGET>::J2Type,
                                   J2OMPTarget<BsplineFunctor<RealType>>>::value,
                      "check consistent type");
        if (targetPtcl.getCoordinates().getKind() != DynamicCoordinateKind::DC_POS_OFFLOAD)
        {
          std::ostringstream msg;
          msg << "Offload enabled Jastrow needs the gpu=\"yes\" attribute in the \"" << targetPtcl.getName()
              << "\" particleset" << std::endl;
          myComm->barrier_and_abort(msg.str());
        }
        app_summary() << "    Running on an accelerator via OpenMP offload." << std::endl;
        return createJ2<BsplineFunctor<RealType>, detail::OMPTARGET>(cur);
      }
      else
#endif
        return createJ2<BsplineFunctor<RealType>>(cur);
#endif
    }
    else if (Jastfunction == "pade")
    {
      guardAgainstPBC();
      return createJ2<PadeFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "user")
    {
      return createJ2<UserFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "rpa" || Jastfunction == "yukawa")
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#else
      guardAgainstOBC();
      return createJ2<RPAFunctor>(cur);
#endif
    }
    else
      app_error() << "Unknown two jastrow body function: " << Jastfunction << ".\n";
  }

  APP_ABORT("RadialJastrowBuilder::buildComponent not able to create Jastrow!\n");
  return nullptr;
}

} // namespace qmcplusplus
