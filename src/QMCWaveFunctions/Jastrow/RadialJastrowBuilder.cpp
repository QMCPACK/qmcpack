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
#include <qmc_common.h>


namespace qmcplusplus
{
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
WaveFunctionComponent* RadialJastrowBuilder::createJ2(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ2(xmlNodePtr)");
  using RT                = typename RadFuncType::real_type;
  using J2OrbitalType     = typename JastrowTypeHelper<RadFuncType>::J2OrbitalType;
  using DiffJ2OrbitalType = typename JastrowTypeHelper<RadFuncType>::DiffJ2OrbitalType;

  std::string j2name = "J2_" + Jastfunction;
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int taskid = is_manager() ? getGroupID() : -1;
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
  J2->setOptimizable(true);  

//  Chiesa KECorrection for SOA
   auto& F = J2->F;
   PtclRef = &targetPtcl;

  #if !defined(QMC_BUILD_SANDBOX_ONLY)
//  if ((!PtclRef->Lattice.SuperCellEnum))
//      return 0.0;
    const int numPoints = 1000;
    RealType vol        = PtclRef->Lattice.Volume;
    int nsp             = PtclRef->groups();
    //FILE *fout=(Write_Chiesa_Correction)?fopen ("uk.dat", "w"):0;
    FILE* fout = 0;
//    if (qmc_common.io_node && taskid > -1) //taskid=-1
    if (myComm->rank() == 0)
    {
      char fname[16];
      sprintf(fname, "uk.g%03d.dat", taskid);
      fout = fopen(fname, "w");
    }
    for (int iG = 0; iG < PtclRef->SK->KLists.ksq.size(); iG++)
    {
      RealType Gmag = std::sqrt(PtclRef->SK->KLists.ksq[iG]);
      RealType sum  = 0.0;
      RealType uk   = 0.0;
      for (int i = 0; i < PtclRef->groups(); i++)
      {
        int Ni = PtclRef->last(i) - PtclRef->first(i);
        RealType aparam = 0.0;
        for (int j = 0; j < PtclRef->groups(); j++)
        {
          int Nj = PtclRef->last(j) - PtclRef->first(j);
          if (F[i * nsp + j])
          {
            auto& ufunc = *(F[i * nsp + j]); 
            RealType radius = ufunc.cutoff_radius;
            RealType k      = Gmag;
            RealType dr     = radius / (RealType)(numPoints - 1);
            for (int ir = 0; ir < numPoints; ir++)
            {
              RealType r = dr * (RealType)ir;
              RealType u = ufunc.evaluate(r);
#if (OHMMS_DIM == 3)
              aparam += (1.0 / 4.0) * k * k * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr;
              uk     += 0.5 * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr * (RealType)Nj / (RealType)(Ni + Nj);
#endif
#if (OHMMS_DIM == 2)
              uk += 0.5 * 2.0 * M_PI * std::sin(k * r) / k * u * dr * (RealType)Nj / (RealType)(Ni + Nj);
#endif
              //aparam += 0.25* 4.0*M_PI*r*r*u*dr;
            }
          }
        }
        //app_log() << "A = " << aparam << std::endl;
        sum += Ni * aparam / vol;
      }
      if (iG == 0)
      {
        RealType a = 1.0;
        for (int iter = 0; iter < 20; iter++)
          a = uk / (4.0 * M_PI * (1.0 / (Gmag * Gmag) - 1.0 / (Gmag * Gmag + 1.0 / a)));
        KEcorr = 4.0 * M_PI * a / (4.0 * vol) * PtclRef->getTotalNum();
      }
      if (fout)
        fprintf(fout, "%1.8f %1.12e %1.12e\n", Gmag, uk, sum);
    }
    if (fout)
      fclose(fout);
#endif

    J2->KEcorr = KEcorr;

  return J2;
}

// specialiation for J2 RPA jastrow.
template<>
WaveFunctionComponent* RadialJastrowBuilder::createJ2<RPAFunctor>(xmlNodePtr cur)
{
  RPAJastrow* rpajastrow = new RPAJastrow(targetPtcl, is_manager());
  rpajastrow->put(cur);
  return rpajastrow;
}

template<class RadFuncType>
WaveFunctionComponent* RadialJastrowBuilder::createJ1(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "createJ1(xmlNodePtr)");
  using RT                = typename RadFuncType::real_type;
  using J1OrbitalType     = typename JastrowTypeHelper<RadFuncType>::J1OrbitalType;
  using DiffJ1OrbitalType = typename JastrowTypeHelper<RadFuncType>::DiffJ1OrbitalType;

  int taskid             = getGroupID();
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
    J1->setOptimizable(Opt);
    return J1;
  }
  else
  {
    PRE.error("BsplineJastrowBuilder failed to add an One-Body Jastrow.");
    delete J1;
    delete dJ1;
    return nullptr;
  }
}

// specialiation for J1 RPA jastrow.  Note that the long range part is not implemented
template<>
WaveFunctionComponent* RadialJastrowBuilder::createJ1<RPAFunctor>(xmlNodePtr cur)
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

  HandlerType* myHandler = nullptr;
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
  J1->setOptimizable(Opt);
  return J1;
}


WaveFunctionComponent* RadialJastrowBuilder::buildComponent(xmlNodePtr cur)
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

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd = species.addAttribute("charge");

  if (TypeOpt.find("one") < TypeOpt.size())
  {
    // it's a one body jastrow factor
    if (Jastfunction == "bspline")
    {
      return createJ1<BsplineFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "pade")
    {
      guardAgainstPBC();
      return createJ1<PadeFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "shortrangecusp")
    {
      //guardAgainstPBC(); // is this needed?
      return createJ1<ShortRangeCuspFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "user")
    {
      return createJ1<UserFunctor<RealType>>(cur);
    }
    else if (Jastfunction == "rpa")
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
#if defined(ENABLE_SOA)
      app_error() << "one body RPA jastrow is not compatible with SOA at the moment\n";
#else
      return createJ1<RPAFunctor>(cur);
#endif
    }
    else
      app_error() << "Unknown one jastrow body function: " << Jastfunction << ".\n";
  }
  else if (TypeOpt.find("two") < TypeOpt.size())
  {
    // it's a two body jastrow factor
    if (Jastfunction == "bspline")
    {
      return createJ2<BsplineFunctor<RealType>>(cur);
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
