//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "RadialJastrowBuilder.h"

namespace qmcplusplus
{

RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
					   ParticleSet& source):
  WaveFunctionComponentBuilder(target, psi),SourcePtcl(source)
{
  ClassName="RadialJastrowBuilder";
  NameOpt="0";
  typeOpt="unknown";
  Jastfunction="unknown";
  SourceOpt=targetPtcl.getName();
  SpinOpt="no";
}

RadialJastrowBuilder::RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi):
  WaveFunctionComponentBuilder(target, psi),SourcePtcl(NULL)
{
  ClassName="RadialJastrowBuilder";
  NameOpt="0";
  typeOpt="unknown";
  Jastfunction="unknown";
  SourceOpt=targetPtcl.getName();
  SpinOpt="no";
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

// quick template helper to allow use of RPA
template <typename> 
class RPAFunctor { };




bool RadialJastrowBuilder::put(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"put(xmlNodePtr)");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(NameOpt,"name");
  aAttrib.add(TypeOpt,"type");
  aAttrib.add(Jastfunction,"function");
  aAttrib.add(SourceOpt, "source");
  aAttrib.add(SpinOpt, "spin");
  aAttrib.put(cur);
  tolower(NameOpt);
  tolower(TypeOpt);
  tolower(Jastfunction);
  tolower(SourceOpt);
  tolower(SpinOpt);

  bool success=false;

  SpeciesSet& species(targetPtcl.getSpeciesSet());
  int chargeInd=species.addAttribute("charge");  

  if (typeOpt.find("one") < typeOpt.size())
  {
    // it's a one body jastrow factor
    if (Jastfunction == "bspline") 
    {
      success = createJ1<BsplineFunctor>(cur);
    }
    else if (Jastfunction == "pade") 
    {
      guardAgainstPBC();
      success = createJ1<PadeFunctor>(cur);
    }
    else if (Jastfunction == "rpa") 
    {
#if !(OHMMS_DIM == 3)
      app_error() << "RPA for one-body jastrow is only available for 3D\n";
#endif
      guardAgainstOBC();
      success = createJ1<RPAFunctor>(cur);
    }
    else
    {
      app_error() << "Unknown one jastrow body function: " << Jastfunction << ".\n";
    }
  }
  else if (typeOpt.find("two") < typeOpt.size())
  {
    // it's a two body jastrow factor
    if (Jastfunction == "bspline") 
    {
      success = createJ2<BsplineFunctor>(cur);
    }
    else if (Jastfunction == "pade") 
    {
      guardAgainstPBC();
      success = createJ2<PadeFunctor>(cur);
    }
    else if (Jastfunction == "rpa") 
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
}      
