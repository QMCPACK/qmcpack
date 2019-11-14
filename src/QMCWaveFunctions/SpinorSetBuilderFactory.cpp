//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/SpinorSetBuilderFactory.h"
#if OHMMS_DIM == 3
#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#else
#endif

#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/RotatedSPOs.h"
#endif

#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#endif
#endif
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/MPIObjectBase.h"


namespace qmcplusplus
{

/** constructor
 * \param els reference to the electrons
 * \param psi reference to the wavefunction
 * \param ions reference to the ions
 */
SpinorSetBuilderFactory::SpinorSetBuilderFactory(Communicate* comm, ParticleSet& els, PtclPoolType& psets)
    : SPOSetBuilderFactory(comm,els,psets)
{
  ClassName = "SpinorSetBuilderFactory";
}

SPOSetBuilder* SpinorSetBuilderFactory::createSPOSetBuilder(xmlNodePtr rootNode)
{
  ReportEngine PRE(ClassName, "createSPOSetBuilder");
  std::string sourceOpt("ion0");
  std::string type("");
  std::string name("");
  std::string keyOpt("NMO");       //gaussian Molecular Orbital
  std::string transformOpt("yes"); //numerical Molecular Orbital
  std::string cuspC("no");         // cusp correction
  std::string cuspInfo("");        // file with precalculated cusp correction info
  std::string MOH5Ref("");         // Path to H5 file for MO calculations
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sourceOpt, "source");
  aAttrib.add(cuspC, "cuspCorrection");
  aAttrib.add(type, "type");
  aAttrib.add(keyOpt, "keyword");
  aAttrib.add(keyOpt, "key");
  aAttrib.add(name, "name");
  aAttrib.add(transformOpt, "transform");
  aAttrib.add(cuspInfo, "cuspInfo");
  aAttrib.add(MOH5Ref, "href");

  if (rootNode != NULL)
    aAttrib.put(rootNode);

  std::string type_in = type;
  tolower(type);

  //when name is missing, type becomes the input
  if (name.empty())
    name = type_in;

  SPOSetBuilder* bb = 0;

  //check if builder can be reused
  std::map<std::string, SPOSetBuilder*>::iterator bbit = spo_builders.find(name);
  if (bbit != spo_builders.end())
  {
    app_log() << "Reuse SpinorSetBuilder \"" << name << "\" type " << type_in << std::endl;
    app_log().flush();
    bb                  = (*bbit).second;
    return last_builder = bb;
  }

  //assign last_builder
  bb = last_builder;

#if OHMMS_DIM == 3
  if (type.find("spline") < type.size())
  {
    name = type_in;
#if defined(HAVE_EINSPLINE)
    PRE << "EinsplineSetBuilder:  using libeinspline for B-spline spinors.\n";
    bb = new EinsplineSpinorSetBuilder(targetPtcl, ptclPool, myComm, rootNode);
#else
    PRE.error("Einspline is missing for B-spline orbitals", true);
#endif
  }
  else if (type == "molecularorbital" || type == "mo")
  {
    APP_ABORT("Spinors in LCAO type basis sets are not currently working.  Contact a developer\n");
//    ParticleSet* ions = 0;
    //initialize with the source tag
//    PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
//    if (pit == ptclPool.end())
//      PRE.error("Missing basisset/@source.", true);
//    else
//      ions = (*pit).second;
#if defined(ENABLE_SOA)
//    bb = new LCAOrbitalBuilder(targetPtcl, *ions, myComm, rootNode);
#else
#endif
  }
#endif //OHMMS_DIM==3
  PRE.flush();

  if (bb == 0)
    APP_ABORT_TRACE(__FILE__, __LINE__, "SpinorSetBuilderFactory::createSPOSetBuilder\n  SpinorSetBuilder creation failed.");

  if (bb == last_builder)
    app_log() << " Missing both \"@name\" and \"@type\". Use the last SpinorSetBuilder." << std::endl;
  else
  {
    app_log() << "  Created SpinorSetBuilder named '" << name << "' of type " << type << std::endl;
    spo_builders[name] = bb; //use name, if missing type is used
  }
  last_builder = bb;

  return bb;
}



} // namespace qmcplusplus
