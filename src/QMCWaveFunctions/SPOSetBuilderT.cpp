//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SPOSetBuilderT.h"
#include "OhmmsData/AttributeSet.h"
#include <Message/UniformCommunicateError.h>

#ifndef QMC_COMPLEX
#include "QMCWaveFunctions/RotatedSPOsT.h"
#endif

namespace qmcplusplus
{
template<typename T>
SPOSetBuilderT<T>::SPOSetBuilderT(const std::string& type_name, Communicate* comm)
    : MPIObjectBase(comm), legacy(true), type_name_(type_name)
{
  reserve_states();
}

template<typename T>
void SPOSetBuilderT<T>::reserve_states(int nsets)
{
  int sets_needed = nsets - states.size();
  if (sets_needed > 0)
    for (int s = 0; s < sets_needed; ++s)
      states.push_back(std::make_unique<SPOSetInfo>());
}

template<typename T>
std::unique_ptr<SPOSetT<T>> SPOSetBuilderT<T>::createSPOSet(xmlNodePtr cur, SPOSetInputInfo& input_info)
{
  myComm->barrier_and_abort("BasisSetBase::createSPOSet(cur,input_info) has not been implemented");
  return 0;
}

template<typename T>
std::unique_ptr<SPOSetT<T>> SPOSetBuilderT<T>::createSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string optimize("no");

  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "id");
  attrib.add(spo_object_name, "name");
  attrib.add(optimize, "optimize");
  attrib.put(cur);

  app_summary() << std::endl;
  app_summary() << "     Single particle orbitals (SPO)" << std::endl;
  app_summary() << "     ------------------------------" << std::endl;
  app_summary() << "      Name: " << spo_object_name << "   Type: " << type_name_
                << "   Builder class name: " << ClassName << std::endl;
  app_summary() << std::endl;

  if (spo_object_name.empty())
    myComm->barrier_and_abort("SPOSet object \"name\" attribute not given in the input!");

  // read specialized sposet construction requests
  //   and translate them into a set of orbital indices
  SPOSetInputInfo input_info(cur);

  // process general sposet construction requests
  //   and preserve legacy interface
  std::unique_ptr<SPOSetT<T>> sposet;

  try
  {
    if (legacy && input_info.legacy_request)
      sposet = createSPOSetFromXML(cur);
    else
      sposet = createSPOSet(cur, input_info);
  }
  catch (const UniformCommunicateError& ue)
  {
    myComm->barrier_and_abort(ue.what());
  }

  if (!sposet)
    myComm->barrier_and_abort("SPOSetBuilderT::createSPOSet sposet creation failed");

  if (optimize == "rotation" || optimize == "yes")
  {
#ifdef QMC_COMPLEX
    app_error() << "Orbital optimization via rotation doesn't support complex wavefunction yet.\n";
    abort();
#else
    app_warning() << "Specifying orbital rotation via optimize tag is deprecated. Use the rotated_spo element instead"
                  << std::endl;

    sposet->storeParamsBeforeRotation();
    // create sposet with rotation
    auto& sposet_ref = *sposet;
    app_log() << "  SPOSet " << sposet_ref.getName() << " is optimizable\n";
    if (!sposet_ref.isRotationSupported())
      myComm->barrier_and_abort("Orbital rotation not supported with '" + sposet_ref.getName() + "' of type '" +
                                sposet_ref.getClassName() + "'.");
    auto rot_spo    = std::make_unique<RotatedSPOsT<T>>(sposet_ref.getName(), std::move(sposet));
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while (tcur != NULL)
    {
      std::string cname((const char*)(tcur->name));
      if (cname == "opt_vars")
      {
        std::vector<RealType> params;
        putContent(params, tcur);
        rot_spo->setRotationParameters(params);
      }
      tcur = tcur->next;
    }
    sposet = std::move(rot_spo);
#endif
  }

  if (sposet->getName().empty())
    app_warning() << "SPOSet object doesn't have a name." << std::endl;
  if (!spo_object_name.empty() && sposet->getName() != spo_object_name)
    app_warning() << "SPOSet object name mismatched! input name: " << spo_object_name
                  << "   object name: " << sposet->getName() << std::endl;

  sposet->checkObject();
  return sposet;
}

template<typename T>
std::unique_ptr<SPOSetT<T>> SPOSetBuilderT<T>::createRotatedSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string method;
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.add(method, "method", {"global", "history"});
  attrib.put(cur);


#ifdef QMC_COMPLEX
  myComm->barrier_and_abort("Orbital optimization via rotation doesn't support complex wavefunctions yet.");
  return nullptr;
#else
  std::unique_ptr<SPOSetT<T>> sposet;
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "sposet")
    {
      sposet = createSPOSet(element);
    }
  });

  if (!sposet)
    myComm->barrier_and_abort("Rotated SPO needs an SPOset");

  if (!sposet->isRotationSupported())
    myComm->barrier_and_abort("Orbital rotation not supported with '" + sposet->getName() + "' of type '" +
                              sposet->getClassName() + "'.");

  sposet->storeParamsBeforeRotation();
  auto rot_spo = std::make_unique<RotatedSPOsT<T>>(spo_object_name, std::move(sposet));

  if (method == "history")
    rot_spo->set_use_global_rotation(false);

  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "opt_vars")
    {
      std::vector<RealType> params;
      putContent(params, element);
      rot_spo->setRotationParameters(params);
    }
  });
  return rot_spo;
#endif
}
template class SPOSetBuilderT<double>;
template class SPOSetBuilderT<float>;
template class SPOSetBuilderT<std::complex<double>>;
template class SPOSetBuilderT<std::complex<float>>;
} // namespace qmcplusplus
