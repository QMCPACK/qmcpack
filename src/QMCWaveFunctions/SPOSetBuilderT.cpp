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
#include "QMCWaveFunctions/RotatedSPOsT.h" // only for real wavefunctions

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


template<>
std::unique_ptr<SPOSetT<float>> SPOSetBuilderT<float>::createRotatedSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string method;
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.add(method, "method", {"global", "history"});
  attrib.put(cur);

  std::unique_ptr<SPOSetT<float>> sposet;
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
  auto rot_spo = std::make_unique<RotatedSPOsT<float>>(spo_object_name, std::move(sposet));

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
}

template<>
std::unique_ptr<SPOSetT<double>> SPOSetBuilderT<double>::createRotatedSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string method;
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.add(method, "method", {"global", "history"});
  attrib.put(cur);

  std::unique_ptr<SPOSetT<double>> sposet;
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
  auto rot_spo = std::make_unique<RotatedSPOsT<double>>(spo_object_name, std::move(sposet));

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
}

template<>
std::unique_ptr<SPOSetT<std::complex<float>>> SPOSetBuilderT<std::complex<float>>::createRotatedSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string method;
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.add(method, "method", {"global", "history"});
  attrib.put(cur);
  myComm->barrier_and_abort("Orbital optimization via rotation doesn't support complex wavefunctions yet.");
  return nullptr;
}

template<>
std::unique_ptr<SPOSetT<std::complex<double>>> SPOSetBuilderT<std::complex<double>>::createRotatedSPOSet(xmlNodePtr cur)
{
  std::string spo_object_name;
  std::string method;
  OhmmsAttributeSet attrib;
  attrib.add(spo_object_name, "name");
  attrib.add(method, "method", {"global", "history"});
  attrib.put(cur);
  myComm->barrier_and_abort("Orbital optimization via rotation doesn't support complex wavefunctions yet.");
  return nullptr;
}


template class SPOSetBuilderT<float>;
template class SPOSetBuilderT<double>;
template class SPOSetBuilderT<std::complex<float>>;
template class SPOSetBuilderT<std::complex<double>>;
} // namespace qmcplusplus
