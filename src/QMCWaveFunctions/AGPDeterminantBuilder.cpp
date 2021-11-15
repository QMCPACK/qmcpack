//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 *@brief definition of three-body jastrow of Geminal functions
 */
#include "AGPDeterminantBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"

namespace qmcplusplus
{
AGPDeterminantBuilder::AGPDeterminantBuilder(Communicate* comm, ParticleSet& els, PtclPoolType& pset)
    : WaveFunctionComponentBuilder(comm, els), ptclPool(pset)
{}

template<class BasisBuilderT>
bool AGPDeterminantBuilder::createAGP(BasisBuilderT* abuilder, xmlNodePtr cur)
{
  bool spinpolarized                             = false;
  typename BasisBuilderT::BasisSetType* basisSet = 0;
  cur                                            = cur->xmlChildrenNode;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "coefficient" || cname == "coefficients")
    {
      if (agpDet == nullptr)
      {
        int nup = targetPtcl.first(1), ndown = 0;
        if (targetPtcl.groups() > 1)
          ndown = targetPtcl.first(2) - nup;
        basisSet->resize(nup + ndown);
        agpDet = std::make_unique<AGPDeterminant>(basisSet);
        agpDet->resize(nup, ndown);
      }
      int offset      = 1;
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while (tcur != NULL)
      {
        if (xmlStrEqual(tcur->name, (const xmlChar*)"lambda"))
        {
          int i                                  = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"i")));
          int j                                  = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"j")));
          double c                               = atof((const char*)(xmlGetProp(tcur, (const xmlChar*)"c")));
          agpDet->Lambda(i - offset, j - offset) = c;
          if (i != j)
          {
            agpDet->Lambda(j - offset, i - offset) = c;
          }
        }
        tcur = tcur->next;
      }
    }
    else if (cname == "unpaired")
    {
      spinpolarized   = true;
      int offset      = 1;
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while (tcur != NULL)
      {
        if (xmlStrEqual(tcur->name, (const xmlChar*)"lambda"))
        {
          int i                                    = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"i")));
          int j                                    = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"j")));
          double c                                 = atof((const char*)(xmlGetProp(tcur, (const xmlChar*)"c")));
          agpDet->LambdaUP(j - offset, i - offset) = c;
        }
        tcur = tcur->next;
      }
    }
    cur = cur->next;
  }
  //app_log() << agpDet->Lambda << std::endl;
  if (spinpolarized)
  {
    app_log() << "  Coefficients for the unpaired electrons " << std::endl;
    app_log() << agpDet->LambdaUP << std::endl;
  }
  return true;
}

std::unique_ptr<WaveFunctionComponent> AGPDeterminantBuilder::buildComponent(xmlNodePtr cur)
{
  if (agpDet)
  {
    APP_ABORT("  AGPDeterminantBuilder::put exits. AGPDeterminant has been already created.\n");
    return std::unique_ptr<AGPDeterminant>();
  }
  app_log() << "  AGPDeterminantBuilder Creating AGPDeterminant." << std::endl;
  xmlNodePtr curRoot = cur;
  bool success       = true;
  std::string cname, tname;
  xmlNodePtr cPtr = NULL;
  xmlNodePtr uPtr = NULL;
  OhmmsAttributeSet oAttrib;
  oAttrib.add(funcOpt, "function");
  oAttrib.add(transformOpt, "transform");
  oAttrib.put(cur);
  cur = cur->children;
  while (cur != NULL)
  {
    getNodeName(cname, cur);
    if (cname.find("coeff") < cname.size())
    {
      cPtr = cur;
    }
    else if (cname.find("un") < cname.size())
    {
      uPtr = cur;
    }
    cur = cur->next;
  }
  if (cPtr == NULL)
  {
    std::ostringstream err_msg;
    err_msg << "  AGPDeterminantBuilder::put Cannot create AGPDeterminant." << std::endl
            << "    Missing <coefficients/>" << std::endl;
    APP_ABORT(err_msg.str());
    return nullptr;
  }
  if (!mySPOSetBuilderFactory)
  {
    mySPOSetBuilderFactory = std::make_unique<SPOSetBuilderFactory>(myComm, targetPtcl, ptclPool);
    mySPOSetBuilderFactory->createSPOSetBuilder(curRoot);
  }
  // mmorales: this needs to be fixed after changes to BasisSetfactory
  //    BasisSetBase<RealType>* myBasisSet=mySPOSetBuilderFactory->getBasisSet();
  BasisSetBase<RealType>* myBasisSet = 0; //=mySPOSetBuilderFactory->getBasisSet();
  int nup = targetPtcl.first(1), ndown = 0;
  if (targetPtcl.groups() > 1)
    ndown = targetPtcl.first(2) - nup;
  myBasisSet->resize(nup + ndown);
  agpDet = std::make_unique<AGPDeterminant>(myBasisSet);
  agpDet->resize(nup, ndown);
  //set Lambda: possible to move to AGPDeterminant
  int offset      = 1;
  xmlNodePtr tcur = cPtr->xmlChildrenNode;
  while (tcur != NULL)
  {
    if (xmlStrEqual(tcur->name, (const xmlChar*)"lambda"))
    {
      int i                                  = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"i")));
      int j                                  = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"j")));
      double c                               = atof((const char*)(xmlGetProp(tcur, (const xmlChar*)"c")));
      agpDet->Lambda(i - offset, j - offset) = c;
      if (i != j)
      {
        agpDet->Lambda(j - offset, i - offset) = c;
      }
    }
    tcur = tcur->next;
  }
  if (uPtr != NULL)
  {
    tcur = uPtr->xmlChildrenNode;
    while (tcur != NULL)
    {
      if (xmlStrEqual(tcur->name, (const xmlChar*)"lambda"))
      {
        int i                                    = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"i")));
        int j                                    = atoi((const char*)(xmlGetProp(tcur, (const xmlChar*)"j")));
        double c                                 = atof((const char*)(xmlGetProp(tcur, (const xmlChar*)"c")));
        agpDet->LambdaUP(j - offset, i - offset) = c;
      }
      tcur = tcur->next;
    }
    app_log() << "  AGPDeterminantBuilder::put Coefficients for the unpaired electrons " << std::endl;
    app_log() << agpDet->LambdaUP << std::endl;
  }
  if (agpDet)
    return std::move(agpDet);

  APP_ABORT("failed to create an AGP determinant!\n");
  return nullptr;
}

} // namespace qmcplusplus
