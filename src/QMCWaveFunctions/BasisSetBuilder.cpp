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
    
    

#include <QMCWaveFunctions/BasisSetBase.h>


namespace qmcplusplus
{

  BasisSetBuilder::BasisSetBuilder()
    : MPIObjectBase(0), legacy(true) 
  {
    reserve_states();
  }


  void BasisSetBuilder::reserve_states(int nsets)
  {
    int sets_needed = nsets - states.size();
    if(sets_needed>0)
      for(int s=0;s<sets_needed;++s)
        states.push_back(new SPOSetInfo());
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input_info)
  { 
    APP_ABORT("BasisSetBase::createSPOSet(cur,input_info) has not been implemented");
    return 0;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(xmlNodePtr cur)
  {
    // read specialized sposet construction requests
    //   and translate them into a set of orbital indices
    SPOSetInputInfo input_info(cur);

    // process general sposet construction requests
    //   and preserve legacy interface 
    SPOSetBase* sposet = 0;
    if(legacy && input_info.legacy_request)
      sposet = createSPOSetFromXML(cur);
    else
      sposet = createSPOSet(cur,input_info);

    // remember created sposets
    if(sposet)
    {
      sposet->put(cur); //initialize C and other internal containers
      sposet->builder_index = sposets.size();
      sposets.push_back(sposet);
    }
    else
      APP_ABORT("BasisSetBuilder::createSPOSet  sposet creation failed");

    return sposet;
  }

}
