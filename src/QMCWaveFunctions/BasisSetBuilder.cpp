//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

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
