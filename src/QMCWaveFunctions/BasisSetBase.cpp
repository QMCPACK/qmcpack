//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/BasisSetBase.h>
#include <QMCWaveFunctions/SPOSetInputInfo.h>


namespace qmcplusplus
{

  BasisSetBuilder::BasisSetBuilder()
  : MPIObjectBase(0), myBasisSet(0) 
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


  SPOSetBase* BasisSetBuilder::createSPOSetFromIndices(indices_t& indices)
  { 
    APP_ABORT("BasisSetBase::createSPOSet(indices) has not been implemented");
    return 0;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(xmlNodePtr cur)
  {
    // setup orbital data for maximal basis set
    Initialize(cur);

    // read specialized sposet construction requests
    //   and translate them into a set of orbital indices
    SPOSetInputInfo input_info(cur);
    indices_t& indices = input_info.get_indices(states);

    // process general sposet construction requests (from indices)
    //   and preserve legacy interface (from xml, may be removed later)
    SPOSetBase* sposet = 0;
    if(indices.size()>0)
      sposet = createSPOSetFromIndices(indices);
    else
      sposet = createSPOSetFromXML(cur);

    // remember created sposets
    if(sposet)
    {
      sposet->builder_index = sposets.size();
      sposets.push_back(sposet);
    }

    return sposet;
  }

}
