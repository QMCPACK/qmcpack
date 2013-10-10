//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/BasisSetBase.h>

namespace qmcplusplus
{

  BasisSetBuilder::energies_t& BasisSetBuilder::get_energies()
  {
    APP_ABORT("BasisSetBase::get_energies has not been implemented");
    return energy_list;
  }


  BasisSetBuilder::degeneracies_t& BasisSetBuilder::get_degeneracies()
  {
    APP_ABORT("BasisSetBase::get_degeneracies has not been implemented");
    return degeneracy_list;
  }


  SPOSetBase* BasisSetBuilder::createSPOSetFromStates(states_t& states)
  { 
    APP_ABORT("BasisSetBase::createSPOSet(states) has not been implemented");
    return 0;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(xmlNodePtr cur)
  {
    SPOSetBase* sposet = createSPOSetFromXML(cur);
    if(sposet)
      sposets.push_back(sposet);
    initialized = true;
    return sposet;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(states_t& states)
  {
    if(!initialized)
      APP_ABORT("BasisSetBuilder::createSPOSetAndSave(states)  an SPOSet must be created from xml to initialize BasisSetBuilder before creating one for a selection of states");
    SPOSetBase* sposet = createSPOSetFromStates(states);
    if(sposet)
      sposets.push_back(sposet);
    return sposet;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(int range_max) 
  {
    return createSPOSet(0,range_max);
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(int range_min,int range_max) 
  { 
    states_t states;
    for(int s=range_min;s<range_max;++s)
      states.push_back(s);
    if(states.size()>0)
      return createSPOSet(states);
    else
      return 0;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(RealType emax,RealType tol) 
  { 
    return createSPOSet(-1e99,emax,tol);
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(RealType emin,RealType emax,RealType tol) 
  { 
    energies_t& energies = get_energies();
    states_t states;
    for(int s=0;s<energies.size();++s)
    {
      RealType& e = energies[s];
      if(e<emax+tol && e>=emin-tol)
        states.push_back(s);
    }
    if(states.size()>0)
      return createSPOSet(states);
    else
      return 0;
  }


  SPOSetBase* BasisSetBuilder::createSPOSet(energies_t& energies_in,RealType tol)
  { 
    sort(energies_in.begin(),energies_in.end());
    energies_t& energies = get_energies();
    degeneracies_t& degeneracies = get_degeneracies();
    int s=0;
    states_t states;
    for(int n=0;n<energies_in.size();++n)
    {
      RealType& e = energies_in[n];
      bool found = false;
      while(s<energies.size())
      {
        if(abs(e-energies[s])<tol)
        {
          found = true;
          break;
        }
        s++;
      }
      if(!found)
        APP_ABORT("BasisSetBuilder::createSPOSet(energies)  energy eigenvalue not found");
      int degeneracy = degeneracies[s];
      for(int i=0;i<degeneracy;++i)
      {
        states.push_back(s);
        s++;
      }
    }
    if(states.size()>0)
      return createSPOSet(states);
    else
      return 0;
  }


}
