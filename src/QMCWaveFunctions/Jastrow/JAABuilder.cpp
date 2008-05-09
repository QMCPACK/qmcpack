//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/JAABuilder.h"
#include "QMCWaveFunctions/Jastrow/ModPadeFunctor.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

   JAABuilder::JAABuilder(ParticleSet& p, TrialWaveFunction& psi):OrbitalBuilderBase(p,psi),
    IgnoreSpin(true){
    }
  /** Create a two-body Jatrow function with a template 
   *@param cur the current xmlNode 
   *@param dummy null pointer used to identify FN
   *
   *The template class JeeType is a functor which handles the
   *evaluation of the function value, gradients and laplacians using
   *distance tables. This is a specialized builder function for
   *spin-dependent Jastrow function,e.g., for electrons, two functions
   *are created for uu(dd) and ud(du).
   */
  template<class FN> bool JAABuilder::createJAA(xmlNodePtr cur, const string& jname) 
  {

    string corr_tag("correlation");
    int ng = targetPtcl.groups();

    int ia=0, ib=0, iab=0;
    int cur_var = targetPsi.VarList.size();
    xmlNodePtr gridPtr=NULL;
    cur = cur->children;
    const SpeciesSet& species(targetPtcl.getSpeciesSet());
    typedef TwoBodyJastrowOrbital<FN> JeeType;
    JeeType *J2 = new JeeType(targetPtcl);

    int pairs=0;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == corr_tag) 
      {
        string spA("u");
        string spB("u");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(spA, "speciesA"); rAttrib.add(spA, "species1");
        rAttrib.add(spB, "speciesB"); rAttrib.add(spA, "species2");
        rAttrib.put(cur);
        if(spA==targetPsi.getName()) //could have used the particle name
        {
          spA=species.speciesName[0];
          spB=species.speciesName[0];
        }
        int ia = species.findSpecies(spA);
        int ib = species.findSpecies(spB);
        if(ia==species.size() || ia == species.size())
        {
          APP_ABORT("JAABuilder::createJAA is trying to use invalid species");
        }
        string pairID=spA+spB;
        FN *j= new FN;
        j->put(cur);
        j->addOptimizables(targetPsi.VarList);
        J2->addFunc(pairID,ia,ib,j);
        ++pairs;
      }
      cur = cur->next;
    } // while cur

    if(pairs)
    {
      //set this jastrow function to be not optimizable
      if(targetPsi.VarList.size() == cur_var) J2->setOptimizable(false);

      string j2name="J2_"+jname;
      targetPsi.addOrbital(J2,j2name);
      return true;
    }
    else
    {//clean up and delete the twobody orbitals
      delete J2;
      return false;
    }
  }

  bool JAABuilder::put(xmlNodePtr cur) {

    string spinOpt("no");
    string typeOpt("Two-Body");
    string jastfunction("pade");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(spinOpt,"spin");
    aAttrib.add(typeOpt,"type");
    aAttrib.add(jastfunction,"function");
    aAttrib.put(cur);

    IgnoreSpin=(spinOpt=="no");
    bool success=false;
    //if(jastfunction == "pade") {
    //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
    //  PadeJastrow<RealType> *dummy = 0;
    //  success = createJAA(cur,dummy);
    //} else 
    if(jastfunction == "short") {
      app_log() << "  Modified Jastrow function Two-Body Jastrow Function = " << jastfunction << endl;
      IgnoreSpin=true;
      //ModPadeFunctor<RealType> *dummy = 0;
      success = createJAA<ModPadeFunctor<RealType> >(cur,jastfunction);
    }
    //} else if(jastfunction == "rpa") {
    //  app_log() << "  Two-Body Jastrow Function = " << jastfunction << endl;
    //  RPAJastrow<RealType> *dummy = 0;
    //  success = createJAA(cur,dummy);
    //}
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
