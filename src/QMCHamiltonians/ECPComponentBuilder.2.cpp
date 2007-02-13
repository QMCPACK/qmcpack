//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"
#include "QMCHamiltonians/FSAtomPseudoPot.h"
#include "Utilities/IteratorUtility.h"

#define QMCPACK_ABORT abort();

namespace qmcplusplus {

  void ECPComponentBuilder::buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr) {
#if defined(__IBM_CPP__)
    app_error() << "  Compiler error with IBM Compilers. Use other pseudo potential format." << endl;
#else
    if(grid_inp==0)
    {
      app_error() << "  Pseudopotential file does not defined a global grid. vps/grid is disabled." << endl;
      OHMMS::Controller->abort();
    }

    //this is abinit/siesta format
    // There should only be one semilocal tag
    bool is_r_times_V(true);

    if (semiPtr.size()> 1) {
      ERRORMSG("We have more than one semilocal sections in the PP xml file.");
      OHMMS::Controller->abort();
    }

    RealType rmax = 0.0;
    //attributes: initailize by defaults
    string units("hartree");
    string format("r*V");
    string lloc;
    int ndown=1; 
    int nup=0;
    int nrule=4;//default quadrature
    int Llocal = -1;

    OhmmsAttributeSet aAttrib;
    aAttrib.add(units,"units");
    aAttrib.add(format,"format");
    aAttrib.add(ndown,"npots-down");
    aAttrib.add(nup,"npots-up");
    aAttrib.add(Llocal,"l-local");
    aAttrib.add(nrule,"nrule");

    xmlNodePtr cur_semilocal = semiPtr[0];
    aAttrib.put(cur_semilocal);

    RealType Vprefactor=1.0;
    if (units == "rydberg") 
      Vprefactor = 0.5;
    else if (units == "hartree")
      Vprefactor = 1.0;
    else {
      ERRORMSG("Unrecognized units """ << units << """ in PP file.");
      OHMMS::Controller->abort();
    }

    if (format == "r*V") 
      is_r_times_V = true;
    else if (format == "V")
      is_r_times_V = false;
    else {
      ERRORMSG("Unrecognized format """ << format << """ in PP file.");
      OHMMS::Controller->abort();
    }

    // Read which quadrature rule to use
    //const char *nrule = (const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"nrule");
    //Nrule = (nrule==NULL) ? 4 : (int)atol(nrule);
    //cerr << "nrule = " << Nrule << endl;

    // We cannot construct the potentials as we construct them since
    // we may not know which one is local yet.
    typedef FSAtomPseudoPot<RealType> InFuncType;
    vector<InFuncType*> Vls;
    int iLocal=-1;
    // Now read vps sections
    xmlNodePtr cur_vps = cur_semilocal->children;
    while (cur_vps != NULL) {
      string vname ((const char*)cur_vps->name);
      if (vname == "vps") 
      {
        OhmmsAttributeSet aAttrib;
        string lstr("s");
        RealType rc=-1.0;
        aAttrib.add(lstr,"l");
        aAttrib.add(rc,"cutoff");
        aAttrib.put(cur_vps);

        InFuncType* avps=new InFuncType(angMon[lstr],rc,grid_inp);

        avps->put(cur_vps);
        Lmax=std::max(Lmax,avps->AngL);
        rmax=std::max(rmax,rc);

        //multiply by R
        if(!is_r_times_V) avps->convert2RV();

        Vls.push_back(avps);
      }
      cur_vps = cur_vps->next;
    }

    if (Llocal == -1) Llocal = Lmax;    

    Lmax=0;
    // Find which channel is the local one
    for (int i=0; i<Vls.size(); i++)
    {
      int l(Vls[i]->AngL);
      if (l == Llocal)
        iLocal = i;
      else
        Lmax = std::max(l,Lmax);
    }

    pp_loc = Vls[iLocal]->getLocalPot(-Vprefactor/Zeff);

    // Now construct the radial potentials
    NumNonLocal=0;
    for (int i=0; i<Vls.size(); i++) {
      if (i == iLocal) continue;
      RadialPotentialType* newpot = Vls[i]->getNonLocalPot(*Vls[iLocal],Vprefactor);
      pp_nonloc->add(Vls[i]->AngL,newpot);
      NumNonLocal++;
    }

    delete_iter(Vls.begin(),Vls.end());
    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;
#endif
  }

} // namespace qmcPlusPlus
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1551 $   $Date: 2006-12-02 09:32:17 -0600 (Sat, 02 Dec 2006) $
 * $Id: ECPComponentBuilder.cpp 1551 2006-12-02 15:32:17Z jnkim $ 
 ***************************************************************************/
