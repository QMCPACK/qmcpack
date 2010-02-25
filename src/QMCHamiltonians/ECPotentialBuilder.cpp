//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim and Simone Chiesa
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
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/CoulombPBCABTemp.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OneDimNumGridFunctor.h"
#ifdef QMC_CUDA
  #include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
  #include "QMCHamiltonians/LocalECPotential_CUDA.h"
  #include "QMCHamiltonians/NonLocalECPotential_CUDA.h"
#endif

namespace qmcplusplus {
  /** constructor
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\param psi trial wavefunction
   */
  ECPotentialBuilder::ECPotentialBuilder(QMCHamiltonian& h,
      ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi,
      Communicate* c): 
    MPIObjectBase(c),
    hasLocalPot(false),hasNonLocalPot(false),
    targetH(h), IonConfig(ions), targetPtcl(els), targetPsi(psi)
  { }

  bool ECPotentialBuilder::put(xmlNodePtr cur) {

    if(localPot.empty()) {
      int ng(IonConfig.getSpeciesSet().getTotalNum());
      localZeff.resize(ng,1);
      localPot.resize(ng,0);
      nonLocalPot.resize(ng,0);
    }

    string ecpFormat("table");
    string pbc("yes");
    string forces("no");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(ecpFormat,"format");
    pAttrib.add(pbc,"pbc");
    pAttrib.add(forces,"forces");
    pAttrib.put(cur);
    bool doForces = (forces == "yes") || (forces == "true");

    //const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"format");
    //if(t != NULL) {
    //  ecpFormat= (const char*)t;
    //} 

    if(ecpFormat == "xml")  
    {
      useXmlFormat(cur);
    } 
    else 
    {
      useSimpleTableFormat();
    } 

    ///create LocalECPotential
    bool usePBC = 
      !(IonConfig.Lattice.SuperCellEnum == SUPERCELL_OPEN || pbc =="no");
    if(hasLocalPot) {
      if(IonConfig.Lattice.SuperCellEnum == SUPERCELL_OPEN || pbc =="no") 
      {
#ifdef QMC_CUDA
        LocalECPotential_CUDA* apot = 
	  new LocalECPotential_CUDA(IonConfig,targetPtcl);
#else
        LocalECPotential* apot = new LocalECPotential(IonConfig,targetPtcl);
#endif
        for(int i=0; i<localPot.size(); i++) {
          if(localPot[i]) apot->add(i,localPot[i],localZeff[i]);
        }
        targetH.addOperator(apot,"LocalECP");
      }
      else
      {
	if (doForces) 
	  app_log() << "  Will compute forces in CoulombPBCABTemp.\n" << endl;
#ifdef QMC_CUDA
        CoulombPBCAB_CUDA* apot=
	  new CoulombPBCAB_CUDA(IonConfig,targetPtcl, doForces);
#else
	CoulombPBCABTemp* apot =
	  new CoulombPBCABTemp(IonConfig,targetPtcl, doForces);
#endif
        for(int i=0; i<localPot.size(); i++) {
          if(localPot[i]) apot->add(i,localPot[i]);
        }
        targetH.addOperator(apot,"LocalECP");
      }
      //if(IonConfig.Lattice.BoxBConds[0]) {
      //  CoulombPBCABTemp* apot=new CoulombPBCABTemp(IonConfig,targetPtcl);
      //  for(int i=0; i<localPot.size(); i++) {
      //    if(localPot[i]) apot->add(i,localPot[i]);
      //  }
      //  targetH.addOperator(apot,"LocalECP");
      //} else {
      //  LocalECPotential* apot = new LocalECPotential(IonConfig,targetPtcl);
      //  for(int i=0; i<localPot.size(); i++) {
      //    if(localPot[i]) apot->add(i,localPot[i],localZeff[i]);
      //  }
      //  targetH.addOperator(apot,"LocalECP");
      //}
    }

    if(hasNonLocalPot) {
      //resize the sphere
      targetPtcl.resizeSphere(IonConfig.getTotalNum());
      RealType rc2=0.0;
#ifdef QMC_CUDA   
      NonLocalECPotential_CUDA* apot = 
	new NonLocalECPotential_CUDA(IonConfig,targetPtcl,targetPsi,usePBC,doForces);
#else
      NonLocalECPotential* apot = 
	new NonLocalECPotential(IonConfig,targetPtcl,targetPsi, doForces);
#endif

      for(int i=0; i<nonLocalPot.size(); i++) 
      {
        if(nonLocalPot[i]) 
        {
          rc2=std::max(rc2,nonLocalPot[i]->Rmax);
          apot->add(i,nonLocalPot[i]);
        }
      }
      targetPtcl.checkBoundBox(2*rc2);
      targetH.addOperator(apot,"NonLocalECP");

      for(int ic=0; ic<IonConfig.getTotalNum(); ic++) 
      {
        int ig=IonConfig.GroupID[ic];
        if(nonLocalPot[ig]) { 
          if(nonLocalPot[ig]->nknot) targetPtcl.Sphere[ic]->resize(nonLocalPot[ig]->nknot);
        }
      }
    }
    return true;
  }

  void ECPotentialBuilder::useXmlFormat(xmlNodePtr cur) {

    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "pseudo") {
        string href("none");
        string ionName("none");
        string format("xml");
        //RealType rc(2.0);//use 2 Bohr
        OhmmsAttributeSet hAttrib;
        hAttrib.add(href,"href");
        hAttrib.add(ionName,"elementType"); hAttrib.add(ionName,"symbol");
        hAttrib.add(format,"format");
        //hAttrib.add(rc,"cutoff");
        hAttrib.put(cur);

        int speciesIndex=IonConfig.getSpeciesSet().findSpecies(ionName);
        if(speciesIndex < IonConfig.getSpeciesSet().getTotalNum()) {
          app_log() << endl << "  Adding pseudopotential for " << ionName << endl;
          ECPComponentBuilder ecp(ionName,myComm);
          bool success=false;
          if(format == "xml") {
            if(href == "none") {
              success=ecp.put(cur);
            } else {
              success=ecp.parse(href,cur);
            }
          } 
          else if(format == "casino")
          {
            //success=ecp.parseCasino(href,rc);
            success=ecp.parseCasino(href,cur);
          }

          if(success) {
            if(OHMMS::Controller->rank()==0) ecp.printECPTable();
            if(ecp.pp_loc) {
              localPot[speciesIndex]=ecp.pp_loc;
              localZeff[speciesIndex]=ecp.Zeff;
              hasLocalPot=true;
            }
            if(ecp.pp_nonloc) {
              nonLocalPot[speciesIndex]=ecp.pp_nonloc;
              hasNonLocalPot=true;
            }
          }
        } else {
          app_error() << "  Ion species " << ionName << " is not found." << endl;
        }
      } 
      cur=cur->next;
    }
  }

  /** reimplement simple table format used by NonLocalPPotential
   */
  void ECPotentialBuilder::useSimpleTableFormat() {

    SpeciesSet& Species(IonConfig.getSpeciesSet());
    int ng(Species.getTotalNum());
    int icharge(Species.addAttribute("charge"));

    for(int ig=0; ig<ng;ig++) {
      vector<RealType> grid_temp, pp_temp;
      string species(Species.speciesName[ig]);
      string fname = species+".psf";

      ifstream fin(fname.c_str(),ios_base::in);
      if(!fin){
	ERRORMSG("Could not open file " << fname)
        exit(-1);
      }      

      // Read Number of potentials (local and non) for this atom
      int npotentials;
      fin >> npotentials;
      RealType r, f1;

      int lmax=-1;
      int numnonloc=0;
      RealType rmax(0.0);

      app_log() << "  ECPotential for " << species << endl;
      NonLocalECPComponent* mynnloc=0;

      typedef OneDimCubicSpline<RealType> CubicSplineFuncType;
      for (int ij=0; ij<npotentials; ij++){
	int angmom,npoints;
	fin >> angmom >> npoints;

        OneDimNumGridFunctor<RealType> inFunc;
        inFunc.put(npoints,fin);

	if(angmom < 0) {//local potential, input is rescale by -r/z

          RealType zinv=-1.0/Species(icharge,ig);
          int ng=npoints-1;
          RealType rf=5.0;
          ng=static_cast<int>(rf*100)+1;//use 1e-2 resolution
          GridType * agrid= new LinearGrid<RealType>;
          agrid->set(0,rf,ng);
          vector<RealType> pp_temp(ng);
          pp_temp[0]=0.0;
          for (int j=1; j<ng; j++){
            RealType r((*agrid)[j]);
            pp_temp[j]=r*zinv*inFunc.splint(r);
          }
          pp_temp[ng-1]=1.0;
          RadialPotentialType *app = new RadialPotentialType(agrid,pp_temp);
	  app->spline();
          localPot[ig]=app;
          app_log() << "    LocalECP l=" << angmom << endl;
          app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << endl;
          hasLocalPot=true; //will create LocalECPotential
        } else {
          hasNonLocalPot=true; //will create NonLocalECPotential
          if(mynnloc == 0) mynnloc = new NonLocalECPComponent;

          RealType rf=inFunc.rmax();

          GridType * agrid= new LinearGrid<RealType>;
          int ng=static_cast<int>(rf*100)+1;
          agrid->set(0.0,rf,ng);
          app_log() << "    NonLocalECP l=" << angmom << " rmax = " << rf << endl;
          app_log() << "      Linear grid=[0," << rf << "] npts=" << ng << endl;

          vector<RealType> pp_temp(ng);
          //get the value
          pp_temp[0]=inFunc(0);
          for (int j=1; j<ng; j++){
            pp_temp[j]=inFunc.splint((*agrid)[j]);
          }

          RadialPotentialType *app = new RadialPotentialType(agrid,pp_temp);
	  app->spline();

	  mynnloc->add(angmom,app);
	  lmax=std::max(lmax,angmom);
	  rmax=std::max(rmax,rf);
          numnonloc++;
	}

        if(mynnloc) {
          mynnloc->lmax=lmax; 
          mynnloc->Rmax=rmax;
          app_log() << "    Maximum cutoff of NonLocalECP " << rmax << endl;
        }
      } 
      fin.close();

      if(mynnloc) {
        nonLocalPot[ig]=mynnloc;
        int numsgridpts=0;

        string fname = species+".sgr";
        ifstream fin(fname.c_str(),ios_base::in);
        if(!fin){
          app_error() << "Could not open file " << fname << endl;
          exit(-1);
        }
        PosType xyz;
        RealType weight;
        while(fin >> xyz >> weight){
          mynnloc->addknot(xyz,weight);
          numsgridpts++;
        }
        //cout << "Spherical grid : " << numsgridpts << " points" <<endl;
        mynnloc->resize_warrays(numsgridpts,numnonloc,lmax);
      }
    }//species
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
