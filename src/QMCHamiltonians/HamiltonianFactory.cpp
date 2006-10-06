/////////////////////////////////////////////////////////////////
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
/**@file HamiltonianFactory.cpp
 *@brief Definition of a HamiltonianFactory 
 */
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/HarmonicPotential.h"
#include "QMCHamiltonians/LocalCorePolPotential.h"
#if !defined(QMC_COMPLEX)
#include "QMCHamiltonians/NonLocalPPotential.h"
#include "QMCHamiltonians/ECPotentialBuilder.h"
#endif
#if defined(HAVE_LIBFFTW)
#include "QMCHamiltonians/ModInsKineticEnergy.h"
#include "QMCHamiltonians/MomentumDistribution.h"
#include "QMCHamiltonians/DispersionRelation.h"
#endif
//CoulombPBC is replaced by CoulombPBCAATemp and CoulombPBCABTemp
//#include "QMCHamiltonians/CoulombPBC.h"
#include "QMCHamiltonians/CoulombPBCAATemp.h"
#include "QMCHamiltonians/CoulombPBCABTemp.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {
  HamiltonianFactory::HamiltonianFactory(ParticleSet* qp, 
    PtclPoolType& pset, OrbitalPoolType& oset): 
    targetPtcl(qp), targetH(0), 
  ptclPool(pset),psiPool(oset), myNode(NULL), psiName("psi0") 
  {
    //PBCType is zero or 1 but should be generalized 
    PBCType=targetPtcl->Lattice.BoxBConds[0];
  }

  /** main hamiltonian build function
   * @param cur element node <hamiltonian/>
   * @param buildtree if true, build xml tree for a reuse
   *
   * A valid hamiltonian node contains
   * \xmlonly
   *  <hamiltonian target="e">
   *    <pairpot type="coulomb" name="ElecElec" source="e"/>
   *    <pairpot type="coulomb" name="IonElec" source="i"/>
   *    <pairpot type="coulomb" name="IonIon" source="i" target="i"/>
   *  </hamiltonian>
   * \endxmlonly
   */
  bool HamiltonianFactory::build(xmlNodePtr cur, bool buildtree) {

    if(cur == NULL) return false;

    string htype("generic"), source("i"), defaultKE("yes");
    OhmmsAttributeSet hAttrib;
    hAttrib.add(htype,"type"); 
    hAttrib.add(source,"source");
    hAttrib.add(defaultKE,"default");
    hAttrib.put(cur);

    renameProperty(source);

    bool attach2Node=false;
    if(buildtree) {
      if(myNode == NULL) {
//#if (LIBXMLD_VERSION < 20616)
//        app_warning() << "   Workaround of libxml2 bug prior to 2.6.x versions" << endl;
//        myNode = xmlCopyNode(cur,2);
//#else
//        app_warning() << "   using libxml2 2.6.x versions" << endl;
//        myNode = xmlCopyNode(cur,1);
//#endif
        myNode = xmlCopyNode(cur,1);
      } else {
        attach2Node=true;
      }
    }

    if(targetH==0) {
      targetH  = new QMCHamiltonian;
      targetH->setName(myName);
      
      if(defaultKE == "yes")
        targetH->addOperator(new BareKineticEnergy,"Kinetic");
    }

    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      const xmlChar* t = xmlGetProp(cur,(const xmlChar*)"type");
      if(t != NULL) { // accept only if it has type
        string pot_type((const char*)t);
        string nuclei("i");
        const xmlChar* sptr = xmlGetProp(cur,(const xmlChar*)"source");
        if(sptr != NULL) nuclei=(const char*)sptr;
        if(cname == "pairpot") {
          if(pot_type == "coulomb") {
            bool sameTarget=true;
            string aNewTarget(targetPtcl->getName());
            const xmlChar* aptr = xmlGetProp(cur,(const xmlChar*)"target");
            if(aptr != NULL) {
              aNewTarget=(const char*)aptr;
              sameTarget= (aNewTarget == targetPtcl->getName());
            } 
            if(sameTarget) 
	      addCoulombPotential(cur);
            else {
              app_log() << "  Creating Coulomb potential " << nuclei << "-" << nuclei << endl;
              addConstCoulombPotential(cur,nuclei);
            }
          } else if(pot_type == "pseudo") {
            addPseudoPotential(cur);
          } else if(pot_type == "cpp") {
            addCorePolPotential(cur);
          }
        } 
        else if(cname == "harmonic") {
          PtclPoolType::iterator pit(ptclPool.find(nuclei));
          if(pit != ptclPool.end()) {
            ParticleSet* ion=(*pit).second;
            targetH->addOperator(new HarmonicPotential(*ion, *targetPtcl),"Harmonic");
            app_log() << "  Adding HarmonicPotential " << endl;
          }
        } else if(cname == "constant") { 
          if(pot_type == "coulomb") { //ugly!!!
            addConstCoulombPotential(cur,nuclei);
          }
        } else if(cname == "modInsKE") {
          addModInsKE(cur);
        }
      }
      if(attach2Node) xmlAddChild(myNode,xmlCopyNode(cur,1));
      cur = cur->next;
    }

    //This should be disabled
    if(targetH->size() == 1) {//no external potential is provided, use type
      WARNMSG("Using pre-determined hamiltonian for molecular systems.")
      PtclPoolType::iterator pit(ptclPool.find(source));
      if(pit == ptclPool.end()) {
        ERRORMSG("No ionic system " << source << " exists.")
        return false;
      }
      ParticleSet* ion=(*pit).second;
      if(PBCType) targetH->addOperator(new CoulombPBCAATemp(*targetPtcl),"ElecElec");
      else targetH->addOperator(new CoulombPotentialAA(*targetPtcl),"ElecElec");
      if(htype == "molecule" || htype=="coulomb"){
        if(PBCType) targetH->addOperator(new CoulombPBCABTemp(*ion,*targetPtcl),"Coulomb");
        else targetH->addOperator(new CoulombPotentialAB(*ion,*targetPtcl),"Coulomb");
      } else {
        ERRORMSG(htype << " is diabled")
      }
      if(ion->getTotalNum()>1) {
        if(PBCType) targetH->addOperator(new CoulombPBCAATemp(*ion),"IonIon");
        else targetH->addOperator(new IonIonPotential(*ion),"IonIon");
      }
    }

    return true;
  }

  void 
  HamiltonianFactory::addCoulombPotential(xmlNodePtr cur) {

    string a("e"),title("ElecElec");
    OhmmsAttributeSet hAttrib;
    hAttrib.add(title,"id"); hAttrib.add(title,"name"); 
    hAttrib.add(a,"source"); 
    hAttrib.put(cur);

    renameProperty(a);

    PtclPoolType::iterator pit(ptclPool.find(a));
    if(pit == ptclPool.end()) {
      ERRORMSG("Missing source ParticleSet" << a)
      return;
    }

    ParticleSet* source = (*pit).second;

    //CHECK PBC and create CoulombPBC for el-el
    if(source == targetPtcl) {
      if(source->getTotalNum()>1)  {
        if(PBCType) {
          //targetH->addOperator(new CoulombPBCAA(*targetPtcl),title);
          targetH->addOperator(new CoulombPBCAATemp(*targetPtcl),title);
        } else {
          targetH->addOperator(new CoulombPotentialAA(*targetPtcl),title);
        }
      }
    } else {
      if(PBCType) {
        //targetH->addOperator(new CoulombPBCAB(*source,*targetPtcl),title);
        targetH->addOperator(new CoulombPBCABTemp(*source,*targetPtcl),title);
      } else {
        targetH->addOperator(new CoulombPotentialAB(*source,*targetPtcl),title);
      }
    }
  }

  void 
  HamiltonianFactory::addPseudoPotential(xmlNodePtr cur) {

#if !defined(QMC_COMPLEX)
    string src("i"),title("PseudoPot"),wfname("invalid"),format("old");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.add(wfname,"wavefunction");
    pAttrib.add(format,"format"); //temperary tag to switch between format
    pAttrib.put(cur);

    renameProperty(src);
    renameProperty(wfname);

    PtclPoolType::iterator pit(ptclPool.find(src));
    if(pit == ptclPool.end()) {
      ERRORMSG("Missing source ParticleSet" << src)
      return;
    }

    ParticleSet* ion=(*pit).second;

    OrbitalPoolType::iterator oit(psiPool.find(wfname));
    TrialWaveFunction* psi=0;
    if(oit == psiPool.end()) {
      if(psiPool.empty()) return;
      app_error() << "  Cannot find " << wfname << " in the Wavefunction pool. Using the first wavefunction."<< endl;
      psi=(*(psiPool.begin())).second->targetPsi;
    } else {
      psi=(*oit).second->targetPsi;
    }

    //remember the TrialWaveFunction used by this pseudopotential
    psiName=wfname;

    if(format == "old") {
      app_log() << "  Using OLD NonLocalPseudopotential "<< endl;
      targetH->addOperator(new NonLocalPPotential(*ion,*targetPtcl,*psi), title);
    }
    else  {
      app_log() << "  Using ECPotential builder "<< endl;
      ECPotentialBuilder ecp(*targetH,*ion,*targetPtcl,*psi);
      ecp.put(cur);
    }
#endif
  }

  void 
  HamiltonianFactory::addCorePolPotential(xmlNodePtr cur) {

    string src("i"),title("CorePol");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.put(cur);

    PtclPoolType::iterator pit(ptclPool.find(src));
    if(pit == ptclPool.end()) {
      ERRORMSG("Missing source ParticleSet" << src)
      return;
    }
    ParticleSet* ion=(*pit).second;

    QMCHamiltonianBase* cpp=(new LocalCorePolPotential(*ion,*targetPtcl));
    cpp->put(cur); 
    targetH->addOperator(cpp, title);
  }

  void 
  HamiltonianFactory::addConstCoulombPotential(xmlNodePtr cur, string& nuclei){
    renameProperty(nuclei);
    PtclPoolType::iterator pit(ptclPool.find(nuclei));
    if(pit != ptclPool.end()) {
      ParticleSet* ion=(*pit).second;
      if(ion->getTotalNum()>1) 
        if(PBCType){
          //targetH->addOperator(new CoulombPBCAA(*ion),"IonIon");
          targetH->addOperator(new CoulombPBCAATemp(*ion),"IonIon");
        } else {
          targetH->addOperator(new IonIonPotential(*ion),"IonIon");
        }
     }
  }

  void
  HamiltonianFactory::addModInsKE(xmlNodePtr cur) {
#if defined(HAVE_LIBFFTW)
    string Dimensions, DispRelType, MomDistType, PtclSelType;
    RealType Cutoff, GapSize(0.0), FermiMomentum(0.0);
    
    OhmmsAttributeSet pAttrib;
    pAttrib.add(Dimensions, "dims");
    pAttrib.add(DispRelType, "dispersion");
    pAttrib.add(PtclSelType, "selectParticle");
    pAttrib.add(Cutoff, "cutoff");
    pAttrib.add(GapSize, "gapSize");
    pAttrib.add(FermiMomentum, "kf");
    pAttrib.put(cur);
    
    TrialWaveFunction* psi;
    psi = (*(psiPool.begin())).second->targetPsi;
    
    Vector<PosType> LocLattice;
    Vector<IndexType> DimSizes;
    if (Dimensions == "3") {
      gen3DLattice(Cutoff, *targetPtcl, LocLattice, DimSizes);
    } else if (Dimensions == "1" || Dimensions == "1averaged") {
      gen1DLattice(Cutoff, *targetPtcl, LocLattice, DimSizes);
    } else {
      ERRORMSG("Dimensions value not recognized!")
    }
    
    Vector<RealType> Dispersion;
    if (DispRelType == "freeParticle") {
      genFreeParticleDispersion(LocLattice, Dispersion);
    } else if (DispRelType == "simpleModel") {
      genSimpleModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
    } else if (DispRelType == "pennModel") {
      genPennModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
    } else {
      ERRORMSG("Dispersion relation not recognized");
    }
    
    QMCHamiltonianBase* modInsKE;
    if (Dimensions == "3") { 
      if (PtclSelType == "random") {
        modInsKE = new ModInsKineticEnergy<ThreeDimMomDist<RandomChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "randomPerWalker") {
        modInsKE = new ModInsKineticEnergy<ThreeDimMomDist<RandomChoicePerWalkerPolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "constant") {
        modInsKE = new ModInsKineticEnergy<ThreeDimMomDist<StaticChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      }
    } else if (Dimensions == "1") {
      if (PtclSelType == "random") {
        modInsKE = new ModInsKineticEnergy<OneDimMomDist<RandomChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "randomPerWalker") {
        modInsKE = new ModInsKineticEnergy<OneDimMomDist<RandomChoicePerWalkerPolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "constant") {
        modInsKE = new ModInsKineticEnergy<OneDimMomDist<StaticChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      }
    } else if (Dimensions == "1averaged") {
      if (PtclSelType == "random") {
        modInsKE = new ModInsKineticEnergy<AveragedOneDimMomDist<RandomChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "randomPerWalker") {
        modInsKE = new ModInsKineticEnergy<AveragedOneDimMomDist<RandomChoicePerWalkerPolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      } else if (PtclSelType == "constant") {
        modInsKE = new ModInsKineticEnergy<AveragedOneDimMomDist<StaticChoicePolicy> >
          (*targetPtcl, *psi, Dispersion, DimSizes);
      }
    }
    modInsKE->put(cur);
    targetH->addOperator(modInsKE, "ModelInsulatorKE");
#else
    app_error() << "  ModelInsulatorKE cannot be used without FFTW " << endl;
#endif
  }
  
  void HamiltonianFactory::renameProperty(const string& a, const string& b){
    RenamedProperty[a]=b;
  }

  void HamiltonianFactory::setCloneSize(int np) {
    myClones.resize(np,0);
  }

  //TrialWaveFunction*
  //HamiltonianFactory::cloneWaveFunction(ParticleSet* qp, int ip) {
  //  HamiltonianFactory* aCopy= new HamiltonianFactory(qp,ptclPool);
  //  aCopy->put(myNode,false);
  //  myClones[ip]=aCopy;
  //  return aCopy->targetPsi;
  //}

  void HamiltonianFactory::renameProperty(string& aname) {
    map<string,string>::iterator it(RenamedProperty.find(aname));
    if(it != RenamedProperty.end()) {
      aname=(*it).second;
    }
  }
  HamiltonianFactory::~HamiltonianFactory() {
    //clean up everything
  }

  HamiltonianFactory*
  HamiltonianFactory::clone(ParticleSet* qp, TrialWaveFunction* psi, 
      int ip, const string& aname) {
    HamiltonianFactory* aCopy=new HamiltonianFactory(qp, ptclPool, psiPool);
    aCopy->setName(aname);

    aCopy->renameProperty("e",qp->getName());
    aCopy->renameProperty(psiName,psi->getName());

    aCopy->build(myNode,false);
    myClones[ip]=aCopy;

    aCopy->get(app_log());
    return aCopy;
  }

  bool HamiltonianFactory::get(std::ostream& os) const {
    targetH->get(os);
    return true;
  }

  bool HamiltonianFactory::put(std::istream& ) {
    return true;
  }

  bool HamiltonianFactory::put(xmlNodePtr cur) {
    return build(cur,true);
  }

  void HamiltonianFactory::reset() { }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
