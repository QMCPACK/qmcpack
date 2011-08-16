/////////////////////////////////////////////////////////////////
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
/**@file HamiltonianFactory.cpp
 *@brief Definition of a HamiltonianFactory 
 */
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/CoulombPotential.h"
#include "QMCHamiltonians/IonIonPotential.h"
#include "QMCHamiltonians/NumericalRadialPotential.h"
#include "QMCHamiltonians/MomentumEstimator.h"
#include "QMCHamiltonians/CoulombPBCAATemp.h"
#include "QMCHamiltonians/CoulombPBCABTemp.h"
#include "QMCHamiltonians/Pressure.h"
#include "QMCHamiltonians/RPAPressure.h"
#include "QMCHamiltonians/PsiValue.h"
#include "QMCHamiltonians/DMCPsiValue.h"
#include "QMCHamiltonians/PsiOverlap.h"
#include "QMCHamiltonians/ForwardWalking.h"
#include "QMCHamiltonians/trialDMCcorrection.h"
#include "QMCHamiltonians/PairCorrEstimator.h"
#include "QMCHamiltonians/LocalMomentEstimator.h"
#include "QMCHamiltonians/DensityEstimator.h"
#include "QMCHamiltonians/SkEstimator.h"
#include "QMCHamiltonians/MomentumEstimator.h"
#if OHMMS_DIM == 3
#include "QMCHamiltonians/LocalCorePolPotential.h"
#include "QMCHamiltonians/ECPotentialBuilder.h"
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/ForceCeperley.h"
#include "QMCHamiltonians/PulayForce.h"
#include "QMCHamiltonians/ZeroVarianceForce.h"
#include "QMCHamiltonians/ChiesaCorrection.h"
#if defined(HAVE_LIBFFTW)
  #include "QMCHamiltonians/MPC.h"
  #include "QMCHamiltonians/VHXC.h"
#endif
#if defined(HAVE_LIBFFTW_LS)
  #include "QMCHamiltonians/ModInsKineticEnergy.h"
  #include "QMCHamiltonians/MomentumDistribution.h"
  #include "QMCHamiltonians/DispersionRelation.h"
#endif
#endif
// #include "QMCHamiltonians/ZeroVarObs.h"
#if QMC_BUILD_LEVEL>2
#include "QMCHamiltonians/HardSphere.h"
#include "QMCHamiltonians/GaussianPot.h"
#include "QMCHamiltonians/OscillatoryPot.h"
#include "QMCHamiltonians/SkPot.h"
#include "QMCHamiltonians/ModPosTelPot.h"
#include "QMCHamiltonians/HFDHE2Potential_tail.h"
#include "QMCHamiltonians/HePressure.h"
#include "QMCHamiltonians/JelliumPotential.h"
#include "QMCHamiltonians/HFDHE2Potential.h"
#include "QMCHamiltonians/HeEPotential.h"
#include "QMCHamiltonians/HeEPotential_tail.h"
#include "QMCHamiltonians/LennardJones_smoothed.h"
#include "QMCHamiltonians/HFDHE2_Moroni1995.h"
//#include "QMCHamiltonians/HFDBHE_smoothed.h"
#include "QMCHamiltonians/HeSAPT_smoothed.h"
#endif

#include "OhmmsData/AttributeSet.h"

#ifdef QMC_CUDA
  #include "QMCHamiltonians/CoulombPBCAA_CUDA.h"
  #include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
  #include "QMCHamiltonians/CoulombPotential_CUDA.h"
  #include "QMCHamiltonians/MPC_CUDA.h"
  #include "QMCHamiltonians/SkEstimator_CUDA.h"
#endif

namespace qmcplusplus {
  HamiltonianFactory::HamiltonianFactory(ParticleSet* qp, 
      PtclPoolType& pset, OrbitalPoolType& oset, Communicate* c)
    : MPIObjectBase(c), targetPtcl(qp), targetH(0)
      , ptclPool(pset),psiPool(oset), myNode(NULL), psiName("psi0") 
  {
    //PBCType is zero or 1 but should be generalized 
    PBCType=targetPtcl->Lattice.SuperCellEnum;

    ClassName="HamiltonianFactory";
    myName="psi0";
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
      
      if(defaultKE == "yes"){
        targetH->addOperator(new BareKineticEnergy ,"Kinetic");
      } else if(defaultKE == "multi"){
        //Multicomponent wavefunction. Designed for 2 species.
	app_log()<<" Multicomponent system. You must add the Kinetic energy term first!"<<endl;
      } else  {
        double mass(1.0);
        string tgt("mass");
        int indx1 = targetPtcl->mySpecies.findSpecies(defaultKE);
        int indx2 = targetPtcl->mySpecies.addAttribute(tgt);
        mass = targetPtcl->mySpecies(indx2,indx1);
        app_log()<<"  Kinetic energy operator:: Mass "<<mass<<endl;
        targetH->addOperator(new BareKineticEnergy(mass),"Kinetic");
      }
    }

    xmlNodePtr cur_saved(cur);
    cur = cur->children;
    while(cur != NULL) 
    {
      string cname((const char*)cur->name);
      string potType("0");
      string potName("any");
      string potUnit("hartree");
      string estType("coulomb");
      string sourceInp(targetPtcl->getName());
      string targetInp(targetPtcl->getName());
      OhmmsAttributeSet attrib;
      attrib.add(sourceInp,"source");
      attrib.add(sourceInp,"sources");
      attrib.add(targetInp,"target");
      attrib.add(potType,"type");
      attrib.add(potName,"name");
      attrib.add(potUnit,"units");
      attrib.add(estType,"potential");
      attrib.put(cur);
      renameProperty(sourceInp);
      renameProperty(targetInp);
      if(cname == "pairpot") 
      {
        if(potType == "coulomb") 
        {
          if(targetInp == targetPtcl->getName())
            addCoulombPotential(cur);
          else 
            addConstCoulombPotential(cur,sourceInp);
        }
#if QMC_BUILD_LEVEL>2
        else if (potType == "hardsphere")
        {
          HardSphere* hs = new HardSphere(*targetPtcl);
          hs->put(cur);
          targetH->addOperator(hs,"HardSphere",true);
        }        
        else if (potType == "gaussian")
        {
          GaussianPot* hs = new GaussianPot(*targetPtcl);
          hs->put(cur);
          targetH->addOperator(hs,"GaussianPot",true);
        }
        else if (potType == "modpostel")
        {
          ModPoschlTeller* hs = new ModPoschlTeller(*targetPtcl);
          hs->put(cur);
          targetH->addOperator(hs,"ModPoschlTeller",true);
        }
        else if (potType == "oscillatory")
        {
          OscillatoryPotential* hs = new OscillatoryPotential(*targetPtcl);
          hs->put(cur);
          targetH->addOperator(hs,"OscillatoryPotential",true);
        }
        else if (potType == "skpot")
        {
          SkPot* hs = new SkPot(*targetPtcl);
          hs->put(cur);
          targetH->addOperator(hs,"SkPot",true);
        }
#endif
#if OHMMS_DIM==3
	/*
	else if (potType == "HFDBHE_smoothed") {
	  HFDBHE_smoothed_phy* HFD = new HFDBHE_smoothed_phy(*targetPtcl);
	  targetH->addOperator(HFD,"HFD-B(He)",true);
	  HFD->addCorrection(*targetH);
	}
	*/
	else if (potType == "MPC" || potType == "mpc")
	  addMPCPotential(cur);
	else if (potType == "VHXC" || potType == "vhxc")
	  addVHXCPotential(cur);
        else if(potType == "pseudo") 
        {
          addPseudoPotential(cur);
        } 
        else if(potType == "cpp") 
        {
          addCorePolPotential(cur);
        }
#if QMC_BUILD_LEVEL>2
        else if (potType == "LJP_smoothed") {
          LennardJones_smoothed_phy* LJP = new LennardJones_smoothed_phy(*targetPtcl);
          targetH->addOperator(LJP,"LJP",true);
          LJP->addCorrection(*targetH);
        }
        else if (potType == "HeSAPT_smoothed") {
          HeSAPT_smoothed_phy* SAPT = new HeSAPT_smoothed_phy(*targetPtcl);
          targetH->addOperator(SAPT,"HeSAPT",true);
          SAPT->addCorrection(*targetH);
        }
        else if (potType == "HFDHE2_Moroni1995") {
          HFDHE2_Moroni1995_phy* HFD = new HFDHE2_Moroni1995_phy(*targetPtcl);
          targetH->addOperator(HFD,"HFD-HE2",true);
          HFD->addCorrection(*targetH);
        }
        else if(potType == "eHe")
        {
          string SourceName = "e";
          OhmmsAttributeSet hAttrib;
          hAttrib.add(SourceName, "source");
          hAttrib.put(cur);

          PtclPoolType::iterator pit(ptclPool.find(SourceName));
          if(pit == ptclPool.end()) 
          {
            APP_ABORT("Unknown source \"" + SourceName + "\" for e-He Potential.");
          }
          ParticleSet* source = (*pit).second;

          HeePotential* eHetype = new HeePotential(*targetPtcl, *source);
          targetH->addOperator(eHetype,potName,true);
          // 	  targetH->addOperator(eHetype->makeDependants(*targetPtcl),potName,false);

        }
        else if(potType == "jellium")
        {
          string SourceName = "e";
          OhmmsAttributeSet hAttrib;
          hAttrib.add(SourceName, "source");
          hAttrib.put(cur);

          PtclPoolType::iterator pit(ptclPool.find(SourceName));
          if(pit == ptclPool.end()) 
          {
            APP_ABORT("Unknown source \"" + SourceName + "\" for e-He Potential.");
          }
          ParticleSet* source = (*pit).second;

          JelliumPotential* JP = new JelliumPotential(*source, *targetPtcl);
          targetH->addOperator(JP,potName,true);
          //    targetH->addOperator(eHetype->makeDependants(*targetPtcl),potName,false);

        }
        else if(potType == "HFDHE2") 
        {
          HFDHE2Potential* HFD = new HFDHE2Potential(*targetPtcl);
          targetH->addOperator(HFD,"HFDHE2",true);
          //HFD->addCorrection(*targetPtcl,*targetH);
          targetH->addOperator(HFD->makeDependants(*targetPtcl),HFD->depName,false);
          app_log() << "  Adding HFDHE2Potential(Au) " << endl;
        }
#endif
#endif
        else if(potType.find("num") < potType.size())
        {
          if(sourceInp == targetInp)//only accept the pair-potential for now
          {
            NumericalRadialPotential* apot=new NumericalRadialPotential(*targetPtcl);
            apot->put(cur);
            targetH->addOperator(apot,potName);
          }
        }
      } 
      else if(cname == "constant") 
      { //ugly!!!
        if(potType == "coulomb")  addConstCoulombPotential(cur,sourceInp);
      } 
      else if(cname == "modInsKE") 
      {
        addModInsKE(cur);
      }
      else if(cname == "estimator")
      { 
        if(potType =="flux")
        {
          targetH->addOperator(new ConservedEnergy,potName,false);
        }
        else if(potType == "Force")
        {
          addForceHam(cur);
        }
	else if(potType == "gofr")
        {
          PairCorrEstimator* apot=new PairCorrEstimator(*targetPtcl,sourceInp);
          apot->put(cur);
          targetH->addOperator(apot,potName,false);
        }
	else if(potType == "localmoment")
        {
          string SourceName = "ion0";
          OhmmsAttributeSet hAttrib;
          hAttrib.add(SourceName, "source");
          hAttrib.put(cur);

          PtclPoolType::iterator pit(ptclPool.find(SourceName));
          if(pit == ptclPool.end()) 
          {
            APP_ABORT("Unknown source \"" + SourceName + "\" for LocalMoment.");
          }
          ParticleSet* source = (*pit).second;

          LocalMomentEstimator* apot=new LocalMomentEstimator(*targetPtcl,*source);
          apot->put(cur);
          targetH->addOperator(apot,potName,false);
        }        
	else if(potType == "density")
        {
	  //          if(PBCType)//only if perioidic 
          {
            DensityEstimator* apot=new DensityEstimator(*targetPtcl);
            apot->put(cur);
            targetH->addOperator(apot,potName,false);
          }
        }
	else if(potType == "sk")
        {
          if(PBCType)//only if perioidic 
          {
#ifdef QMC_CUDA
            SkEstimator_CUDA* apot=new SkEstimator_CUDA(*targetPtcl);
#else
            SkEstimator* apot=new SkEstimator(*targetPtcl);
#endif
            apot->put(cur);
            targetH->addOperator(apot,potName,false);
            app_log()<<"Adding S(k) estimator"<<endl;
          }
        }
#if OHMMS_DIM==3
	else if(potType == "chiesa")
	{
	  string PsiName="psi0";
	  string SourceName = "e";
	  OhmmsAttributeSet hAttrib;
	  hAttrib.add(PsiName,"psi"); 
	  hAttrib.add(SourceName, "source");
	  hAttrib.put(cur);

	  PtclPoolType::iterator pit(ptclPool.find(SourceName));
	  if(pit == ptclPool.end()) 
	  {
            APP_ABORT("Unknown source \""+SourceName+"\" for Chiesa correction.");
	  }
	  ParticleSet &source = *pit->second;

	  OrbitalPoolType::iterator psi_it(psiPool.find(PsiName));
	  if(psi_it == psiPool.end()) 
          {
            APP_ABORT("Unknown psi \""+PsiName+"\" for Chiesa correction.");
          }

	  const TrialWaveFunction &psi = *psi_it->second->targetPsi;
	  ChiesaCorrection *chiesa = new ChiesaCorrection (source, psi);
	  targetH->addOperator(chiesa,"KEcorr",false);
	}  
#endif
        else if(potType == "Pressure")
        {
          if(estType=="coulomb")
          {
            Pressure* BP = new Pressure(*targetPtcl);
            BP-> put(cur);
            targetH->addOperator(BP,"Pressure",false);

            int nlen(100);
            attrib.add(nlen,"truncateSum");
            attrib.put(cur);
            //             DMCPressureCorr* DMCP = new DMCPressureCorr(*targetPtcl,nlen);
            //             targetH->addOperator(DMCP,"PressureSum",false);

          } 
#if defined(QMC_BUILD_COMPLETE)
	  else if (estType=="HFDHE2")
	  {
            HePressure* BP = new HePressure(*targetPtcl);
            BP-> put(cur);
            targetH->addOperator(BP,"HePress",false);
          } 
	  else if (estType=="RPAZVZB")
	  {
            RPAPressure* BP= new RPAPressure(*targetPtcl);
            
            ParticleSet* Isource;
            bool withSource=false;
            xmlNodePtr tcur = cur->children;
            while(tcur != NULL) {
              string cname((const char*)tcur->name);
              if(cname == "OneBody") 
              {
                string PsiName="psi0";
                withSource=true;
//                 string in0("ion0");
                OhmmsAttributeSet hAttrib;
//                 hAttrib.add(in0,"source");
                hAttrib.add(PsiName,"psi"); 
                hAttrib.put(tcur);
//                 renameProperty(a);
                PtclPoolType::iterator pit(ptclPool.find(sourceInp));
                if(pit == ptclPool.end()) 
		{
                  ERRORMSG("Missing source ParticleSet" << sourceInp)
                }
                Isource = (*pit).second;
                BP-> put(cur, *targetPtcl,*Isource,*(psiPool[PsiName]->targetPsi));
              }
              tcur = tcur->next; 
            }
            if (!withSource) BP-> put(cur, *targetPtcl);
            targetH->addOperator(BP,BP->MyName,false);
            
            int nlen(100);
            attrib.add(nlen,"truncateSum");
            attrib.put(cur);
//             DMCPressureCorr* DMCP = new DMCPressureCorr(*targetPtcl,nlen);
//             targetH->addOperator(DMCP,"PressureSum",false);
          }
#endif
        }
	else if(potType=="psi")
	{
	  int pwr=2;
	  OhmmsAttributeSet hAttrib;
	  hAttrib.add(pwr,"power"); 
	  hAttrib.put(cur);
	  PsiValue* PV = new PsiValue(pwr);
	  PV->put(cur,targetPtcl,ptclPool,myComm);
	  targetH->addOperator(PV,"PsiValue",false);
	}
	else if(potType=="overlap")
	{
	  int pwr=1;
	  OhmmsAttributeSet hAttrib;
	  hAttrib.add(pwr,"power"); 
	  hAttrib.put(cur);
	  
	  PsiOverlapValue* PV = new PsiOverlapValue(pwr);
	  PV->put(cur,targetPtcl,ptclPool,myComm);
	  targetH->addOperator(PV,"PsiRatio",false);
	}
	else if(potType=="DMCoverlap")
	{
	  DMCPsiValue* PV = new DMCPsiValue( );
	  PV->put(cur,targetPtcl,ptclPool,myComm);
	  targetH->addOperator(PV,"DMCPsiRatio",false);
	}
	else if(potType=="momentum")
	{
	  app_log()<<"  Adding Momentum Estimator"<<endl;
	  
	  string PsiName="psi0";
	  OhmmsAttributeSet hAttrib;
	  hAttrib.add(PsiName,"wavefunction"); 
	  hAttrib.put(cur);

	  OrbitalPoolType::iterator psi_it(psiPool.find(PsiName));
	  if(psi_it == psiPool.end()) 
          {
            APP_ABORT("Unknown psi \""+PsiName+"\" for momentum.");
          }
          TrialWaveFunction *psi=(*psi_it).second->targetPsi;
	  MomentumEstimator* ME = new MomentumEstimator(*targetPtcl, *psi);
    bool rt(myComm->rank()==0);
	  ME->putSpecial(cur,*targetPtcl,rt);
	  targetH->addOperator(ME,"MomentumEstimator",false);
	}
      } 
      else if (cname == "Kinetic")
      {
      	  string TargetName="e";
	  string SourceName = "I";
	  OhmmsAttributeSet hAttrib;
	  hAttrib.add(TargetName,"Dependant"); 
	  hAttrib.add(SourceName, "Independant");
	  hAttrib.put(cur);
      }
      if(attach2Node) xmlAddChild(myNode,xmlCopyNode(cur,1));
      cur = cur->next;
    }

    //add observables with physical and simple estimators
    int howmany=targetH->addObservables(*targetPtcl);

    //do correction
    bool dmc_correction=false;
    cur = cur_saved->children;
    while(cur != NULL) 
    {
      string cname((const char*)cur->name);
      string potType("0");
      OhmmsAttributeSet attrib;
      attrib.add(potType,"type");
      attrib.put(cur);
      if(cname == "estimator")
      {
        if(potType=="ZeroVarObs")
        {
          app_log()<<"  Adding ZeroVarObs Operator"<<endl;
          //         ZeroVarObs* FW=new ZeroVarObs();
          //         FW->put(cur,*targetH,*targetPtcl);
          //         targetH->addOperator(FW,"ZeroVarObs",false);
        }
        else if(potType == "DMCCorrection")
        {
          TrialDMCCorrection* TE = new TrialDMCCorrection();
          TE->putSpecial(cur,*targetH,*targetPtcl);
          targetH->addOperator(TE,"DMC_CORR",false);
          dmc_correction=true;
        }
        else if(potType == "ForwardWalking")
        {
          app_log()<<"  Adding Forward Walking Operator"<<endl;
          ForwardWalking* FW=new ForwardWalking();
          FW->putSpecial(cur,*targetH,*targetPtcl);
          targetH->addOperator(FW,"ForwardWalking",false);
          dmc_correction=true;
        }
      }
      cur = cur->next;
    }

    //evaluate the observables again
    if(dmc_correction) howmany=targetH->addObservables(*targetPtcl);

    return true;
  }

  void
  HamiltonianFactory::addMPCPotential(xmlNodePtr cur, bool isphysical) 
  {
#if OHMMS_DIM ==3 && defined(HAVE_LIBFFTW)
    string a("e"), title("MPC"), physical("no");
    OhmmsAttributeSet hAttrib;
    double cutoff = 30.0;
    hAttrib.add(title,"id"); 
    hAttrib.add(title,"name"); 
    hAttrib.add(cutoff,"cutoff");
    hAttrib.add(physical,"physical");
    hAttrib.put(cur);

    renameProperty(a);

    isphysical = (physical=="yes" || physical == "true");

#ifdef QMC_CUDA
    MPC_CUDA *mpc = new MPC_CUDA (*targetPtcl, cutoff);
#else
    MPC *mpc = new MPC (*targetPtcl, cutoff);
#endif
    targetH->addOperator(mpc, "MPC", isphysical);
#else
    APP_ABORT("HamiltonianFactory::addMPCPotential MPC is disabled because FFTW3 was not found during the build process.");
#endif // defined(HAVE_LIBFFTW)
  }

  void
  HamiltonianFactory::addVHXCPotential(xmlNodePtr cur) 
  {
#if OHMMS_DIM==3 && defined(HAVE_LIBFFTW)
    string a("e"), title("VHXC");
    OhmmsAttributeSet hAttrib;
    bool physical = true;
    hAttrib.add(title,"id"); 
    hAttrib.add(title,"name"); 
    hAttrib.add(physical,"physical");
    hAttrib.put(cur);

    renameProperty(a);

    VHXC *vhxc = new VHXC (*targetPtcl);
    cerr << "physical = " << physical << endl;
    targetH->addOperator(vhxc, "VHXC", physical);
#else
    APP_ABORT("HamiltonianFactory::addVHXCPotential VHXC is disabled because FFTW3 was not found during the build process.");
#endif // defined(HAVE_LIBFFTW)
  }



  void 
  HamiltonianFactory::addCoulombPotential(xmlNodePtr cur) {

    string a("e"),title("ElecElec"),pbc("yes");
    bool physical = true;
    OhmmsAttributeSet hAttrib;
    hAttrib.add(title,"id"); hAttrib.add(title,"name"); 
    hAttrib.add(a,"source"); 
    hAttrib.add(pbc,"pbc"); 
    hAttrib.add(physical,"physical");
    hAttrib.put(cur);
    

    renameProperty(a);

    PtclPoolType::iterator pit(ptclPool.find(a));
    if(pit == ptclPool.end()) 
    {
      ERRORMSG("Missing source ParticleSet" << a)
      return;
    }

    ParticleSet* source = (*pit).second;

    bool applyPBC= (PBCType && pbc=="yes");

    //CHECK PBC and create CoulombPBC for el-el
    if(source == targetPtcl) {
      if(applyPBC) 
      {
        //targetH->addOperator(new CoulombPBCAA(*targetPtcl),title);
#ifdef QMC_CUDA
	targetH->addOperator(new CoulombPBCAA_CUDA(*targetPtcl,true),title,physical);
#else
	targetH->addOperator(new CoulombPBCAATemp(*targetPtcl,true),title,physical);
#endif
      } 
      else 
      {
        if(source->getTotalNum()>1) 
#ifdef QMC_CUDA
          targetH->addOperator(new CoulombPotentialAA_CUDA(*targetPtcl), title,physical);
#else
          targetH->addOperator(new CoulombPotentialAA(*targetPtcl), title,physical);
#endif

      }
    } else {
      if(applyPBC) {
        //targetH->addOperator(new CoulombPBCAB(*source,*targetPtcl),title);
#ifdef QMC_CUDA
        targetH->addOperator(new CoulombPBCAB_CUDA(*source,*targetPtcl),title);
#else
        targetH->addOperator(new CoulombPBCABTemp(*source,*targetPtcl),title);
#endif
      } else {
        targetH->addOperator(new CoulombPotentialAB(*source,*targetPtcl),title);
      }
    }
  }

  // void
  // HamiltonianFactory::addPulayForce (xmlNodePtr cur) {
  //   string a("ion0"),targetName("e"),title("Pulay");
  //   OhmmsAttributeSet hAttrib;
  //   hAttrib.add(a,"source"); 
  //   hAttrib.add(targetName,"target"); 

  //   PtclPoolType::iterator pit(ptclPool.find(a));
  //   if(pit == ptclPool.end()) {
  //     ERRORMSG("Missing source ParticleSet" << a)
  //     return;
  //   }

  //   ParticleSet* source = (*pit).second;
  //   pit = ptclPool.find(targetName);
  //   if(pit == ptclPool.end()) {
  //     ERRORMSG("Missing target ParticleSet" << targetName)
  //     return;
  //   }
  //   ParticleSet* target = (*pit).second;
    
  //   targetH->addOperator(new PulayForce(*source, *target), title, false);

  // }

  void 
  HamiltonianFactory::addForceHam(xmlNodePtr cur) {
#if OHMMS_DIM==3
    string a("ion0"),targetName("e"),title("ForceBase"),pbc("yes"),
      PsiName="psi0";
    OhmmsAttributeSet hAttrib;
    string mode("bare");
    //hAttrib.add(title,"id");
    //hAttrib.add(title,"name"); 
    hAttrib.add(a,"source"); 
    hAttrib.add(targetName,"target"); 
    hAttrib.add(pbc,"pbc"); 
    hAttrib.add(mode,"mode"); 
    hAttrib.add(PsiName, "psi");
    hAttrib.put(cur);
    cerr << "HamFac forceBase mode " << mode << endl;
    renameProperty(a);

    PtclPoolType::iterator pit(ptclPool.find(a));
    if(pit == ptclPool.end()) {
      ERRORMSG("Missing source ParticleSet" << a)
      return;
    }
    ParticleSet* source = (*pit).second;
    pit = ptclPool.find(targetName);
    if(pit == ptclPool.end()) {
      ERRORMSG("Missing target ParticleSet" << targetName)
      return;
    }
    ParticleSet* target = (*pit).second;

    //bool applyPBC= (PBCType && pbc=="yes");
    if(mode=="bare") 
      targetH->addOperator(new BareForce(*source, *target), title, false);
    else if(mode=="cep") 
      targetH->addOperator(new ForceCeperley(*source, *target), title, false);
    else if(mode=="pulay") {
      OrbitalPoolType::iterator psi_it(psiPool.find(PsiName));
      if(psi_it == psiPool.end()) {
	APP_ABORT("Unknown psi \""+PsiName+"\" for Pulay force.");
      }
      TrialWaveFunction &psi = *psi_it->second->targetPsi;
      targetH->addOperator(new PulayForce(*source, *target, psi), 
			   "PulayForce", false);
    }
    else if(mode=="zero_variance") {
      app_log() << "Adding zero-variance force term.\n";
      OrbitalPoolType::iterator psi_it(psiPool.find(PsiName));
      if(psi_it == psiPool.end()) {
	APP_ABORT("Unknown psi \""+PsiName+"\" for zero-variance force.");
      }
      TrialWaveFunction &psi = *psi_it->second->targetPsi;
      targetH->addOperator
	(new ZeroVarianceForce(*source, *target, psi), "ZVForce", false);
    }
    else {
      ERRORMSG("Failed to recognize Force mode " << mode);
      //} else if(mode=="FD") {
      //  targetH->addOperator(new ForceFiniteDiff(*source, *target), title, false);
    }
#endif
  }

  void 
  HamiltonianFactory::addPseudoPotential(xmlNodePtr cur) 
  {

#if OHMMS_DIM == 3
    string src("i"),title("PseudoPot"),wfname("invalid"),format("xml");

    OhmmsAttributeSet pAttrib;
    pAttrib.add(title,"name");
    pAttrib.add(src,"source");
    pAttrib.add(wfname,"wavefunction");
    pAttrib.add(format,"format"); //temperary tag to switch between format
    pAttrib.put(cur);

    if(format == "old")
    {
      APP_ABORT("pseudopotential Table format is not supported.");
    }

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

    app_log() << endl << "  ECPotential builder for pseudopotential "<< endl;
    ECPotentialBuilder ecp(*targetH,*ion,*targetPtcl,*psi,myComm);
    ecp.put(cur);
#else
    APP_ABORT("HamiltonianFactory::addPseudoPotential\n pairpot@type=\"pseudo\" is invalid if DIM != 3");
#endif
  }

  void 
  HamiltonianFactory::addCorePolPotential(xmlNodePtr cur) 
  {
#if OHMMS_DIM == 3
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
#else
    APP_ABORT("HamiltonianFactory::addCorePolPotential\n pairpot@type=\"cpp\" is invalid if DIM != 3");
#endif
  }

  void 
  HamiltonianFactory::addConstCoulombPotential(xmlNodePtr cur, string& nuclei)
  {
    OhmmsAttributeSet hAttrib;
    string hname("IonIon");
    string forces("no");
    hAttrib.add(forces,"forces");
    hAttrib.add(hname,"name");
    hAttrib.put(cur);
    bool doForces = (forces == "yes") || (forces == "true");

    app_log() << "  Creating Coulomb potential " << nuclei << "-" << nuclei << endl;
    renameProperty(nuclei);
    PtclPoolType::iterator pit(ptclPool.find(nuclei));
    if(pit != ptclPool.end()) {
      ParticleSet* ion=(*pit).second;
      if(PBCType)
      {
#ifdef QMC_CUDA
	targetH->addOperator(new CoulombPBCAA_CUDA(*ion,false,doForces),hname);
#else
	targetH->addOperator(new CoulombPBCAATemp(*ion,false,doForces),hname);
#endif
      } else {
        if(ion->getTotalNum()>1) 
          targetH->addOperator(new IonIonPotential(*ion),hname);
      }
    }
  }

  void
  HamiltonianFactory::addModInsKE(xmlNodePtr cur) {
#if defined(HAVE_LIBFFTW_LS)
    typedef QMCTraits::RealType    RealType;
    typedef QMCTraits::IndexType   IndexType;
    typedef QMCTraits::PosType     PosType;
    
    string Dimensions, DispRelType, PtclSelType, MomDistType;
    RealType Cutoff, GapSize(0.0), FermiMomentum(0.0);
    
    OhmmsAttributeSet pAttrib;
    pAttrib.add(Dimensions, "dims");
    pAttrib.add(DispRelType, "dispersion");
    pAttrib.add(PtclSelType, "selectParticle");
    pAttrib.add(Cutoff, "cutoff");
    pAttrib.add(GapSize, "gapSize");
    pAttrib.add(FermiMomentum, "kf");
    pAttrib.add(MomDistType, "momdisttype");
    pAttrib.put(cur);
    
    if (MomDistType == "") MomDistType = "FFT";
    
    TrialWaveFunction* psi;
    psi = (*(psiPool.begin())).second->targetPsi;
    
    Vector<PosType> LocLattice;
    Vector<IndexType> DimSizes;
    Vector<RealType> Dispersion;

    if (Dimensions == "3") {
      gen3DLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
    } else if (Dimensions == "1" || Dimensions == "1averaged") {
      gen1DLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
    } else if (Dimensions == "homogeneous") {
      genDegenLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
    } else {
      ERRORMSG("Dimensions value not recognized!")
    }
    
    if (DispRelType == "freeParticle") {
      genFreeParticleDispersion(LocLattice, Dispersion);
    } else if (DispRelType == "simpleModel") {
      genSimpleModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
    } else if (DispRelType == "pennModel") {
      genPennModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
    } else if (DispRelType == "debug") {
      genDebugDispersion(LocLattice, Dispersion);  
    } else {
      ERRORMSG("Dispersion relation not recognized");
    }
    
    PtclChoiceBase* pcp;
    if (PtclSelType == "random") {
      pcp = new RandomChoice(*targetPtcl);
    } else if (PtclSelType == "randomPerWalker") {
      pcp = new RandomChoicePerWalker(*targetPtcl);
    } else if (PtclSelType == "constant") {
      pcp = new StaticChoice(*targetPtcl);
    } else {
      ERRORMSG("Particle choice policy not recognized!");
    }
    
    MomDistBase* mdp;
    if (MomDistType == "direct") {
      mdp = new RandomMomDist(*targetPtcl, LocLattice, pcp);
    } else if (MomDistType == "FFT" || MomDistType =="fft") { 
      if (Dimensions == "3") {
	mdp = new ThreeDimMomDist(*targetPtcl, DimSizes, pcp);
      } else if (Dimensions == "1") {
	mdp = new OneDimMomDist(*targetPtcl, DimSizes, pcp);
      } else if (Dimensions == "1averaged") {
	mdp = new AveragedOneDimMomDist(*targetPtcl, DimSizes, pcp);
      } else {
	ERRORMSG("Dimensions value not recognized!");
      }
    } else {
      ERRORMSG("MomDistType value not recognized!");
    }
    delete pcp;
    
    QMCHamiltonianBase* modInsKE = new ModInsKineticEnergy(*psi, Dispersion, mdp);
    modInsKE->put(cur);
    targetH->addOperator(modInsKE, "ModelInsKE");
    
    delete mdp;
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
    HamiltonianFactory* aCopy=new HamiltonianFactory(qp, ptclPool, psiPool, myComm);
    aCopy->setName(aname);

    aCopy->renameProperty("e",qp->getName());
    aCopy->renameProperty(psiName,psi->getName());
    aCopy->build(myNode,false);
    myClones[ip]=aCopy;
    //aCopy->get(app_log());
    return aCopy;
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
