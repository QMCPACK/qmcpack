//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    D.C. Yang, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file HamiltonianFactory.cpp
 *@brief Definition of a HamiltonianFactory
 */
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCHamiltonians/BareKineticEnergy.h"
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCHamiltonians/SpeciesKineticEnergy.h"
#include "QMCHamiltonians/LatticeDeviationEstimator.h"
#include "QMCHamiltonians/NumericalRadialPotential.h"
#include "QMCHamiltonians/MomentumEstimator.h"
#include "QMCHamiltonians/Pressure.h"
#include "QMCHamiltonians/ForwardWalking.h"
#include "QMCHamiltonians/NumberFluctuations.h"
#include "QMCHamiltonians/PairCorrEstimator.h"
#include "QMCHamiltonians/LocalMomentEstimator.h"
#include "QMCHamiltonians/DensityEstimator.h"
#include "QMCHamiltonians/SkEstimator.h"
#include "QMCHamiltonians/model/HarmonicExternalPotential.h"
#include "QMCHamiltonians/StaticStructureFactor.h"
#include "QMCHamiltonians/SpinDensity.h"
#include "QMCHamiltonians/OrbitalImages.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "QMCHamiltonians/EnergyDensityEstimator.h"
#include "QMCHamiltonians/NearestNeighborsEstimator.h"
#include "QMCHamiltonians/DensityMatrices1B.h"
#endif
#if OHMMS_DIM == 3
#include "QMCHamiltonians/ChiesaCorrection.h"
#include "QMCHamiltonians/SkAllEstimator.h"
#if defined(HAVE_LIBFFTW_LS)
#include "QMCHamiltonians/model/ModInsKineticEnergy.h"
#include "QMCHamiltonians/MomentumDistribution.h"
#include "QMCHamiltonians/DispersionRelation.h"
#endif
#endif
// #include "QMCHamiltonians/ZeroVarObs.h"
#if !defined(QMC_CUDA) && QMC_BUILD_LEVEL>2
#include "QMCHamiltonians/SkPot.h"
#include "QMCHamiltonians/model/HardSphere.h"
#include "QMCHamiltonians/model/GaussianPot.h"
#include "QMCHamiltonians/model/HusePot.h"
#include "QMCHamiltonians/model/OscillatoryPot.h"
#include "QMCHamiltonians/model/ModPosTelPot.h"
#include "QMCHamiltonians/model/HFDHE2Potential_tail.h"
#include "QMCHamiltonians/model/HePressure.h"
#include "QMCHamiltonians/model/JelliumPotential.h"
#include "QMCHamiltonians/model/HFDHE2Potential.h"
#include "QMCHamiltonians/model/HeEPotential.h"
#include "QMCHamiltonians/model/HeEPotential_tail.h"
#include "QMCHamiltonians/model/LennardJones_smoothed.h"
#include "QMCHamiltonians/model/HFDHE2_Moroni1995.h"
//#include "QMCHamiltonians/HFDBHE_smoothed.h"
#include "QMCHamiltonians/model/HeSAPT_smoothed.h"
#endif
#include "OhmmsData/AttributeSet.h"
#ifdef QMC_CUDA
#include "QMCHamiltonians/SkEstimator_CUDA.h"
#endif

namespace qmcplusplus
{
HamiltonianFactory::HamiltonianFactory(ParticleSet* qp,
                                       PtclPoolType& pset, OrbitalPoolType& oset, Communicate* c)
  : MPIObjectBase(c), targetPtcl(qp), targetH(0)
  , ptclPool(pset),psiPool(oset), myNode(NULL), psiName("psi0")
{
  //PBCType is zero or 1 but should be generalized
  PBCType=targetPtcl->Lattice.SuperCellEnum;
  ClassName="HamiltonianFactory";
  myName="psi0";
  targetPtcl->set_quantum();
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
bool HamiltonianFactory::build(xmlNodePtr cur, bool buildtree)
{
  if(cur == NULL)
    return false;
  std::string htype("generic"), source("i"), defaultKE("yes");
  OhmmsAttributeSet hAttrib;
  hAttrib.add(htype,"type");
  hAttrib.add(source,"source");
  hAttrib.add(defaultKE,"default");
  hAttrib.put(cur);
  renameProperty(source);
  bool attach2Node=false;
  if(buildtree)
  {
    if(myNode == NULL)
    {
//#if (LIBXMLD_VERSION < 20616)
//        app_warning() << "   Workaround of libxml2 bug prior to 2.6.x versions" << std::endl;
//        myNode = xmlCopyNode(cur,2);
//#else
//        app_warning() << "   using libxml2 2.6.x versions" << std::endl;
//        myNode = xmlCopyNode(cur,1);
//#endif
      myNode = xmlCopyNode(cur,1);
    }
    else
    {
      attach2Node=true;
    }
  }
  if(targetH==0)
  {
    targetH  = new QMCHamiltonian;
    targetH->setName(myName);
    targetH->addOperator(new BareKineticEnergy<OHMMS_PRECISION_FULL>(*targetPtcl),"Kinetic");
  }
  OrbitalPoolType::iterator psi_it(psiPool.find(psiName));
  if(psi_it == psiPool.end())
    APP_ABORT("Unknown psi \""+psiName+"\" for target Psi");
  TrialWaveFunction* targetPsi = psi_it->second->targetPsi;
  xmlNodePtr cur_saved(cur);
  cur = cur->children;
  while(cur != NULL)
  {
    std::string notype = "0";
    std::string noname = "any";
    std::string cname((const char*)cur->name);
    std::string potType(notype);
    std::string potName(noname);
    std::string potUnit("hartree");
    std::string estType("coulomb");
    std::string sourceInp(targetPtcl->getName());
    std::string targetInp(targetPtcl->getName());
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

    int nham=targetH->total_size();
    if(cname == "pairpot")
    {
      if(potType == "coulomb")
      {
        addCoulombPotential(cur);
      }
#if !defined(QMC_CUDA) && QMC_BUILD_LEVEL>2
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
      else if (potType == "huse")
      {
        HusePot* hs = new HusePot(*targetPtcl);
        hs->put(cur);
        targetH->addOperator(hs,"HusePot",true);
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
        addPseudoPotential(cur);
#if !defined(QMC_CUDA) && QMC_BUILD_LEVEL>2
      else if(potType == "cpp")
      {
        addCorePolPotential(cur);
      }
      else if (potType == "LJP_smoothed")
      {
        LennardJones_smoothed_phy* LJP = new LennardJones_smoothed_phy(*targetPtcl);
        targetH->addOperator(LJP,"LJP",true);
        LJP->addCorrection(*targetH);
      }
      else if (potType == "HeSAPT_smoothed")
      {
        HeSAPT_smoothed_phy* SAPT = new HeSAPT_smoothed_phy(*targetPtcl);
        targetH->addOperator(SAPT,"HeSAPT",true);
        SAPT->addCorrection(*targetH);
      }
      else if (potType == "HFDHE2_Moroni1995")
      {
        HFDHE2_Moroni1995_phy* HFD = new HFDHE2_Moroni1995_phy(*targetPtcl);
        targetH->addOperator(HFD,"HFD-HE2",true);
        HFD->addCorrection(*targetH);
      }
      else if(potType == "eHe")
      {
        std::string SourceName = "e";
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
        std::string SourceName = "e";
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
        app_log() << "  Adding HFDHE2Potential(Au) " << std::endl;
      }
      else if(cname == "modInsKE")
      {
        addModInsKE(cur);
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
    {
      //just to support old input
      if(potType == "coulomb")
        addCoulombPotential(cur);
    }
    else if(cname == "extpot" )
    {
      if (potType == "harmonic_ext" || potType=="HarmonicExt")
      {
        HarmonicExternalPotential* hs = new HarmonicExternalPotential(*targetPtcl);
        hs->put(cur);
        targetH->addOperator(hs,"HarmonicExt",true);
      }
    }
    else if(cname == "estimator")
    {
      if(potType =="flux")
      {
        targetH->addOperator(new ConservedEnergy,potName,false);
      }
      else if(potType =="specieskinetic")
      {
        SpeciesKineticEnergy* apot = new SpeciesKineticEnergy(*targetPtcl);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
      else if(potType =="latticedeviation")
      {
        // find target particle set
        PtclPoolType::iterator pit(ptclPool.find(targetInp));
        if(pit == ptclPool.end())
        {
          APP_ABORT("Unknown target \"" + targetInp + "\" for LatticeDeviation.");
        }
        ParticleSet* target_particle_set = (*pit).second;

        // find source particle set
        PtclPoolType::iterator spit(ptclPool.find(sourceInp));
        if(spit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + sourceInp + "\" for LatticeDeviation.");
        }
        ParticleSet* source_particle_set = (*spit).second;

        // read xml node
        OhmmsAttributeSet local_attrib;
        std::string target_group,source_group;
        local_attrib.add(target_group,"tgroup");
        local_attrib.add(source_group,"sgroup");
        local_attrib.put(cur);

        LatticeDeviationEstimator* apot = new LatticeDeviationEstimator(*target_particle_set
          ,*source_particle_set,target_group,source_group);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
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
        std::string SourceName = "ion0";
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
      else if(potType == "numberfluctuations")
      {
        app_log()<<" Adding Number Fluctuation estimator"<< std::endl;
        NumberFluctuations* apot=new NumberFluctuations(*targetPtcl);
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
      else if(potType == "spindensity")
      {
        app_log()<<"  Adding SpinDensity"<< std::endl;
        SpinDensity* apot=new SpinDensity(*targetPtcl);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
      else if(potType == "structurefactor")
      {
        app_log()<<"  Adding StaticStructureFactor"<< std::endl;
        StaticStructureFactor* apot=new StaticStructureFactor(*targetPtcl);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
      else if(potType == "orbitalimages")
      {
        app_log()<<"  Adding OrbitalImages"<< std::endl;
        OrbitalImages* apot=new OrbitalImages(*targetPtcl,ptclPool,myComm);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
#if !defined(REMOVE_TRACEMANAGER)
      else if(potType == "energydensity" || potType == "EnergyDensity")
      {
        app_log()<<"  Adding EnergyDensityEstimator"<< std::endl;
        EnergyDensityEstimator* apot=new EnergyDensityEstimator(ptclPool,defaultKE);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
      else if(potType == "nearestneighbors" || potType == "NearestNeighbors")
      {
        app_log()<<"  Adding NearestNeighborsEstimator"<< std::endl;
        NearestNeighborsEstimator* apot=new NearestNeighborsEstimator(ptclPool);
        apot->put(cur);
        targetH->addOperator(apot,"nearest_neighbors",false);
      }
      else if(potType == "dm1b")
      {
        app_log()<<"  Adding DensityMatrices1B"<< std::endl;
        std::string source = "";
        OhmmsAttributeSet attrib;
        attrib.add(source,"source");
        attrib.put(cur);
        PtclPoolType::iterator pit(ptclPool.find(source));
        ParticleSet* Pc;
        if(source=="")
          Pc = NULL;
        else if(pit!=ptclPool.end())
          Pc = pit->second;
        else
        {
          APP_ABORT("Unknown source \"" + source + "\" for DensityMatrices1B");
        }
        DensityMatrices1B* apot=new DensityMatrices1B(*targetPtcl,*targetPsi,Pc);
        apot->put(cur);
        targetH->addOperator(apot,potName,false);
      }
#endif
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
          app_log()<<"Adding S(k) estimator"<< std::endl;
#if defined(USE_REAL_STRUCT_FACTOR)
          app_log()<<"S(k) estimator using Real S(k)"<< std::endl;
#endif
        }
      }
#if OHMMS_DIM==3
      else if(potType == "chiesa")
      {
#ifdef ENABLE_SOA
        app_warning() << "Skip Chiesa estimator due to the lack of support with SoA."
                      << " Access the correction via AoS at the moment." << std::endl;
#else
        std::string PsiName="psi0";
        std::string SourceName = "e";
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
#endif
      }
      else if(potType == "skall")
      {
        std::string SourceName = "";
        OhmmsAttributeSet attrib;
        attrib.add(SourceName,"source");
        attrib.put(cur);

        PtclPoolType::iterator pit(ptclPool.find(SourceName));
        if(pit == ptclPool.end())
        {
          APP_ABORT("Unknown source \"" + SourceName + "\" for LocalMoment.");
        }
        ParticleSet* source = (*pit).second;

        if(PBCType)
        {

          SkAllEstimator* apot=new SkAllEstimator(*source, *targetPtcl);
          apot->put(cur);
          targetH->addOperator(apot,potName,false);
          app_log()<<"Adding S(k) ALL estimator"<<std::endl;
        }
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
#endif
      }
      else if(potType=="momentum")
      {
        app_log()<<"  Adding Momentum Estimator"<< std::endl;
        std::string PsiName="psi0";
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
      std::string TargetName="e";
      std::string SourceName = "I";
      OhmmsAttributeSet hAttrib;
      hAttrib.add(TargetName,"Dependant");
      hAttrib.add(SourceName, "Independant");
      hAttrib.put(cur);
    }

    if(nham<targetH->total_size()) //if(cname!="text" && cname !="comment")
    {
      if(potName==noname) 
      {
        potName=potType;
        app_log() << "Provide name for hamiltonian element for type " << potType << std::endl;
      }
      //APP_ABORT("HamiltonianFactory::build\n  a name for operator of type "+cname+" "+potType+" must be provided in the xml input");
      targetH->addOperatorType(potName,potType);
    }

    if(attach2Node)
      xmlAddChild(myNode,xmlCopyNode(cur,1));
    cur = cur->next;
  }
  //add observables with physical and simple estimators
  int howmany=targetH->addObservables(*targetPtcl);
  //do correction
  bool dmc_correction=false;
  cur = cur_saved->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    std::string potType("0");
    OhmmsAttributeSet attrib;
    attrib.add(potType,"type");
    attrib.put(cur);
    if(cname == "estimator")
    {
      if(potType=="ZeroVarObs")
      {
        app_log()<<"  Not Adding ZeroVarObs Operator"<< std::endl;
        //         ZeroVarObs* FW=new ZeroVarObs();
        //         FW->put(cur,*targetH,*targetPtcl);
        //         targetH->addOperator(FW,"ZeroVarObs",false);
      }
      //         else if(potType == "DMCCorrection")
      //         {
      //           TrialDMCCorrection* TE = new TrialDMCCorrection();
      //           TE->putSpecial(cur,*targetH,*targetPtcl);
      //           targetH->addOperator(TE,"DMC_CORR",false);
      //           dmc_correction=true;
      //         }
      else if(potType == "ForwardWalking")
      {
        app_log()<<"  Adding Forward Walking Operator"<< std::endl;
        ForwardWalking* FW=new ForwardWalking();
        FW->putSpecial(cur,*targetH,*targetPtcl);
        targetH->addOperator(FW,"ForwardWalking",false);
        dmc_correction=true;
      }
    }
    cur = cur->next;
  }
  //evaluate the observables again
  if(dmc_correction)
    howmany=targetH->addObservables(*targetPtcl);
  return true;
}

void
HamiltonianFactory::addModInsKE(xmlNodePtr cur)
{
#if defined(HAVE_LIBFFTW_LS)
  typedef QMCTraits::RealType    RealType;
  typedef QMCTraits::IndexType   IndexType;
  typedef QMCTraits::PosType     PosType;
  std::string Dimensions, DispRelType, PtclSelType, MomDistType;
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
  if (MomDistType == "")
    MomDistType = "FFT";
  TrialWaveFunction* psi;
  psi = (*(psiPool.begin())).second->targetPsi;
  Vector<PosType> LocLattice;
  Vector<IndexType> DimSizes;
  Vector<RealType> Dispersion;
  if (Dimensions == "3")
  {
    gen3DLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
  }
  else if (Dimensions == "1" || Dimensions == "1averaged")
  {
    gen1DLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
  }
  else if (Dimensions == "homogeneous")
  {
    genDegenLattice(Cutoff, *targetPtcl, LocLattice, Dispersion, DimSizes);
  }
  else
  {
    ERRORMSG("Dimensions value not recognized!")
  }
  if (DispRelType == "freeParticle")
  {
    genFreeParticleDispersion(LocLattice, Dispersion);
  }
  else if (DispRelType == "simpleModel")
  {
    genSimpleModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
  }
  else if (DispRelType == "pennModel")
  {
    genPennModelDispersion(LocLattice, Dispersion, GapSize, FermiMomentum);
  }
  else if (DispRelType == "debug")
  {
    genDebugDispersion(LocLattice, Dispersion);
  }
  else
  {
    ERRORMSG("Dispersion relation not recognized");
  }
  PtclChoiceBase* pcp;
  if (PtclSelType == "random")
  {
    pcp = new RandomChoice(*targetPtcl);
  }
  else if (PtclSelType == "randomPerWalker")
  {
    pcp = new RandomChoicePerWalker(*targetPtcl);
  }
  else if (PtclSelType == "constant")
  {
    pcp = new StaticChoice(*targetPtcl);
  }
  else
  {
    ERRORMSG("Particle choice policy not recognized!");
  }
  MomDistBase* mdp;
  if (MomDistType == "direct")
  {
    mdp = new RandomMomDist(*targetPtcl, LocLattice, pcp);
  }
  else if (MomDistType == "FFT" || MomDistType =="fft")
  {
    if (Dimensions == "3")
    {
      mdp = new ThreeDimMomDist(*targetPtcl, DimSizes, pcp);
    }
    else if (Dimensions == "1")
    {
      mdp = new OneDimMomDist(*targetPtcl, DimSizes, pcp);
    }
    else if (Dimensions == "1averaged")
    {
      mdp = new AveragedOneDimMomDist(*targetPtcl, DimSizes, pcp);
    }
    else
    {
      ERRORMSG("Dimensions value not recognized!");
    }
  }
  else
  {
    ERRORMSG("MomDistType value not recognized!");
  }
  delete pcp;
  QMCHamiltonianBase* modInsKE = new ModInsKineticEnergy(*psi, Dispersion, mdp);
  modInsKE->put(cur);
  targetH->addOperator(modInsKE, "ModelInsKE");
  delete mdp;
#else
  app_error() << "  ModelInsulatorKE cannot be used without FFTW " << std::endl;
#endif
}

void HamiltonianFactory::renameProperty(const std::string& a, const std::string& b)
{
  RenamedProperty[a]=b;
}

void HamiltonianFactory::setCloneSize(int np)
{
  myClones.resize(np,0);
}

//TrialWaveFunction*
//HamiltonianFactory::cloneWaveFunction(ParticleSet* qp, int ip) {
//  HamiltonianFactory* aCopy= new HamiltonianFactory(qp,ptclPool);
//  aCopy->put(myNode,false);
//  myClones[ip]=aCopy;
//  return aCopy->targetPsi;
//}

void HamiltonianFactory::renameProperty( std::string& aname)
{
  std::map<std::string,std::string>::iterator it(RenamedProperty.find(aname));
  if(it != RenamedProperty.end())
  {
    aname=(*it).second;
  }
}
HamiltonianFactory::~HamiltonianFactory()
{
  //clean up everything
}

HamiltonianFactory*
HamiltonianFactory::clone(ParticleSet* qp, TrialWaveFunction* psi,
                          int ip, const std::string& aname)
{
  HamiltonianFactory* aCopy=new HamiltonianFactory(qp, ptclPool, psiPool, myComm);
  aCopy->setName(aname);
  aCopy->renameProperty("e",qp->getName());
  aCopy->renameProperty(psiName,psi->getName());
  aCopy->build(myNode,false);
  myClones[ip]=aCopy;
  //aCopy->get(app_log());
  return aCopy;
}

bool HamiltonianFactory::put(xmlNodePtr cur)
{
  bool success=build(cur,false);
  
  //if(targetH->EnableVirtualMoves)
  //{
  //  app_log() << "  Hamiltonian enables VirtualMoves" << std::endl;
  //  targetPtcl->enableVirtualMoves();
  //}
  //else
  //  app_log() << "  Hamiltonian disables VirtualMoves" << std::endl;
  
  return success;
}

void HamiltonianFactory::reset() { }
}
