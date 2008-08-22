//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/DMC/WalkerControlFactory.h"
#include <numeric>
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"
#include "QMCDrivers/WalkerControlBase.h"
#include "Estimators/EstimatorManager.h"
//#include "Estimators/DMCEnergyEstimator.h"
#include "QMCDrivers/BranchIO.h"
#include "HDFVersion.h"

//#include <boost/archive/text_oarchive.hpp>

namespace qmcplusplus 
{

  ///enum to yes/no options saved in sParam 
  enum {COMBOPT, USETAUOPT, DUMMYOPT};

  SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
    vParam(1.0), WalkerController(0), MyEstimator(0), PopHist(5), DMCEnergyHist(5)
  {

    BranchMode.set(B_DMCSTAGE,0); //warmup stage
    BranchMode.set(B_POPCONTROL,1); //use standard DMC
    BranchMode.set(B_USETAUEFF,1); //use taueff
    BranchMode.set(B_CLEARHISTORY,0); //clear history and start with the current average

    vParam[B_TAU]=tau;
    vParam[B_TAUEFF]=tau;
    Feedback=1.0;
    R2Accepted(1.0e-10); 
    R2Proposed(1.0e-10);

    //set the default values for integer parameters
    iParam[B_WARMUPSTEPS]=1000;
    iParam[B_ENERGYUPDATEINTERVAL]=1;
    iParam[B_BRANCHINTERVAL]=1;
    iParam[B_TARGETWALKERS]=nideal;
    iParam[B_MAXWALKERS]=nideal;
    iParam[B_MINWALKERS]=nideal;
    iParam[B_COUNTER]=-1;

    //default is no
    sParam.resize(DUMMYOPT,"no");

    registerParameters();
    reset();
  }

  /** copy constructor
   *
   * Copy only selected data members and WalkerController is never copied.
   */
  SimpleFixedNodeBranch::SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch):
    BranchMode(abranch.BranchMode),
    iParam(abranch.iParam),
    vParam(abranch.vParam),
    Feedback(abranch.Feedback),
    WalkerController(0), MyEstimator(0),
    sParam(abranch.sParam)
  {
    registerParameters();
    reset();
  }

  void SimpleFixedNodeBranch::registerParameters() 
  {
    m_param.add(iParam[B_WARMUPSTEPS],"warmupSteps","int"); 
    m_param.add(iParam[B_ENERGYUPDATEINTERVAL],"energyUpdateInterval","int"); 
    m_param.add(iParam[B_BRANCHINTERVAL],"branchInterval","int"); 

    m_param.add(iParam[B_TARGETWALKERS],"targetWalkers","int"); 
    m_param.add(iParam[B_TARGETWALKERS],"targetwalkers","int"); 
    m_param.add(iParam[B_TARGETWALKERS],"target_walkers","int"); 

    //trial energy
    m_param.add(vParam[B_EREF],"refEnergy","AU"); 
    m_param.add(vParam[B_EREF],"ref_energy","AU"); 
    m_param.add(vParam[B_EREF],"en_ref","AU");

    //feed back parameter for population control
    m_param.add(Feedback,"feedback","double"); 

    //turn on/off effective tau onl for time-step error comparisons
    m_param.add(sParam[USETAUOPT],"useBareTau","option");
  }

  void SimpleFixedNodeBranch::start(const string& froot, bool append) 
  {
    RootName=froot;
    MyEstimator->RootName=froot;
    MyEstimator->reset();
  }

  void SimpleFixedNodeBranch::initWalkerController(RealType tau, bool fixW) 
  {

    vParam[B_TAU]=tau;
    if(!BranchMode[B_DMCSTAGE])
      vParam[B_TAUEFF]=tau*R2Accepted.result()/R2Proposed.result();

    if(WalkerController == 0) 
    {
      BranchMode.set(B_DMC,1);//set DMC
      BranchMode.set(B_POPCONTROL,!fixW);//fixW -> 0 
      WalkerController = createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode);
      iParam[B_MAXWALKERS]=WalkerController->Nmax;
      iParam[B_MINWALKERS]=WalkerController->Nmin;
      WalkerController->start();

      PopHist.clear();
      PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],5));
    }

    //save the BranchMode in anticipating state changes in reset
    bitset<B_MODE_MAX> bmode(BranchMode);
    //reset Feedback pararmeter
    this->reset();

    MyEstimator->reset();
    //update the simulation parameters
    WalkerController->put(myNode);
    //assign current Eref and a large number for variance
    WalkerController->setEnergyAndVariance(vParam[B_EREF],vParam[B_SIGMA]);

    //determine the branch cutoff to limit wild weights based on the sigma and sigmaBound
    RealType sigma=std::max(std::sqrt(static_cast<RealType>(iParam[B_TARGETWALKERS]))*vParam[B_SIGMA]*WalkerController->targetSigma,100.0);
    vParam[B_BRANCHCUTOFF]=std::min(sigma,5.0/tau);
    //vParam[B_BRANCHCUTOFF]=vParam[B_SIGMA]*WalkerController->targetSigma;
    vParam[B_BRANCHMAX]=vParam[B_BRANCHCUTOFF]*1.5;
    vParam[B_BRANCHFILTER]=1.0/(vParam[B_BRANCHMAX]-vParam[B_BRANCHCUTOFF]);

    //reset controller 
    WalkerController->reset();
    BranchMode=bmode;

    app_log() << "  QMC counter      = " << iParam[B_COUNTER] << endl;
    app_log() << "  time step        = " << vParam[B_TAU] << endl;
    app_log() << "  effective time step = " << vParam[B_TAUEFF] << endl;
    app_log() << "  trial energy     = " << vParam[B_ETRIAL] << endl;
    app_log() << "  reference energy = " << vParam[B_EREF] << endl;
    app_log() << "  Feedback = " << Feedback <<  endl;
    app_log() << "  reference variance = " << vParam[B_SIGMA] << endl;
    app_log() << "  branch cutoff = " <<  vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << endl;
    app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << endl;
    app_log() << "  Max and mimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS] << endl;
    app_log() << "  QMC Status (BranchMode) = " << BranchMode << endl;
  }

  void SimpleFixedNodeBranch::flush(int counter) 
  {
    if(counter==0 && WalkerController) WalkerController->reset();
  }

  void SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& Walkers) 
  {
    //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance
    RealType pop_now = WalkerController->branch(iter,Walkers,0.1);

    //current energy
    vParam[B_ENOW]=WalkerController->EnsembleProperty.Energy;
    EnergyHist(vParam[B_ENOW]);
    R2Accepted(WalkerController->EnsembleProperty.R2Accepted);
    R2Proposed(WalkerController->EnsembleProperty.R2Proposed);
    //PopHist(pop_now);
    vParam[B_EREF]=EnergyHist.mean();//current mean

    if(BranchMode[B_DMCSTAGE]) // main stage
    { 
      if(BranchMode[B_POPCONTROL])
      {
        if(ToDoSteps>0)
          --ToDoSteps;
        else
        {
          vParam[B_ETRIAL]=vParam[B_EREF]+Feedback*(logN-std::log(pop_now)); 
          ToDoSteps=iParam[B_ENERGYUPDATEINTERVAL]-1;
        }
      }
    }
    else//warmup
    {
      if(BranchMode[B_USETAUEFF]) vParam[B_TAUEFF]=vParam[B_TAU]*R2Accepted.result()/R2Proposed.result();
      if(BranchMode[B_POPCONTROL])
      {
        //using enow for the first 100 step
        //RealType emix=(ToDoSteps>900)?vParam[B_ENOW]:vParam[B_EREF];
        //vParam[B_ETRIAL]=emix+Feedback*(logN-std::log(pop_now));
        
        //using eref
        vParam[B_ETRIAL]=vParam[B_EREF]+Feedback*(logN-std::log(pop_now));
      }
      --ToDoSteps;
      if(ToDoSteps==0)  //warmup is done
      {

        RealType sigma_eq=std::sqrt(iParam[B_TARGETWALKERS]*EnergyHist.variance());
        RealType sigma=std::max(sigma_eq*WalkerController->targetSigma,10.0);
        vParam[B_BRANCHCUTOFF]=std::min(sigma,5.0/vParam[B_TAU]);
        vParam[B_BRANCHMAX]=vParam[B_BRANCHCUTOFF]*1.5;
        vParam[B_BRANCHFILTER]=1.0/(vParam[B_BRANCHMAX]-vParam[B_BRANCHCUTOFF]);


        app_log() << "\n Warmup is completed after " << iParam[B_WARMUPSTEPS] << endl;
        if(BranchMode[B_USETAUEFF])
          app_log() << "\n  TauEff     = " << vParam[B_TAUEFF] << "\n TauEff/Tau = " << vParam[B_TAUEFF]/vParam[B_TAU];
        else
          app_log() << "\n  TauEff proposed   = " << vParam[B_TAUEFF]*R2Accepted.result()/R2Proposed.result();

        app_log()  << "\n  Etrial     = " << vParam[B_ETRIAL] << endl;
        app_log() << " Running average of energy = " << EnergyHist.mean() << endl;
        //app_log() << "            sqrt(variance) = " << std::sqrt(EnergyHist.variance()) << endl;
        app_log() << "            sqrt(variance) = " << sigma_eq << endl;
        app_log() << "branch cutoff = " <<  vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << endl;

        ToDoSteps = iParam[B_ENERGYUPDATEINTERVAL]-1;
        BranchMode.set(B_DMCSTAGE,1); //set BranchModex to main stage

        //reset the histogram
        EnergyHist.clear();
        EnergyHist(vParam[B_EREF]);
       
        //This is not necessary
        //EnergyHist(DMCEnergyHist.mean());
      }
    }

    WalkerController->setTrialEnergy(vParam[B_ETRIAL]);

    //evaluate everything else
    MyEstimator->accumulate(Walkers);

  }

  /** perform branching
   *
   * Set the trial energy of clones
   */
  void SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones) 
  {
    branch(iter,w);
    //synchronize it
    for(int i=0; i<clones.size(); i++) clones[i]->vParam=vParam;
  }

  void SimpleFixedNodeBranch::reset() 
  {
    //use effective time step of BranchInterval*Tau
    //Feed = 1.0/(static_cast<RealType>(NumGeneration*BranchInterval)*Tau);
    //logN = Feed*std::log(static_cast<RealType>(Nideal));
    if(WalkerController)
    {
      //this is to compare the time step errors
      BranchMode.set(B_USETAUEFF,sParam[USETAUOPT]=="no");

      if(BranchMode[B_DMCSTAGE]) //
        ToDoSteps = iParam[B_ENERGYUPDATEINTERVAL]-1;
      else
        ToDoSteps = iParam[B_WARMUPSTEPS];

      if(BranchMode[B_POPCONTROL])
      {
        //logN = Feedback*std::log(static_cast<RealType>(iParam[B_TARGETWALKERS]));
        logN = std::log(static_cast<RealType>(iParam[B_TARGETWALKERS]));
      }
      else
      {
        //may set Eref to a safe value
        //if(EnergyHistory.count()<5) Eref -= vParam[EnergyWindowIndex];
        vParam[B_ETRIAL]=0.0;Feedback=0.0;logN=0.0;
      }

      WalkerController->start();
    }
  }

  void SimpleFixedNodeBranch::finalize() 
  {

    if(!WalkerController)
    {//running VMC
      RealType e, w,sigma2;
      MyEstimator->getEnergyAndWeight(e,w,sigma2);
      
      vParam[B_ETRIAL]=vParam[B_EREF]=e/w;
      vParam[B_SIGMA]=std::sqrt(sigma2);

      //this is just to avoid diving by n-1  == 0
      EnergyHist(vParam[B_EREF]);
      EnergyHist(vParam[B_EREF]);

      //add Eref to the DMCEnergyHistory
      DMCEnergyHist(vParam[B_EREF]);
    }

    //write to a file
    write(RootName,true);
  }

  /**  Parse the xml file for parameters
   *@param cur current xmlNode 
   *@param LogOut ostream to which the run-time report is sent
   *
   * Few important parameters are:
   * <ul>
   * <li> en_ref: a reference energy
   * <li> num_gen: number of generations \f$N_G\f$ to reach  equilibrium, used in the feedback parameter
   * \f$ feed = \frac{1}{N_G \tau} \f$ 
   * </ul>
   */
  bool SimpleFixedNodeBranch::put(xmlNodePtr cur)
  {
    //save it
    myNode=cur;
    //check dmc/vmc and decide to create WalkerControllerBase
    m_param.put(cur);
    reset();
    MyEstimator->setCollectionMode(true); //always collect
    return true;
  }

  void SimpleFixedNodeBranch::write(const string& fname, bool overwrite) 
  {
    RootName=fname;
    if(MyEstimator->is_manager())
    {
      //\since 2008-06-24
      vParam[B_ACC_ENERGY]=EnergyHist.result();
      vParam[B_ACC_SAMPLES]=EnergyHist.count();
      BranchIO hh(*this,MyEstimator->getCommunicator());
      bool success= hh.write(fname);
    }
  }

  void SimpleFixedNodeBranch::read(const string& fname) 
  {
    RootName=fname;
    if(RootName.find(hdf::config_ext)>=RootName.size())
    {
      RootName.append(hdf::config_ext);
    }

    vParam[B_ACC_ENERGY]=EnergyHist.result();
    vParam[B_ACC_SAMPLES]=EnergyHist.count();

    BranchIO hh(*this,MyEstimator->getCommunicator());
    bool success=hh.read(RootName);

    if(success)
    {
      EnergyHist.clear();
      if(BranchMode[B_DMC] && BranchMode[B_DMCSTAGE])
      {
        app_log() << "  Restarting main DMC" << endl;
        EnergyHist(vParam[B_ACC_ENERGY]/vParam[B_ACC_SAMPLES],vParam[B_ACC_SAMPLES]);
      }
      else
      {
        app_log() << "  Restarting VMC or warm-up DMC" << endl;
        EnergyHist(vParam[B_ACC_ENERGY]/vParam[B_ACC_SAMPLES]);
      }
      app_log() << "    Cummulative mean energy = " << EnergyHist.mean() << endl;
      app_log() << "    Cummaltive samples = " << EnergyHist.count() << endl;
    }
  }
}

/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

