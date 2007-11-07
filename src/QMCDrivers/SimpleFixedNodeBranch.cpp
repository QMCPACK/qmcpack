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
#include "Estimators/DMCEnergyEstimator.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/HDFSTLAttrib.h"

//#include <boost/archive/text_oarchive.hpp>

namespace qmcplusplus {

  SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
    FixedNumWalkers(false), QMCCounter(-1), SwapMode(0), Counter(0), 
  Nideal(nideal), NumGeneration(-1), BranchInterval(1),
  Tau(tau), Feedback(1.0), 
  Eref(0.0),Etrial(0.0), DeltaE(1.0),EavgSum(0.0), WgtSum(0.0), PopControl(0.1), 
  WalkerController(0), MyEstimator(0), SwapWalkers("yes")
  {
    registerParameters();
    reset();
  }

  /** copy constructor
   *
   * Copy only selected data members and WalkerController is never copied.
   */
  SimpleFixedNodeBranch::SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch):
    Counter(0), WalkerController(0),MyEstimator(0)
    {
      //copy the values
      FixedNumWalkers=abranch.FixedNumWalkers;
      QMCCounter=abranch.QMCCounter;
      SwapMode=abranch.SwapMode;
      Nideal=abranch.Nideal;
      NumGeneration=abranch.NumGeneration;
      BranchInterval=abranch.BranchInterval;
      Tau=abranch.Tau;
      Feedback=abranch.Feedback;
      Eref=abranch.Eref;
      Etrial=abranch.Etrial;
      DeltaE=abranch.DeltaE;
      EavgSum=abranch.EavgSum;
      WgtSum=abranch.WgtSum;
      SwapMode=abranch.SwapMode;
      branchCutoff=abranch.branchCutoff;
      branchMax=abranch.branchMax;
      branchFilter=abranch.branchFilter;

      registerParameters();
      reset();
    }

  void SimpleFixedNodeBranch::registerParameters() {
    m_param.add(Feedback,"feedback","double"); 
    m_param.add(Etrial,"refEnergy","AU");
    m_param.add(DeltaE,"refVariance","AU");
    m_param.add(NumGeneration,"popControl","int"); 
    m_param.add(Nideal,"targetWalkers","int"); 
    m_param.add(EavgSum,"energySum","AU"); 
    m_param.add(WgtSum,"weightSum","none");
    m_param.add(PopControl,"swapTrigger","none");
    m_param.add(SwapWalkers,"swapWalkers","string");
    m_param.add(SwapWalkers,"collect","string");
    m_param.add(BranchInterval,"branchInterval","int");

    //backward compatability
    m_param.add(Etrial,"ref_energy","AU"); m_param.add(Etrial,"en_ref","AU");
    m_param.add(NumGeneration,"pop_control","int"); 
    m_param.add(NumGeneration,"num_gen","int"); 
    m_param.add(Nideal,"target_walkers","int");
    m_param.add(PopControl,"swap_trigger","none");
    m_param.add(SwapWalkers,"swap_walkers","string");
  }

  void 
    SimpleFixedNodeBranch::start(const string& froot, bool append) 
    {
      RootName=froot;
      MyEstimator->RootName=froot;
      MyEstimator->reset();
    }

  void SimpleFixedNodeBranch::initWalkerController(RealType tau, bool fixW) {
    Tau=tau;

    ////reset Feedback pararmeter
    //reset();

    if(WalkerController == 0) 
    {
      FixedNumWalkers=fixW;
      //moved these to reset
      //if(fixW) 
      //{
      //  if(WgtSum<5) Eref-= DeltaE;
      //  Etrial=0.0;Feedback=0.0;logN=0.0;
      //}

      WalkerController = createWalkerController(Nideal, MyEstimator->getCommunicator(), myNode);
      //WalkerController 
      //  = CreateWalkerController(FixedNumWalkers, SwapMode, 
      //      Nideal, Nmax, Nmin, WalkerController,MyEstimator->getCommunicator());

      Nmax=WalkerController->Nmax;
      Nmin=WalkerController->Nmin;

      WalkerController->start();

      //use DMCEnergyEstimator
      //DMCEnergyEstimator* dmcE=new DMCEnergyEstimator;
      //dmcE->setWalkerControl(WalkerController);
      //MyEstimator->add(dmcE);
    }

    //reset Feedback pararmeter
    this->reset();

    MyEstimator->reset();
    //update the simulation parameters
    WalkerController->put(myNode);
    //assign current Eref and a large number for variance
    WalkerController->setEnergyAndVariance(Eref,DeltaE);

    //determine the branch cutoff to limit wild weights
    branchCutoff=DeltaE*WalkerController->targetSigma;
    branchMax=branchCutoff*1.5;
    branchFilter=1.0/(branchMax-branchCutoff);

    //reset
    WalkerController->reset();

    //if(!FixedNumWalkers)
    //{
    //  EtrialIndex = MyEstimator->addColumn("Etrial");
    //  MyEstimator->addColumn("Popupation");
    //}

    app_log() << "  QMC counter      = " << QMCCounter << endl;
    app_log() << "  reference energy = " << Eref << endl;
    app_log() << "  trial energy     = " << Etrial << endl;
    app_log() << "  reference variance = " << DeltaE << endl;
    app_log() << "  branch cutoff = " << branchCutoff << " " << branchMax << endl;
    app_log() << "  target_walkers = " << Nideal << endl;
    app_log() << "  Max and mimum walkers per node= " << Nmax << " " << Nmin << endl;
    app_log() << "  Feedback = " << Feedback <<  endl;
  }

  void SimpleFixedNodeBranch::flush(int counter) 
  {
    if(counter==0 && WalkerController) 
      WalkerController->reset();
    Counter=counter;
  }

  void
    SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& Walkers) 
    {
      int pop_old= Walkers.getGlobalNumWalkers();
      //collect the total weights and redistribute the walkers accordingly
      int pop_now = WalkerController->branch(iter,Walkers,PopControl);

      //a smarter average may work better but now simply add one
      EavgSum += WalkerController->EnsembleProperty.Energy; 
      WgtSum += 1.0;

      if(FixedNumWalkers)
      {
        WalkerController->setTrialEnergy(EavgSum/WgtSum);
      }
      else
      {//use an average: instantenous energy is good, too
        Etrial=Eref= EavgSum/WgtSum-Feedback*std::log(static_cast<RealType>(pop_now))+logN;
        WalkerController->setTrialEnergy(Etrial);
        //MyEstimator->setColumn(EtrialIndex,Etrial);
        //MyEstimator->setColumn(EtrialIndex+1,wgt_cur);
      }


      //evaluate everything else
      MyEstimator->accumulate(Walkers);
    }

  /** perform branching
   *
   * Set the trial energy of clones
   */
  void 
    SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones) 
    {
      branch(iter,w);
      for(int i=0; i<clones.size(); i++) 
      { 
        clones[i]->Eref=Eref;
        clones[i]->Etrial=Etrial; 
      }
    }

  void SimpleFixedNodeBranch::reset() 
  {
    //use effective time step of BranchInterval*Tau
    //Feed = 1.0/(static_cast<RealType>(NumGeneration*BranchInterval)*Tau);
    //logN = Feed*std::log(static_cast<RealType>(Nideal));
    if(WalkerController)
    {
      if(FixedNumWalkers)
      {
        if(WgtSum<5) Eref-= DeltaE;
        Etrial=0.0;Feedback=0.0;logN=0.0;
      }
      else
        logN = Feedback*std::log(static_cast<RealType>(Nideal));

      WalkerController->start();
    }
  }

  void SimpleFixedNodeBranch::finalize() {

    //Estimator now is handled by Mover classs
    //MyEstimator->stop();
    //MyEstimator->finalize();
    if(!WalkerController) // || EtrialIndex<0) 
    {//running VMC
      MyEstimator->getEnergyAndWeight(EavgSum,WgtSum);
      EavgSum/=WgtSum;
      WgtSum=1;
      Etrial=Eref=EavgSum;
      //E_T=EavgSum/WgtSum;
      //EavgSum=E_T;
      //WgtSum=1;
    }
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
  bool SimpleFixedNodeBranch::put(xmlNodePtr cur){
    //save it
    myNode=cur;
    //check dmc/vmc and decide to create WalkerControllerBase
    m_param.put(cur);

    //set the SwapMode
    SwapMode = (SwapWalkers=="yes");

    /* \warning{backward compatability qmc/@collect="yes|no"}
    */
    const xmlChar* t=xmlGetProp(cur,(const xmlChar*)"collect");
    if(t != NULL) 
    {
      SwapMode = xmlStrEqual(t,(const xmlChar*)"yes");
    } 

    //overwrite the SwapMode with the number of contexts
    SwapMode = (SwapMode && (OHMMS::Controller->ncontexts()>1));

    if(SwapMode) {
      app_log() << "  Collective handling of Average/walkers" << endl;
    } else {
      app_log() << "  Node-parallel handling of Average/walkers" << endl;
    }

    reset();
    //flush(0);
    //MyEstimator->put(myNode);
    MyEstimator->setCollectionMode(SwapMode);
    return true;
  }

#if defined(HAVE_LIBHDF5)
  void SimpleFixedNodeBranch::write(hid_t grp, bool append) {
    TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
    HDFAttribIO<TinyVector<RealType,3> > eh(esave,append);
    eh.write(grp,"Summary");
    if(LogNorm.size())
    {
      HDFAttribIO<vector<RealType> > lh(LogNorm,append);
      lh.write(grp,"LogNorm");
    }
  }

  void SimpleFixedNodeBranch::read(hid_t grp) {
    herr_t status = H5Eset_auto(NULL, NULL);
    status = H5Gget_objinfo (grp, "Summary", 0, NULL);

    if(status == 0) {
      TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
      HDFAttribIO<TinyVector<RealType,3> > eh(esave);
      eh.read(grp,"Summary");
      Etrial=Eref=esave[0]; EavgSum=esave[1]; WgtSum=esave[2];
      bool normFound(0);
      if(LogNorm.size()) {
        HDFAttribIO<vector<RealType> > lh(LogNorm);
        lh.read(grp,"LogNorm");
        normFound=1;
      }
      if(normFound)app_log() << " Normalization factor defined from previous calculation " << endl;
    } else {
      app_log() << "  Summary is not found. Starting from scratch" << endl;
    }
  }

  void SimpleFixedNodeBranch::read(const string& fname) {
    string h5file = fname;
    string ext=getExtension(h5file);
    if(ext != "h5") { //if the filename does not h5 extension, add the extension
      h5file.append(".config.h5");
    }
    hid_t h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t h_config = H5Gopen(h_file,"config_collection");
    read(h_config);
    H5Gclose(h_config);
    H5Fclose(h_file);
  }
#else
  void SimpleFixedNodeBranch::write(hid_t grp, bool append) { }
  void SimpleFixedNodeBranch::read(hid_t grp) {}
  void SimpleFixedNodeBranch::read(const string& fname) {}
#endif
}

/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

