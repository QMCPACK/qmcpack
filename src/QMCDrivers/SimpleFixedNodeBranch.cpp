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
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Message/CommUtilities.h"
#include "HDFVersion.h"

//#include <boost/archive/text_oarchive.hpp>

namespace qmcplusplus 
{

  ///enum to yes/no options saved in sParam 
  enum {COMBOPT, USETAUOPT, DUMMYOPT};

  SimpleFixedNodeBranch::SimpleFixedNodeBranch(RealType tau, int nideal): 
    vParam(1.0), WalkerController(0), MyEstimator(0), PopHist(10)
  {

    BranchMode.set(B_DMCSTAGE,0); //warmup stage
    BranchMode.set(B_POPCONTROL,1); //use standard DMC
    BranchMode.set(B_USETAUEFF,1); //use taueff
    BranchMode.set(B_CLEARHISTORY,0); //clear history and start with the current average

    Tau=tau;
    TauEff=tau;
    Feedback=1.0;
    R2Accepted(1.0e-10); 
    R2Proposed(1.0e-10);

    //set the default values for integer parameters
    iParam[B_WARMUPSTEPS]=100;
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
    Tau(abranch.Tau), TauEff(abranch.TauEff), Feedback(abranch.Feedback),
    Eref(abranch.Eref), Etrial(abranch.Etrial),
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
    m_param.add(Etrial,"refEnergy","AU"); 
    m_param.add(Etrial,"ref_energy","AU"); 
    m_param.add(Etrial,"en_ref","AU");

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

    Tau=tau;
    TauEff=tau*R2Accepted.result()/R2Proposed.result();

    if(WalkerController == 0) 
    {
      BranchMode.set(B_POPCONTROL,!fixW);//fixW -> 0 
      WalkerController = createWalkerController(iParam[B_TARGETWALKERS], MyEstimator->getCommunicator(), myNode);
      iParam[B_MAXWALKERS]=WalkerController->Nmax;
      iParam[B_MINWALKERS]=WalkerController->Nmin;
      WalkerController->start();

      PopHist.clear();
      PopHist.reserve(std::max(iParam[B_ENERGYUPDATEINTERVAL],10));
    }

    //reset Feedback pararmeter
    this->reset();

    MyEstimator->reset();
    //update the simulation parameters
    WalkerController->put(myNode);
    //assign current Eref and a large number for variance
    WalkerController->setEnergyAndVariance(Eref,vParam[B_ENERGYWINDOW]);

    //determine the branch cutoff to limit wild weights
    //this can be done with average
    vParam[B_BRANCHCUTOFF]=vParam[B_ENERGYWINDOW]*WalkerController->targetSigma;
    vParam[B_BRANCHMAX]=vParam[B_BRANCHCUTOFF]*1.5;
    vParam[B_BRANCHFILTER]=1.0/(vParam[B_BRANCHMAX]-vParam[B_BRANCHCUTOFF]);

    //reset controller 
    WalkerController->reset();

    app_log() << "  QMC counter      = " << iParam[B_COUNTER] << endl;
    app_log() << "  time step        = " << Tau << endl;
    app_log() << "  effective time step = " << TauEff << endl;
    app_log() << "  Feedback = " << Feedback <<  endl;
    app_log() << "  reference energy = " << Eref << endl;
    app_log() << "  trial energy     = " << Etrial << endl;
    app_log() << "  reference variance = " << vParam[B_ENERGYWINDOW] << endl;
    app_log() << "  branch cutoff = " <<  vParam[B_BRANCHCUTOFF] << " " << vParam[B_BRANCHMAX] << endl;
    app_log() << "  target walkers = " << iParam[B_TARGETWALKERS] << endl;
    app_log() << "  Max and mimum walkers per node= " << iParam[B_MAXWALKERS] << " " << iParam[B_MINWALKERS] << endl;

  }

  void SimpleFixedNodeBranch::flush(int counter) 
  {
    if(counter==0 && WalkerController) WalkerController->reset();
  }

  void SimpleFixedNodeBranch::branch(int iter, MCWalkerConfiguration& Walkers) 
  {
    //collect the total weights and redistribute the walkers accordingly, using a fixed tolerance
    int pop_now = WalkerController->branch(iter,Walkers,0.1);

    EnergyHist(WalkerController->EnsembleProperty.Energy);
    R2Accepted(WalkerController->EnsembleProperty.R2Accepted);
    R2Proposed(WalkerController->EnsembleProperty.R2Proposed);
    PopHist(pop_now);

    Eref=EnergyHist.mean();//current mean

    //EavgSum += ene_now;
    //WgtSum += 1.0;
    if(BranchMode[B_DMCSTAGE]) // main stage
    { 
      if(BranchMode[B_POPCONTROL])
      {
        if(ToDoSteps>0)
          --ToDoSteps;
        else
        {
          //Etrial=Eref-Feedback*std::log(static_cast<RealType>(pop_now))+logN; 
          Etrial=Eref-Feedback*std::log(static_cast<RealType>(PopHist.mean()))+logN; 
          ToDoSteps=iParam[B_ENERGYUPDATEINTERVAL]-1;
          //app_log() << "  UPDATE " << iter << " "  << Etrial << " " << Eref << " " << PopHist.mean() << endl;
        }
      }
    }
    else//warmup
    {
      if(BranchMode[B_USETAUEFF]) TauEff=Tau*R2Accepted.result()/R2Proposed.result();

      if(BranchMode[B_POPCONTROL])
        Etrial=Eref-Feedback*std::log(static_cast<RealType>(pop_now))+logN;

      --ToDoSteps;
      if(ToDoSteps==0)  //warmup is done
      {
        app_log() << "\n Warmup is completed. ";
        if(BranchMode[B_USETAUEFF])
          app_log() << "\n TauEff     = " << TauEff << "\n TauEff/Tau = " << TauEff/Tau ;
        else
          app_log() << "\n TauEff proposed   = " << Tau*R2Accepted.result()/R2Proposed.result();

        app_log()  << "\n Etrial     = " << Etrial << endl;

        ToDoSteps = iParam[B_ENERGYUPDATEINTERVAL]-1;
        BranchMode.set(B_DMCSTAGE,1); //set BranchModex to main stage
      }
    }

    WalkerController->setTrialEnergy(Etrial);

    //app_log() << "  HIST " << EnergyHist.average() << " " << EavgSum/WgtSum << " "  << PopHist.average() << endl;
    //if(FixedNumWalkers)
    //{
    //  WalkerController->setTrialEnergy(EavgSum/WgtSum);
    //}
    //else
    //{//use an average: instantenous energy is good, too
    //  Etrial=Eref= Emix-Feedback*std::log(static_cast<RealType>(pop_now))+logN;
    //  WalkerController->setTrialEnergy(Etrial);
    //  //MyEstimator->setColumn(EtrialIndex,Etrial);
    //  //MyEstimator->setColumn(EtrialIndex+1,wgt_cur);
    //}

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

    for(int i=0; i<clones.size(); i++) 
    { 
      clones[i]->Eref=Eref;
      clones[i]->Etrial=Etrial; 
      clones[i]->TauEff=TauEff; 
    }

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
        logN = Feedback*std::log(static_cast<RealType>(iParam[B_TARGETWALKERS]));
      }
      else
      {
        //may set Eref to a safe value
        //if(EnergyHistory.count()<5) Eref -= vParam[EnergyWindowIndex];
        Etrial=0.0;Feedback=0.0;logN=0.0;
      }

      WalkerController->start();
    }
  }

  void SimpleFixedNodeBranch::finalize() 
  {

    if(!WalkerController)
    {//running VMC
      RealType e, w;
      MyEstimator->getEnergyAndWeight(e,w);
      Etrial=Eref=e/w;

      //this is just to avoid diving by n-1  == 0
      EnergyHist(Eref);
      EnergyHist(Eref);
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

#if defined(HAVE_LIBHDF5)
  void SimpleFixedNodeBranch::write(const string& fname, bool overwrite) 
  {
    RootName=fname;
    if(MyEstimator->is_manager())
    {
      hid_t fid =  H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      hid_t h1 =  H5Gopen(fid,hdf::main_state);
      herr_t status = H5Eset_auto(NULL, NULL);
      status = H5Gget_objinfo(h1,hdf::energy_history,0,NULL);
      //TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
      TinyVector<RealType,3> esave(Eref,EnergyHist.result(),EnergyHist.count());
      overwrite=(status == 0);
      HDFAttribIO<TinyVector<RealType,3> > eh(esave,overwrite);
      eh.write(h1,hdf::energy_history);
      if(LogNorm.size())//check if collection is done correctly
      {
        HDFAttribIO<vector<RealType> > lh(LogNorm,overwrite);
        lh.write(h1,hdf::norm_history);
      }
      H5Gclose(h1);
      H5Fclose(fid);
    }
  }

  //  void SimpleFixedNodeBranch::write(hid_t grp, bool overwrite) {
  //    TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
  //    HDFAttribIO<TinyVector<RealType,3> > eh(esave,overwrite);
  //    eh.write(grp,hdf::energy_history);
  //    if(LogNorm.size())
  //    {
  //      HDFAttribIO<vector<RealType> > lh(LogNorm,overwrite);
  //      lh.write(grp,hdf::norm_history);
  //    }
  //  }


  void SimpleFixedNodeBranch::read(const string& fname) {

    //esave is used for communication
    //TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
    TinyVector<RealType,3> esave(Eref,EnergyHist.result(),EnergyHist.count());

    RootName=fname;
    if(RootName.find(hdf::config_ext)>=RootName.size())
    {
      RootName.append(hdf::config_ext);
    }

    //string ext=getExtension(fname);
    //if(ext != "h5") { //if the filename does not h5 extension, add the extension
    //  RootName.append(".config.h5");
    //}
    if(MyEstimator->is_manager())
    {
      hid_t h_file =  H5Fopen(RootName.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
      if(h_file<0)
      {
        app_error() << "  Failed to open " << RootName << endl;
        return;
      }

      //start using major=0 and minor=4
      HDFVersion res_version(0,4);
      HDFVersion in_version(0,1);
      herr_t status = H5Eset_auto(NULL, NULL);
      in_version.read(h_file,hdf::version);
      //status = H5Gget_objinfo (h_file, hdf::version, 0, NULL);
      //if(status == 0)//version exists after (0,4)
      if(in_version>=res_version)
      {
        //in_version.read(h_file,hdf::version);
        hid_t h1=H5Gopen(h_file,hdf::main_state);
        //get the history
        HDFAttribIO<TinyVector<RealType,3> > eh(esave);
        eh.read(h1,hdf::energy_history);
        if(LogNorm.size()) 
        {
          HDFAttribIO<vector<RealType> > lh(LogNorm);
          lh.read(h1,hdf::norm_history);
        }
        H5Gclose(h1);
      }
      else 
      { 
        app_log() << "  Missing version. Using old format " << endl;
        hid_t h1 = H5Gopen(h_file,"config_collection");
        herr_t status = H5Eset_auto(NULL, NULL);
        status = H5Gget_objinfo (h1, "Summary", 0, NULL);
        if(status == 0) {
          HDFAttribIO<TinyVector<RealType,3> > eh(esave);
          eh.read(h1,"Summary");
          if(LogNorm.size()) 
          {
            HDFAttribIO<vector<RealType> > lh(LogNorm);
            lh.read(h1,"LogNorm");
            app_log() << " Normalization factor defined from previous calculation " << endl;
          }
        } else {
          app_log() << "  Summary is not found. Starting from scratch" << endl;
        }
        H5Gclose(h1);
      }
      H5Fclose(h_file);
    }

    //broadcast to the nodes : need to add a namespace mpi::
    bcast(esave,MyEstimator->getCommunicator());

    Etrial=Eref=esave[0];
    //EavgSum=esave[1];
    //WgtSum=esave[2];
  }

  //  void SimpleFixedNodeBranch::read(hid_t grp) {
  //    //Disable this
  //    app_log() << "  Missing version. Using old format " << endl;
  //    hid_t h_config = H5Gopen(grp,"config_collection");
  //    herr_t status = H5Eset_auto(NULL, NULL);
  //    status = H5Gget_objinfo (h_config, "Summary", 0, NULL);
  //    if(status == 0) {
  //      TinyVector<RealType,3> esave(Eref,EavgSum,WgtSum);
  //      HDFAttribIO<TinyVector<RealType,3> > eh(esave);
  //      eh.read(h_config,"Summary");
  //      Etrial=Eref=esave[0]; EavgSum=esave[1]; WgtSum=esave[2];
  //      bool normFound(0);
  //      if(LogNorm.size()) {
  //        HDFAttribIO<vector<RealType> > lh(LogNorm);
  //        lh.read(h_config,"LogNorm");
  //        normFound=1;
  //      }
  //      if(normFound)app_log() << " Normalization factor defined from previous calculation " << endl;
  //    } else {
  //      app_log() << "  Summary is not found. Starting from scratch" << endl;
  //    }
  //    H5Gclose(h_config);
  //  }

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

