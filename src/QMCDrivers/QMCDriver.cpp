/////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#include "QMCDrivers/QMCDriver.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "HDFVersion.h"
#include <qmc_common.h>
#include <limits>

namespace qmcplusplus
{

QMCDriver::QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool)
  : MPIObjectBase(0), branchEngine(0)
  , W(w), Psi(psi), H(h), psiPool(ppool), Estimators(0), qmcNode(NULL), wOut(0)
{
  //set defaults
  ResetRandom=false;
  AppendRun=false;
  DumpConfig=false;
  ConstPopulation=true; //default is a fixed population method
  MyCounter=0;
  //<parameter name=" "> value </parameter>
  //accept multiple names for the same value
  //recommend using all lower cases for a new parameter
  RollBackBlocks=0;
  m_param.add(RollBackBlocks,"rewind","int");
  Period4CheckPoint=-1;
  storeConfigs=0;
  m_param.add(storeConfigs,"storeConfigs","int");
  m_param.add( storeConfigs,"storeconfigs","int");
  m_param.add( storeConfigs,"store_configs","int");
  Period4CheckProperties=100;
  m_param.add(Period4CheckProperties,"checkProperties","int");
  m_param.add(Period4CheckProperties,"checkproperties","int");
  m_param.add(Period4CheckProperties,"check_properties","int");
  Period4WalkerDump=0;
  m_param.add(Period4WalkerDump,"recordWalkers","int");
  m_param.add(Period4WalkerDump,"record_walkers","int");
  m_param.add(Period4WalkerDump,"recordwalkers","int");
  Period4ConfigDump=0;
  m_param.add(Period4ConfigDump,"recordConfigs","int");
  m_param.add(Period4ConfigDump,"recordconfigs","int");
  m_param.add(Period4ConfigDump,"record_configs","int");
  CurrentStep=0;
  m_param.add(CurrentStep,"current","int");
  nBlocks=1;
  m_param.add(nBlocks,"blocks","int");
  nSteps=10;
  m_param.add(nSteps,"steps","int");
  nSubSteps=1;
  m_param.add(nSubSteps,"substeps","int");
  m_param.add(nSubSteps,"subSteps","int");
  m_param.add(nSubSteps,"sub_steps","int");
  nWarmupSteps=0;
  m_param.add(nWarmupSteps,"warmupsteps","int");
  m_param.add(nWarmupSteps,"warmupSteps","int");
  m_param.add(nWarmupSteps,"warmup_steps","int");
  nAccept=0;
  nReject=0;
  nTargetWalkers=W.getActiveWalkers();
  m_param.add(nTargetWalkers,"walkers","int");
  //sample-related parameters
  //samples will set nTargetPopulation
  nTargetSamples=0;
  nStepsBetweenSamples=1;
  m_param.add(nStepsBetweenSamples,"stepsbetweensamples","int");
  nSamplesPerThread=0;
  m_param.add(nSamplesPerThread,"samplesperthread","real");
  m_param.add(nSamplesPerThread,"dmcwalkersperthread","real");
  nTargetPopulation=0;
  m_param.add(nTargetPopulation,"samples","real");
  Tau=0.1;
  m_param.add(Tau,"timeStep","AU");
  m_param.add(Tau,"timestep","AU");
  m_param.add(Tau,"time_step","AU");
  m_param.add(Tau,"Tau","AU");
  m_param.add(Tau,"tau","AU");
  MaxCPUSecs=360000; //100 hours
  m_param.add(MaxCPUSecs,"maxcpusecs","real");
  QMCType="invalid";
  ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
}

QMCDriver::~QMCDriver()
{
  delete_iter(Rng.begin(),Rng.end());
}

void QMCDriver::add_H_and_Psi(QMCHamiltonian* h, TrialWaveFunction* psi)
{
  H1.push_back(h);
  Psi1.push_back(psi);
}

/** process a <qmc/> element
 * @param cur xmlNode with qmc tag
 *
 * This function is called before QMCDriver::run and following actions are taken:
 * - Initialize basic data to execute run function.
 * -- distance tables
 * -- resize deltaR and drift with the number of particles
 * -- assign cur to qmcNode
 * - process input file
 *   -- putQMCInfo: <parameter/> s for generic QMC
 *   -- put : extra data by derived classes
 * - initialize branchEngine to accumulate energies
 * - initialize Estimators
 * - initialize Walkers
 */
void QMCDriver::process(xmlNodePtr cur)
{
  deltaR.resize(W.getTotalNum());
  drift.resize(W.getTotalNum());
  qmcNode=cur;
  //process common parameters
  putQMCInfo(cur);
  ////set the Tau parameter inside the Hamiltonian
  //H.setTau(Tau);
  //need to initialize properties
  int numCopies= (H1.empty())?1:H1.size();
  W.resetWalkerProperty(numCopies);
  //create branchEngine first
  if(branchEngine==0)
  {
    branchEngine = new BranchEngineType(Tau,W.getGlobalNumWalkers());
  }
  //execute the put function implemented by the derived classes
  put(cur);
  //create and initialize estimator
  Estimators = branchEngine->getEstimatorManager();
  if(Estimators==0)
  {
    Estimators = new EstimatorManager(myComm);
    branchEngine->setEstimatorManager(Estimators);
    branchEngine->read(h5FileRoot);
  }
  branchEngine->put(cur);
  Estimators->put(W,H,cur);
  if(wOut==0)
    wOut = new HDFWalkerOutput(W,RootName,myComm);
  branchEngine->start(RootName);
  branchEngine->write(RootName);
  //use new random seeds
  if(ResetRandom)
  {
    app_log() << "  Regenerate random seeds." << endl;
    RandomNumberControl::make_seeds();
    ResetRandom=false;
  }
  //flush the ostreams
  OhmmsInfo::flush();
  //increment QMCCounter of the branch engine
  branchEngine->advanceQMCCounter();
}

void QMCDriver::setStatus(const string& aname, const string& h5name, bool append)
{
  RootName = aname;
  app_log() << "\n========================================================="
            << "\n  Start " << QMCType
            << "\n  File Root " << RootName;
  if(append)
    app_log() << " append = yes ";
  else
    app_log() << " append = no ";
  app_log() << "\n=========================================================" << endl;
  if(h5name.size())
    h5FileRoot = h5name;
  AppendRun = append;
}


/** Read walker configurations from *.config.h5 files
 * @param wset list of xml elements containing mcwalkerset
 */
void QMCDriver::putWalkers(vector<xmlNodePtr>& wset)
{
  if(wset.empty())
    return;
  int nfile=wset.size();
  HDFWalkerInputManager W_in(W,myComm);
  for(int i=0; i<wset.size(); i++)
    if(W_in.put(wset[i]))
      h5FileRoot = W_in.getFileRoot();
  //clear the walker set
  wset.clear();

  int nwtot=W.getActiveWalkers();
  myComm->bcast(nwtot);

  if(nwtot)
  {
    int np=myComm->size();
    vector<int> nw(np,0), nwoff(np+1,0);
    nw[myComm->rank()]=W.getActiveWalkers();
    myComm->allreduce(nw);
    for(int ip=0; ip<np; ++ip)
      nwoff[ip+1]=nwoff[ip]+nw[ip];
    W.setGlobalNumWalkers(nwoff[np]);
    W.setWalkerOffsets(nwoff);
    qmc_common.is_restart=true;
  }
  else
    qmc_common.is_restart=false;
}

void QMCDriver::recordBlock(int block)
{
  ////first dump the data for restart
  if(DumpConfig &&block%Period4CheckPoint == 0)
  {
    wOut->dump(W);
    branchEngine->write(RootName,true); //save energy_history
    RandomNumberControl::write(RootName,myComm);
//       if (storeConfigs) wOut->dump( ForwardWalkingHistory);
  }
  //save positions for optimization: this is done within VMC
  //if(QMCDriverMode[QMC_OPTIMIZE]) W.saveEnsemble();
  //if(Period4WalkerDump>0) wOut->append(W);
  //flush the ostream
  //OhmmsInfo::flush();
}

bool QMCDriver::finalize(int block, bool dumpwalkers)
{
  TimerManager.print(myComm);
  TimerManager.reset();
  if(DumpConfig && dumpwalkers)
    wOut->dump(W);
  branchEngine->finalize(W);
  RandomNumberControl::write(RootName,myComm);
  delete wOut;
  wOut=0;
  //Estimators->finalize();
  //set the target walkers
  nTargetWalkers = W.getActiveWalkers();
  //increment MyCounter
  MyCounter++;
  //flush the ostream
  OhmmsInfo::flush();
  return true;
}

/** Add walkers to the end of the ensemble of walkers.
 * @param nwalkers number of walkers to add
 */
void
QMCDriver::addWalkers(int nwalkers)
{
  if(nwalkers>0)
  {
    //add nwalkers walkers to the end of the ensemble
    int nold = W.getActiveWalkers();
    app_log() << "  Adding " << nwalkers << " walkers to " << nold << " existing sets" << endl;
    W.createWalkers(nwalkers);
    if(nold)
    {
      int iw=nold;
      for(MCWalkerConfiguration::iterator it=W.begin()+nold; it != W.end(); ++it,++iw)
        (*it)->R=W[iw%nold]->R;//assign existing walker configurations when the number of walkers change
    }
  }
  else
    if(nwalkers<0)
    {
      W.destroyWalkers(-nwalkers);
      app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers() << endl;
    }
    else
    {
      app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." <<  endl;
    }
  //update the global number of walkers
  //int nw=W.getActiveWalkers();
  //myComm->allreduce(nw);
  vector<int> nw(myComm->size(),0),nwoff(myComm->size()+1,0);
  nw[myComm->rank()]=W.getActiveWalkers();
  myComm->allreduce(nw);
  for(int ip=0; ip<myComm->size(); ip++)
    nwoff[ip+1]=nwoff[ip]+nw[ip];
  W.setGlobalNumWalkers(nwoff[myComm->size()]);
  W.setWalkerOffsets(nwoff);
  app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples  <<  endl;
  app_log() << "  Total weight: " << W.EnsembleProperty.Weight  <<  endl;
}


/** Parses the xml input file for parameter definitions for a single qmc simulation.
 *
 * Basic parameters are handled here and each driver will perform its own initialization with the input
 * attribute list
 * - checkpoint="-1|0|n" default=-1
 *   -- 1 = do not write anything
 *   -- 0 = dump after the completion of a qmc section
 *   -- n = dump after n blocks
 */
bool QMCDriver::putQMCInfo(xmlNodePtr cur)
{
  //SpeciesSet tspecies(W.getSpeciesSet());
  //RealType mass = tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]));
  //if (mass < 1e-12) {
  //  mass=1.0;
  //  tspecies(tspecies.addAttribute("mass"),tspecies.addSpecies(tspecies.speciesName[W.GroupID[0]]))=1.0;
  //}
  //oneovermass = 1.0/mass;
  //set the default walker to the number of threads times 10
  int defaultw = omp_get_max_threads();
  OhmmsAttributeSet aAttrib;
  aAttrib.add(Period4CheckPoint,"checkpoint");
  aAttrib.put(cur);
  if(cur != NULL)
  {
    //initialize the parameter set
    m_param.put(cur);
    xmlNodePtr tcur=cur->children;
    //determine how often to print walkers to hdf5 file
    while(tcur != NULL)
    {
      string cname((const char*)(tcur->name));
      if(cname == "record")
      {
        //dump walkers for optimization
        OhmmsAttributeSet rAttrib;
        rAttrib.add(Period4WalkerDump,"stride");
        rAttrib.add(Period4WalkerDump,"period");
        rAttrib.put(tcur);
      }
      else if(cname == "checkpoint")
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(Period4CheckPoint,"stride");
        rAttrib.add(Period4CheckPoint,"period");
        rAttrib.put(tcur);
        //DumpConfig=(Period4CheckPoint>0);
      }
      else if(cname == "dumpconfig")
      {
        OhmmsAttributeSet rAttrib;
        rAttrib.add(Period4ConfigDump,"stride");
        rAttrib.add(Period4ConfigDump,"period");
        rAttrib.put(tcur);
      } else if(cname == "random")
        {
          ResetRandom = true;
        }
      tcur=tcur->next;
    }
  }

  int oldStepsBetweenSamples=nStepsBetweenSamples;
  //set the minimum blocks
  if (nBlocks<1)
    nBlocks=1;
  if(qmc_common.is_restart  || qmc_common.qmc_counter || !ConstPopulation)
  {
    app_log() << "Using the driver from the previous qmc section. Not resetting any variables concerning samples or walkers" << endl;
  }
  else
  {
    //compute samples and overwrite steps for the given samples
    int Nthreads = omp_get_max_threads();
    int Nprocs=myComm->size();
    //nTargetWalkers is a local quantity, always set to multiple of number of threads
    nTargetWalkers=(std::max(Nthreads,nTargetWalkers)/Nthreads)*Nthreads;
    //target samples set by samples or samplesperthread/dmcwalkersperthread
    nTargetPopulation=std::max(nTargetPopulation,nSamplesPerThread*Nprocs*Nthreads);
    nTargetSamples=static_cast<int>(std::ceil(nTargetPopulation));
    if(nBlocks==1) nBlocks=nSamplesPerThread;
    if(nTargetSamples)
    {
      int nwtot=nTargetWalkers*Nprocs;  //total number of walkers used by this qmcsection
      nTargetSamples=std::max(nwtot,nTargetSamples);
      nTargetSamples=((nTargetSamples+nwtot-1)/nwtot)*nwtot; // nTargetSamples are always multiples of total number of walkers
      nSamplesPerThread=nTargetSamples/Nprocs/Nthreads;
      int ns_target=nTargetSamples*nStepsBetweenSamples; //total samples to generate
      int ns_per_step=Nprocs*nTargetWalkers;  //total samples per step
      nSteps=std::max(nSteps,(ns_target/ns_per_step+nBlocks-1)/nBlocks);
      Period4WalkerDump=nStepsBetweenSamples=ns_per_step*nSteps*nBlocks/nTargetSamples;
    }
    else
    {
      Period4WalkerDump = nStepsBetweenSamples=(nBlocks+1)*nSteps; //some positive number, not used
      nSamplesPerThread=0;
      //nTargetPopulation = nTargetWalkers*Nprocs;
    }
  }
  DumpConfig=(Period4CheckPoint>=0);
  if(Period4CheckPoint<1)
    Period4CheckPoint=nBlocks;
  //reset CurrentStep to zero if qmc/@continue='no'
  if(!AppendRun)
    CurrentStep=0;
  //if walkers are initialized via <mcwalkerset/>, use the existing one
  if(!(qmc_common.qmc_counter || qmc_common.is_restart))
  {
    //always reset the walkers
    int nw  = W.getActiveWalkers();
    int ndiff = 0;
    if(nw)
    {
      // nTargetWalkers == 0, if it is not set by the input file
      ndiff = (nTargetWalkers)? nTargetWalkers-nw: 0;
    }
    else
    {
      ndiff= (nTargetWalkers)? nTargetWalkers:defaultw;
    }
    addWalkers(ndiff);
  }
  app_log() << "\n QMCDriver::putQMCInfo " << endl;
  app_log() << "  timestep       = " << Tau << endl;
  app_log() << "  blocks         = " << nBlocks << endl;
  app_log() << "  steps          = " << nSteps << endl;
  app_log() << "  substeps       = " << nSubSteps << endl;
  app_log() << "  current        = " << CurrentStep << endl;
  app_log() << "  target samples = " << nTargetPopulation << endl;
  app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << endl << endl;
  if(nStepsBetweenSamples != oldStepsBetweenSamples)
    app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << "  (input="<<oldStepsBetweenSamples<<")" << endl;
  else
    app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << endl;
  
  if(DumpConfig)
    app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint << " blocks" << endl;
  else
    app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << endl;
  if (Period4WalkerDump>0)
    app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << endl;
  app_log().flush();
  //always true
  return (W.getActiveWalkers()>0);
}

xmlNodePtr QMCDriver::getQMCNode()
{
  xmlNodePtr newqmc = xmlCopyNode(qmcNode,1);
  xmlNodePtr current_ptr=NULL;
  xmlNodePtr cur=newqmc->children;
  while(cur != NULL && current_ptr == NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "parameter")
    {
      const xmlChar* aptr= xmlGetProp(cur, (const xmlChar *) "name");
      if(aptr)
      {
        if(xmlStrEqual(aptr,(const xmlChar*)"current"))
          current_ptr=cur;
      }
    }
    cur=cur->next;
  }
  if(current_ptr == NULL)
  {
    current_ptr = xmlNewTextChild(newqmc,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0");
    xmlNewProp(current_ptr,(const xmlChar*)"name",(const xmlChar*)"current");
  }
  getContent(CurrentStep,current_ptr);
  return newqmc;
}

}


/***************************************************************************
 * $RCSfile: QMCDriver.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
