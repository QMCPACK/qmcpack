//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCDriver.h"
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
#include <typeinfo>

#include <ADIOS/ADIOS_config.cpp>
#ifdef HAVE_ADIOS
#include <adios.h>
#include "Particle/AdiosWalkerInput.h"
#include <ADIOS/ADIOS_profile.h>
#endif
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif


namespace qmcplusplus
{

QMCDriver::QMCDriver(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, WaveFunctionPool& ppool)
  : MPIObjectBase(0), branchEngine(0), W(w), Psi(psi), H(h), psiPool(ppool),
    Estimators(0),Traces(0), qmcNode(NULL), wOut(0)
{
  ResetRandom=false;
  AppendRun=false;
  DumpConfig=false;
  ConstPopulation=true; //default is a fixed population method
  IsQMCDriver=true;
  allow_traces = false;
  MyCounter=0;
  //<parameter name=" "> value </parameter>
  //accept multiple names for the same value
  //recommend using all lower cases for a new parameter
  RollBackBlocks=0;
  m_param.add(RollBackBlocks,"rewind","int");
  Period4CheckPoint=-1;
  storeConfigs=0;
  //m_param.add(storeConfigs,"storeConfigs","int");
  m_param.add( storeConfigs,"storeconfigs","int");
  m_param.add( storeConfigs,"store_configs","int");
  Period4CheckProperties=100;
  m_param.add(Period4CheckProperties,"checkProperties","int");
  m_param.add(Period4CheckProperties,"checkproperties","int");
  m_param.add(Period4CheckProperties,"check_properties","int");
  Period4WalkerDump=0;
  //m_param.add(Period4WalkerDump,"recordWalkers","int");
  m_param.add(Period4WalkerDump,"record_walkers","int");
  m_param.add(Period4WalkerDump,"recordwalkers","int");
  Period4ConfigDump=0;
  //m_param.add(Period4ConfigDump,"recordConfigs","int");
  m_param.add(Period4ConfigDump,"recordconfigs","int");
  m_param.add(Period4ConfigDump,"record_configs","int");
  CurrentStep=0;
  m_param.add(CurrentStep,"current","int");
  nBlocks=1;
  m_param.add(nBlocks,"blocks","int");
  nSteps=1;
  m_param.add(nSteps,"steps","int");
  nSubSteps=1;
  m_param.add(nSubSteps,"substeps","int");
  //m_param.add(nSubSteps,"subSteps","int");
  m_param.add(nSubSteps,"sub_steps","int");
  nWarmupSteps=0;
  m_param.add(nWarmupSteps,"warmupsteps","int");
  //m_param.add(nWarmupSteps,"warmupSteps","int");
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
  //m_param.add(Tau,"timeStep","AU");
  m_param.add(Tau,"timestep","AU");
  m_param.add(Tau,"time_step","AU");
  //m_param.add(Tau,"Tau","AU");
  m_param.add(Tau,"tau","AU");
  MaxCPUSecs=360000; //100 hours
  m_param.add(MaxCPUSecs,"maxcpusecs","real");
  // by default call recompute at the end of each block in the mixed precision case.
#ifdef QMC_CUDA
  if (typeid(CudaRealType) == typeid(float))
  {
    // gpu mixed precision
    nBlocksBetweenRecompute = 1;
  }
  else if (typeid(CudaRealType) == typeid(double))
  {
    // gpu double precision
    nBlocksBetweenRecompute = 0;
  }
#else
#ifdef MIXED_PRECISION
  // cpu mixed precision
  nBlocksBetweenRecompute = 1;
#else
  // cpu double precision
  nBlocksBetweenRecompute = 0;
#endif
#endif
  m_param.add(nBlocksBetweenRecompute,"blocks_between_recompute","int");
  QMCType="invalid";
  ////add each QMCHamiltonianBase to W.PropertyList so that averages can be taken
  //H.add2WalkerProperty(W);
  //if (storeConfigs) ForwardWalkingHistory.storeConfigsForForwardWalking(w);
  rotation = 0;

  checkpointTimer = TimerManager.createTimer("checkpoint::recordBlock", timer_level_medium);
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
    Estimators = new EstimatorManagerBase(myComm);
    branchEngine->setEstimatorManager(Estimators);
    branchEngine->read(h5FileRoot);
  }
#if !defined(REMOVE_TRACEMANAGER)
  //create and initialize traces
  if(Traces==0)
  {
    Traces = new TraceManager(myComm);
  }
  Traces->put(traces_xml,allow_traces,RootName);
#endif
  branchEngine->put(cur);
  Estimators->put(W,H,cur);
  if(wOut==0)
    wOut = new HDFWalkerOutput(W,RootName,myComm);
  branchEngine->start(RootName);
  branchEngine->write(RootName);
  //use new random seeds
  if(ResetRandom)
  {
    app_log() << "  Regenerate random seeds." << std::endl;
    RandomNumberControl::make_seeds();
    ResetRandom=false;
  }
  //flush the std::ostreams
  infoSummary.flush();
  infoLog.flush();
  //increment QMCCounter of the branch engine
  branchEngine->advanceQMCCounter();
}

void QMCDriver::setStatus(const std::string& aname, const std::string& h5name, bool append)
{
  RootName = aname;
  app_log() << "\n========================================================="
            << "\n  Start " << QMCType
            << "\n  File Root " << RootName;
  if(append)
    app_log() << " append = yes ";
  else
    app_log() << " append = no ";
  app_log() << "\n=========================================================" << std::endl;
  if(h5name.size())
    h5FileRoot = h5name;
  AppendRun = append;
}


/** Read walker configurations from *.config.h5 files
 * @param wset list of xml elements containing mcwalkerset
 */
void QMCDriver::putWalkers(std::vector<xmlNodePtr>& wset)
{
  if(wset.empty()) return;
  int nfile=wset.size();
  HDFWalkerInputManager W_in(W,myComm);
  for(int i=0; i<wset.size(); i++)
    if(W_in.put(wset[i])) h5FileRoot = W_in.getFileRoot();
  //clear the walker set
  wset.clear();
  int nwtot=W.getActiveWalkers();
  myComm->bcast(nwtot);
  if(nwtot)
  {
    int np=myComm->size();
    std::vector<int> nw(np,0), nwoff(np+1,0);
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

std::string QMCDriver::getRotationName( std::string RootName)
{
  std::string r_RootName;
  if(rotation % 2 == 0)
  {
    r_RootName = RootName;
  }
  else
  {
    r_RootName = RootName + ".bk";
  }
  rotation++;
  return r_RootName;
}

std::string QMCDriver::getLastRotationName( std::string RootName)
{
  std::string r_RootName;
  if((rotation -1)%2 == 0)
  {
    r_RootName = RootName;
  }
  else
  {
    r_RootName = RootName + ".bk";
  }
  return r_RootName;
}

void QMCDriver::adiosCheckpoint(int block)
{
#ifdef HAVE_ADIOS
  int64_t adios_handle;
  uint64_t adios_groupsize, adios_totalsize;
  EstimatorManagerBase* myEstimator = branchEngine->getEstimatorManager();
  if (sizeof(OHMMS_PRECISION) == sizeof(double))
  {
    adios_open(&adios_handle, "checkpoint_double", (getRotationName(RootName) + ".config.bp").c_str(), "w", myComm->getMPI());
  }
  else
  {
    adios_open(&adios_handle, "checkpoint_float", (getRotationName(RootName) + ".config.bp").c_str(), "w", myComm->getMPI());
  }

  adios_groupsize = wOut->get_group_size(W);
  adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
  wOut->adios_checkpoint(W, adios_handle, block);

  //if (myEstimator->is_manager())
  //{
  //  BranchIO hh(*branchEngine,myEstimator->getCommunicator());
  //  adios_groupsize = hh.get_Checkpoint_size();
  //  adios_groupsize += RandomNumberControl::get_group_size();
  //  adios_groupsize += wOut->get_group_size(W);
  //  adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
  //  hh.adios_checkpoint(adios_handle);
  //  branchEngine->save_energy();
  //  wOut->adios_checkpoint(W, adios_handle, block);
  //  RandomNumberControl::adios_checkpoint(adios_handle);
  //}
  //else
  //{
  //  adios_groupsize = RandomNumberControl::get_group_size();
  //  adios_groupsize += wOut->get_group_size(W);
  //  adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
  //  wOut->adios_checkpoint(W, adios_handle, block);
  //  RandomNumberControl::adios_checkpoint(adios_handle);
  //}
  adios_close(adios_handle);

#ifdef IO_PROFILE
  ADIOS_PROFILE::profile_adios_size(myComm, ADIOS_PROFILE::CKPOINT, adios_groupsize, adios_totalsize);
#endif
#ifdef ADIOS_VERIFY
  ADIOS_FILE *fp = adios_read_open_file((getLastRotationName(RootName) + ".config.bp").c_str(),
      ADIOS_READ_METHOD_BP,
      myComm->getMPI());
  if (fp == NULL)
    app_error() << "Fail to open adios file "<<(getLastRotationName(RootName) + ".config.bp").c_str()<<" Abort. "<< std::endl;
  //if (myEstimator->is_manager())
  //{
  //  BranchIO hh(*branchEngine,myEstimator->getCommunicator());
  //  hh.adios_checkpoint_verify(fp);
  //}
  //RandomNumberControl::adios_checkpoint_verify(fp);
  wOut->adios_checkpoint_verify(W, fp);
  adios_read_close(fp);
#endif
#endif
}

void QMCDriver::adiosCheckpointFinal(int block, bool dumpwalkers)
{
#ifdef HAVE_ADIOS
  int64_t adios_handle;
  //string group_name;
  //get the size of walker related data that we are writing to disk
  //adios_groupsize += RandomNumberControl::get_group_size();
  //EstimatorManager* myEstimator = branchEngine->getEstimatorManager();
  if(DumpConfig && dumpwalkers){
    if (sizeof(OHMMS_PRECISION) == sizeof(double))
    {
      adios_open(&adios_handle, "checkpoint_double", (getRotationName(RootName) + ".config.bp").c_str(), "w", myComm->getMPI());
    }
    else
    {
      adios_open(&adios_handle, "checkpoint_float", (getRotationName(RootName) + ".config.bp").c_str(), "w", myComm->getMPI());
    }
    //if (myEstimator->is_manager())
    //{
    //  //Since we are in the main process we need to write out some more information
    //  BranchIO hh(*branchEngine,myEstimator->getCommunicator());
    //  //Get the size of the data we are writing out for qmc_status
    //  adios_groupsize += hh.get_Checkpoint_size();
    //  //Tell adios how much space we are using for this write out.
    //  adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
    //  //Checkpoint qmc status related data
    //  branchEngine->save_energy(); //save energy_history
    //  hh.adios_checkpoint(adios_handle);
    //}
    //else
    //  adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
    //Checkpoint the data for RandomNumber Control
    
    //adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
    //RandomNumberControl::adios_checkpoint(adios_handle);
    wOut->adios_checkpoint(W, adios_handle, block);
    adios_close(adios_handle);
  }
#ifdef ADIOS_VERIFY
  if(DumpConfig && dumpwalkers){
    ADIOS_FILE *fp = adios_read_open_file((getLastRotationName(RootName) + ".config.bp").c_str(),
                                           ADIOS_READ_METHOD_BP,
                                           myComm->getMPI());
    //if (myEstimator->is_manager())
    //{
    //  BranchIO hh(*branchEngine,myEstimator->getCommunicator());
    //  hh.adios_checkpoint_verify(fp);
    //}
    //RandomNumberControl::adios_checkpoint_verify(fp);
    wOut->adios_checkpoint_verify(W, fp);
    adios_read_close(fp);
  }
#endif
  //if(!ADIOS::useHDF5()){
  //  branchEngine->finalize(W);
  //  delete wOut;
  //  wOut=0;
  //  nTargetWalkers = W.getActiveWalkers();
  //  MyCounter++;
  //  OhmmsInfo::flush();
  //}
#endif
}

void QMCDriver::recordBlock(int block)
{
  if(DumpConfig && block % Period4CheckPoint == 0)
  {
    checkpointTimer->start();
    if(ADIOS::useADIOS())
    {
      adiosCheckpoint(block);
    }
    if(ADIOS::useHDF5()) 
    {
      wOut->dump(W, block);
    }
    branchEngine->write(RootName,true); //save energy_history
    RandomNumberControl::write(RootName,myComm);
    checkpointTimer->stop();
  }
}

bool QMCDriver::finalize(int block, bool dumpwalkers)
{
  if(ADIOS::useADIOS())
  {
    adiosCheckpointFinal(block, dumpwalkers);
  }

  if(ADIOS::useHDF5())
  {
    if(DumpConfig && dumpwalkers) wOut->dump(W, block);
    delete wOut;
    wOut=0;
    //Estimators->finalize();
    nTargetWalkers = W.getActiveWalkers();
    MyCounter++;
    infoSummary.flush();
    infoLog.flush();
  }

  branchEngine->finalize(W);
  
  if(DumpConfig) RandomNumberControl::write(RootName,myComm);

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
    app_log() << "  Adding " << nwalkers << " walkers to " << nold << " existing sets" << std::endl;
    W.createWalkers(nwalkers);
    if(nold)
    {
      int iw=nold;
      for(MCWalkerConfiguration::iterator it=W.begin()+nold; it != W.end(); ++it,++iw)
        (*it)->R=W[iw%nold]->R;//assign existing walker configurations when the number of walkers change
    }
  }
  else if(nwalkers<0)
  {
    W.destroyWalkers(-nwalkers);
    app_log() << "  Removed " << -nwalkers << " walkers. Current number of walkers =" << W.getActiveWalkers() << std::endl;
  }
  else
  {
    app_log() << "  Using the current " << W.getActiveWalkers() << " walkers." <<  std::endl;
  }
  setWalkerOffsets();
  ////update the global number of walkers
  ////int nw=W.getActiveWalkers();
  ////myComm->allreduce(nw);
}

void QMCDriver::setWalkerOffsets()
{
  std::vector<int> nw(myComm->size(),0),nwoff(myComm->size()+1,0);
  nw[myComm->rank()]=W.getActiveWalkers();
  myComm->allreduce(nw);
  for(int ip=0; ip<myComm->size(); ip++)
    nwoff[ip+1]=nwoff[ip]+nw[ip];
  W.setGlobalNumWalkers(nwoff[myComm->size()]);
  W.setWalkerOffsets(nwoff);
  long id=nwoff[myComm->rank()];
  for(int iw=0; iw<nw[myComm->rank()]; ++iw,++id)
  {
    W[iw]->ID       = id;
    W[iw]->ParentID = id;
  }
  app_log() << "  Total number of walkers: " << W.EnsembleProperty.NumSamples  <<  std::endl;
  app_log() << "  Total weight: " << W.EnsembleProperty.Weight  <<  std::endl;
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
  if(!IsQMCDriver)
  {
    app_log() << getName() << "  Skip QMCDriver::putQMCInfo " << std::endl;
    return true;
  }

  ////store the current nSteps and nStepsBetweenSamples 
  //int oldStepsBetweenSamples=nStepsBetweenSamples;
  //int oldSteps=nSteps;

  //set the default walker to the number of threads times 10
  Period4CheckPoint=-1;
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
      std::string cname((const char*)(tcur->name));
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
      }
      else if(cname == "random")
      {
        ResetRandom = true;
      }
      tcur=tcur->next;
    }
  }
  //set the minimum blocks
  if (nBlocks<1) nBlocks=1;

  DumpConfig=(Period4CheckPoint>=0);
  if(Period4CheckPoint<1)
    Period4CheckPoint=nBlocks;
  //reset CurrentStep to zero if qmc/@continue='no'
  if(!AppendRun) CurrentStep=0;

  //if walkers are initialized via <mcwalkerset/>, use the existing one
  if(qmc_common.qmc_counter || qmc_common.is_restart)
  {
    app_log() << "Using existing walkers " << std::endl;
  }
  else
  { //always reset the walkers
    int nths=(qmc_common.compute_device)?1:omp_get_max_threads();
    nTargetWalkers=(std::max(nths,(nTargetWalkers/nths)*nths));
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

  return (W.getActiveWalkers()>0);
}

xmlNodePtr QMCDriver::getQMCNode()
{
  xmlNodePtr newqmc = xmlCopyNode(qmcNode,1);
  xmlNodePtr current_ptr=NULL;
  xmlNodePtr cur=newqmc->children;
  while(cur != NULL && current_ptr == NULL)
  {
    std::string cname((const char*)(cur->name));
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


