#include<tuple>
#include<map>
#include<string>
#include<iomanip>

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include<Message/MPIObjectBase.h>
#include "Message/OpenMP.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "OhmmsData/libxmldefs.h"
#include "Configuration.h"
#include <qmc_common.h>

#include "AFQMC/config.h"
#include "AFQMC/Drivers/AFQMCDriver.h"

namespace qmcplusplus {

bool AFQMCDriver::run()
{

  if(compare_libraries)
  {
    prop0->benchmark();
    return true;
  }

  if(!restarted) {
    Eshift=wlkBucket->getEloc(0).real();
    Etav=wlkBucket->getEloc(0).real();
    step0=block0=0;
  }

  RealType w0 = wlkBucket->GlobalWeight();
  int nwalk_ini = wlkBucket->GlobalPopulation();
  estim0->setTargetWeight(w0);

  app_log()<<"Initial weight and number of walkers: " <<w0 <<" " <<nwalk_ini <<std::endl;

  std::ofstream out_timers;
  if(print_timers > 0 && myComm->rank()==0) {
    out_timers.open("timers.dat");//,std::ios_base::app | std::ios_base::out); 
  }

  // problems with using step_tot to do ortho and load balance
  int time = step0*nSubstep, step_tot=step0, iBlock ; 
  for (iBlock=block0; iBlock<nBlock; ++iBlock) {

    LocalTimer.start("Block::TOTAL");
    for (int iStep=0; iStep<nStep; ++iStep, ++step_tot) {

      // propagate
      for (int iSubstep=0; iSubstep<nSubstep; ++iSubstep,++time) {

        LocalTimer.start("SubStep::Propagate");
        prop0->Propagate(time,wlkBucket,Eshift,Eshift);
        LocalTimer.stop("SubStep::Propagate");        

        estim0->accumulate_substep(wlkBucket);

      }  // iSubstep

      // quantities that are measured once per step 
      estim0->accumulate_step(wlkBucket);

      if (step_tot != 0 && step_tot % nPopulationControl == 0) {
        LocalTimer.start("Step::PopControl");
        wlkBucket->popControl();
        LocalTimer.stop("Step::PopControl");
      }

      if (step_tot != 0 && step_tot % nloadBalance == 0) {
        LocalTimer.start("Step::loadBalance");
        wlkBucket->loadBalance();
        LocalTimer.stop("Step::loadBalance");
      }    
 
      if (step_tot != 0 && step_tot % nStabalize == 0) { // && it->alive) {
        LocalTimer.start("Step::Orthogonalize");
        wlkBucket->Orthogonalize();
        wfn0->evaluateOverlap("ImportanceSampling",-1,wlkBucket);
        LocalTimer.stop("Step::Orthogonalize");
      }

      //Etav += estim0->getEloc_step();
      if(time*dt < 1.0)  
        Etav = estim0->getEloc_step();
      else
        Etav += dShift*0.1*(estim0->getEloc_step()-Etav);

      if(time*dt < 1.0)  
        Eshift = estim0->getEloc_step();
      else
        Eshift += dShift*(estim0->getEloc_step()-Eshift);

    }

    // checkpoint 
    if(nCheckpoint > 0 && iBlock != 0 && iBlock % nCheckpoint == 0) 
      if(!checkpoint(iBlock,step_tot)) {
        app_error()<<" Error in AFQMCDriver::checkpoint(). \n" <<std::endl;
        return false;
      }

    // write samples
    if(samplePeriod > 0 && iBlock != 0 && iBlock % samplePeriod == 0)
      if(!writeSamples()) {
        app_error()<<" Error in AFQMCDriver::writeSamples(). \n" <<std::endl;
        return false;
      }

    // quantities that are measured once per block
    estim0->accumulate_block(wlkBucket);

    LocalTimer.stop("Block::TOTAL");

    estim0->print(iBlock+1,time*dt,Eshift,Etav,wlkBucket);

    if(print_timers > 0 && myComm->rank()==0 && iBlock%print_timers == 0) output_timers(out_timers,iBlock); 

  }

  checkpoint(iBlock,step_tot);

  app_log()<<"----------------------------------------------------------------\n";
  app_log()<<" LocalTimer: \n";
  LocalTimer.print_average_all(app_log());
  app_log()<<" Timer: \n";
  Timer.print_average_all(app_log());
  app_log()<<"----------------------------------------------------------------\n";

  return true;
}

bool AFQMCDriver::parse(xmlNodePtr cur)
{
  if(cur==NULL) return false;

  std::string str,str1;

  int cmp_lib=0;
  int deb_prop=0;
  ncores_per_TG=1; 
  ParameterSet m_param;
  m_param.add(nBlock,"blocks","int");
  m_param.add(nStep,"steps","int");
  m_param.add(nSubstep,"substeps","int");
  m_param.add(nPopulationControl,"popControl","int");
  m_param.add(nWalkers,"nWalkers","int");
  m_param.add(nStabalize,"ortho","int");
  m_param.add(nloadBalance,"loadBalance","int");
  m_param.add(nCheckpoint,"checkpoint","int");
  m_param.add(samplePeriod,"samplePeriod","int");
  m_param.add(print_timers,"timers","int");
  m_param.add(print_timers,"timer","int");
  m_param.add(print_timers,"print_timers","int");
  m_param.add(print_timers,"print_timer","int");
  m_param.add(ncores_per_TG,"ncores_per_TG","int");
  m_param.add(ncores_per_TG,"ncores","int");
  m_param.add(ncores_per_TG,"cores","int");
  m_param.add(dt,"dt","double");
  m_param.add(dt,"timestep","double");
  m_param.add(dShift,"dshift","double");
  m_param.add(cmp_lib,"test_library","int");
  m_param.add(deb_prop,"debug","int");
  m_param.add(min_total_weight,"min_total_weight","double");
  m_param.add(str1,"set_nWalker_to_target","std::string");
  m_param.add(str1,"set_nwalker_to_target","std::string");
  m_param.add(hdf_read_tag,"hdf_read_tag","std::string");
  m_param.add(hdf_read_restart,"hdf_read_file","std::string");
  m_param.add(hdf_write_tag,"hdf_write_tag","std::string");
  m_param.add(hdf_write_restart,"hdf_write_file","std::string");
  m_param.put(cur);

  min_total_weight = std::max( std::min(min_total_weight , 1.0), 0.1 );
  if(cmp_lib != 0) compare_libraries=true; 
  if(deb_prop != 0) debug=true;

  std::transform(str1.begin(),str1.end(),str1.begin(),(int (*)(int)) tolower);
  if(str1 == "yes" || str1 == "true")
    set_nWalker_target = true;

  estim0 = new EstimatorHandler(myComm); 
  estim0->copyInfo(*this);
  estim0->parse(cur);

  return true;
}


bool AFQMCDriver::setup(HamPtr h0, WSetPtr w0, PropPtr p0, WfnPtr wf0)
{

  ham0=h0;
  wlkBucket=w0;
  prop0=p0;
  wfn0=wf0;
  restarted=false;

  // temporary  
/*  
  localwlkBucket = dynamic_cast<LocalWalkerHandler*>(wlkBucket);
  if(!localwlkBucket) {
    app_error()<<" Error in AFQMCDriver::setup() \n"
               <<" Conversion to LocalWalkerHandler was unsuccessful. \n\n"; 
    return false;
  }
*/
  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Beginning Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  app_log()<<" Using " <<ncores_per_TG <<" cores per node in a TaskGroup. \n";
  // right now this TG is not used. It is needed for setup purposes and to  
  // get a unique TG number for every group of cores on a node (used in the WalkerSet)
  TG.setup(ncores_per_TG,1,false);  
  std::vector<int> TGdata(5);
  TG.getSetupInfo(TGdata); 
   
  // setup local-to-node MPI Comm
  // TGdata[0]: node_number
  myComm->split_comm(TGdata[0],MPI_COMM_NODE_LOCAL);
  TG.setNodeCommLocal(MPI_COMM_NODE_LOCAL);
  int key = TG.getTGNumber(); // This works because the TG used has nnodes_per_TG=1 
  myComm->split_comm(key,MPI_COMM_TG_LOCAL);
  TG.setTGCommLocal(MPI_COMM_TG_LOCAL);
  key = TG.getCoreRank();  
  myComm->split_comm(key,MPI_COMM_TG_LOCAL_HEADS);

  CommBuffer.setup(TG.getCoreRank()==0,std::string("COMMBuffer_")+std::to_string(TG.getTGNumber()),MPI_COMM_TG_LOCAL);
  TG.setBuffer(&CommBuffer);

  hdf_archive read(myComm);
  if(myComm->rank() == 0) {
    if(hdf_read_restart != std::string("")) {

      if(read.open(hdf_read_restart,H5F_ACC_RDONLY)) 
        restarted = restart(read);

      if(!restarted) {
        read.close();
        app_log()<<" WARNING: Problems restarting simulation. Starting from default settings. \n";
      }
    
    }
  }
  myComm->bcast(restarted);
  if(restarted) {
    app_log()<<" Restarted from file. Block, step: " <<block0 <<" " <<step0 <<std::endl; 
    app_log()<<"                      Eshift, Etav: " <<Eshift <<" " <<Etav <<std::endl;
    myComm->bcast(Eshift);
    myComm->bcast(Etav);
    myComm->bcast(block0);
    myComm->bcast(step0);
  }

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Hamiltonian \n"
           <<"****************************************************\n"
           <<std::endl;

  // hamiltonian
  if(!ham0->init(TGdata,&CommBuffer,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL)) {
    app_error()<<"Error initializing Hamiltonian in AFQMCDriver::setup" <<std::endl; 
    return false; 
  }   

  app_log()<<"\n****************************************************\n"
           <<"               Initializating Wavefunction \n"
           <<"****************************************************\n"
           <<std::endl;

  if(!wfn0->init(TGdata,&CommBuffer,read,hdf_read_tag,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL)) {
    app_error()<<"Error initializing Wavefunction in AFQMCDriver::setup" <<std::endl; 
    return false; 
  }   
  if(!wfn0->setup(ham0)) {
    app_error()<<"Error in WavefunctionHandler::setup in AFQMCDriver::setup" <<std::endl; 
    return false; 
  }   

  app_log()<<"\n****************************************************\n"
           <<"             Initializating Walker Handler \n"
           <<"****************************************************\n"
           <<std::endl;

  // walker set
  wlkBucket->setup(TG.getCoreRank(),ncores_per_TG,TG.getTGNumber(),MPI_COMM_TG_LOCAL_HEADS,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,&LocalTimer);
  wlkBucket->setHF(wfn0->getHF());
  if(restarted) { 
    wlkBucket->restartFromHDF5(nWalkers,read,hdf_read_tag,set_nWalker_target);
    app_log()<<"Number of walkers after restart: " <<wlkBucket->GlobalPopulation() <<std::endl;
    wfn0->evaluateLocalEnergyAndOverlap("ImportanceSampling",-1,wlkBucket);
  } else {
    wlkBucket->initWalkers(nWalkers);
    wfn0->evaluateLocalEnergyAndOverlap("ImportanceSampling",-1,wlkBucket);
  }

  app_log()<<"\n****************************************************\n"
           <<"              Initializating Propagator \n"
           <<"****************************************************\n"
           <<std::endl;

  // propagator
  if(!prop0->setup(TGdata,&CommBuffer,ham0,wfn0,dt,read,hdf_read_tag,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL)) { 
    app_error()<<"Error in PropagatorBase::setup in AFQMCDriver::setup" <<std::endl; 
    return false; 
  }   

  if(myComm->rank() == 0) 
    read.close();

  app_log()<<"\n****************************************************\n"
           <<"              Initializating Estimators \n"
           <<"****************************************************\n"
           <<std::endl;

  // estimator setup
  estim0->setup(TGdata,&CommBuffer,ham0,wfn0,&LocalTimer,MPI_COMM_HEAD_OF_NODES,MPI_COMM_TG_LOCAL,MPI_COMM_NODE_LOCAL,MPI_COMM_TG_LOCAL_HEADS);

  app_log()<<"\n****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"          Finished Driver initialization.\n" 
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<"****************************************************\n"   
           <<std::endl;

  myComm->barrier(); 

  return true;
}

// writes checkpoint file
bool AFQMCDriver::checkpoint(int block, int step) 
{

  hdf_archive dump(myComm,false);
  if(myComm->rank() == 0) {
    std::string file;
    char fileroot[128];
    int nproc = myComm->size();
    if(hdf_write_restart != std::string("")) 
      file = hdf_write_restart;
    else
      file = myComm->getName()+std::string(".chk.h5"); 

    if(!dump.create(file)) {
      app_error()<<" Error opening checkpoint file for write. \n";
      return false;
    }

    std::vector<RealType> Rdata(2);
    Rdata[0]=Eshift;
    Rdata[1]=Etav;

    std::vector<IndexType> Idata(2);
    Idata[0]=block;
    Idata[1]=step;

    // always write driver data and walkers 
    dump.push("AFQMCDriver"); 
    if(hdf_write_tag != std::string("")) dump.push(hdf_write_tag);
    dump.write(Idata,"DriverInts");
    dump.write(Rdata,"DriverReals");
    if(hdf_write_tag != std::string("")) dump.pop();
    dump.pop();
  }

  if(!wlkBucket->dumpToHDF5(dump,hdf_write_tag) ) {
    app_error()<<" Problems writting checkpoint file in Driver/AFQMCDriver::checkpoint(). \n";
    return false;
  }

  if(myComm->rank() == 0) {
    dump.close();
  }

  return true;
}


// writes samples
bool AFQMCDriver::writeSamples()
{

  hdf_archive dump(myComm,false);
  if(myComm->rank() == 0) {
    std::string file;
    char fileroot[128];
    int nproc = myComm->size();
    file = myComm->getName()+std::string(".confg.h5");

    if(!dump.create(file)) {
      app_error()<<" Error opening checkpoint file for write. \n";
      return false;
    }

  }

  int nwtowrite=-1;
  if(!wlkBucket->dumpSamplesHDF5(dump,nwtowrite) ) {
    app_error()<<" Problems writting checkpoint file in Driver/AFQMCDriver::writeSample(). \n";
    return false;
  }

  if(myComm->rank() == 0) {
    dump.close();
  }

  return true;

}

// sets up restart archive and reads  
bool AFQMCDriver::restart(hdf_archive& read)
{

  // always write driver data and walkers 
  if(!read.push("AFQMCDriver",false)) return false; 
  if(hdf_read_tag != std::string("")) 
    if(!read.push(hdf_read_tag,false)) return false;

  std::vector<IndexType> Idata(2);
  std::vector<RealType> Rdata(2);

  if(!read.read(Idata,"DriverInts")) return false;
  if(!read.read(Rdata,"DriverReals")) return false;

  Eshift = Rdata[0];
  Etav = Rdata[1];

  block0=Idata[0];
  step0=Idata[1];

  if(hdf_read_tag != std::string("")) read.pop();
  read.pop();

  return true;
}

bool AFQMCDriver::clear()
{
  return true;
}

void AFQMCDriver::output_timers(std::ofstream& out_timers, int n)
{

  if(n==0) out_timers<<"Propagate::applyHSPropagator  Propagate::calculateMixedMatrixElementOfOneBodyOperators  Propagate::eloc  Propagate::product_SD  Propagate::sampleGaussianFields  Propagate::apply_expvHS_Ohmms  Propagate::build_vHS  PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators  PureSingleDeterminant:evaluateLocalEnergy  PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix " <<std::endl;
  out_timers<<Timer.average("Propagate::applyHSPropagator") <<" "
   <<Timer.average("Propagate::calculateMixedMatrixElementOfOneBodyOperators") <<" "
   <<Timer.average("Propagate::eloc") <<" "
   <<Timer.average("Propagate::eloc2") <<" "
   <<Timer.average("Propagate::eloc3") <<" "
   <<Timer.average("Propagate::product_SD") <<" "
   <<Timer.average("Propagate::sampleGaussianFields") <<" "
   <<Timer.average("Propagate::apply_expvHS_Ohmms") <<" "
   <<Timer.average("Propagate::build_vHS") <<" "
   <<Timer.average("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators") <<" "
   <<Timer.average("PureSingleDeterminant:evaluateLocalEnergy") <<" "
   <<Timer.average("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix") <<std::endl;
   Timer.reset("Propagate::applyHSPropagator");
   Timer.reset("Propagate::calculateMixedMatrixElementOfOneBodyOperators");
   Timer.reset("Propagate::eloc");
   Timer.reset("Propagate::eloc2");
   Timer.reset("Propagate::eloc3");
   Timer.reset("Propagate::product_SD");
   Timer.reset("Propagate::sampleGaussianFields");
   Timer.reset("Propagate::apply_expvHS_Ohmms");
   Timer.reset("Propagate::build_vHS");
   Timer.reset("PureSingleDeterminant:calculateMixedMatrixElementOfOneBodyOperators");
   Timer.reset("PureSingleDeterminant:evaluateLocalEnergy");
   Timer.reset("PureSingleDeterminant:local_evaluateOneBodyMixedDensityMatrix");

}

}
