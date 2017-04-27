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
#include "AFQMC/Drivers/VMCDriver.h"
#include "AFQMC/Estimators/SlaterDetOperations.h"

namespace qmcplusplus {

bool VMCDriver::run()
{
/*
  if(!restarted) {
    step0=block0=0;
  }
  ComplexType enume=0.0,edeno=0.0;
  std::vector<RealType> data(10);
  RealType accept=0.0; 

  ComplexType exactEstimatorEnergy;
  SlaterDetOperations SDetOps(myComm);
  SDetOps.copyInfo(*this);
  std::vector<RealType> diagEnergies(diagHam);
  ComplexMatrix diagEigVec(1); 
  if(diagHam > 0 && myComm->size() > 1 ) {
    app_error()<<" Error: Diagonalization of hamiltonian in space of walkers only implemented in serial. \n";
    return false;
  }  
  SDetOps.setup(ham0,&LocalTimer);

  prop0->SDetOps = &SDetOps; 

  std::ofstream out("vmc.dat", std::ios_base::app | std::ios_base::out);
  if(out.fail()) {
    app_error()<<"Problems opening output file.\n";
    return false;
  }
  out<<fixed;
 
  if(!restarted) { // assume a field of zero 
    ComplexMatrix Mat;
    for(WalkerHandler::WalkerIterator it=wlkBucket->begin(); it!=wlkBucket->end(); it++) 
      SDetOps.green_function((it->SlaterMat).data(),(it->SlaterMat).data()+2*NMO*NAEA,it->weight,Mat,false);
  }  

  // problems with using step_tot to do ortho and load balance
  int time = step0*nSubstep; 
  for (int iBlock=block0, step_tot=step0; iBlock<nBlock; ++iBlock) {

    LocalTimer.start("Block::TOTAL");
    enume=0.0;
    edeno=0.0;
    for (int iStep=0; iStep<nStep; ++iStep, ++step_tot) {

      // propagate
      for (int iSubstep=0; iSubstep<nSubstep; ++iSubstep,++time) {

        for(WalkerHandler::WalkerIterator it=wlkBucket->begin(); it!=wlkBucket->end(); it++) {

          if(!it->alive || std::abs(it->weight) <= 1e-6) continue;

          LocalTimer.start("SubStep::Propagate");
          prop0->Propagate(time,*it,accept);
          LocalTimer.stop("SubStep::Propagate");

        }// walkers

      }  // iSubstep

      for(WalkerHandler::WalkerIterator it=wlkBucket->begin(); it!=wlkBucket->end(); it++) {
        ComplexType hamME,ovlp;
        SDetOps.matrix_element_and_overlap((it->SlaterMat).data(),(it->SlaterMat).data()+2*NMO*NAEA,ovlp,hamME);
        register ComplexType w = ovlp/std::abs(ovlp);
        enume += w*hamME/ovlp; 
        edeno += w;
      }
    }

    data[0]=enume.real();
    data[1]=edeno.real();
    data[2]=accept;
    myComm->allreduce(data);

    if(diagHam > 0 && iBlock % diagHam_freq == 0 && iBlock > 200 ) {
      SDetOps.diag(wlkBucket->begin(),wlkBucket->end(),diagHam,diagEnergies,diagEigVec,exactEstimatorEnergy,&wfn0->getHF()); 
    }

    if(myComm->rank() == 0) {
      out<<iBlock <<" " <<std::setprecision(6) <<accept/((iBlock-block0+1)*nSubstep*nStep*myComm->size()) <<" " <<data[0]/data[1];
      if(diagHam > 0) {
        for(int i=0; i<diagHam; i++) out<<"  " <<diagEnergies[i];
        out<<"  " <<exactEstimatorEnergy.real();  
      }
      out<<" " <<LocalTimer.average("Block::TOTAL") <<std::endl;
    }

    // add estimators here
   
    // checkpoint 
    if(iBlock != 0 && iBlock % nCheckpoint == 0) 
      if(!checkpoint(iBlock,step_tot)) {
        app_error()<<" Error in VMCDriver::checkpoint(). \n" <<std::endl;
        return false;
      }
 
    LocalTimer.stop("Block::TOTAL");
  }

  app_log()<<"----------------------------------------------------------------\n";
  app_log()<<" LocalTimer: \n";
  LocalTimer.print_average_all(app_log());
  app_log()<<" Timer: \n";
  Timer.print_average_all(app_log());
  app_log()<<"----------------------------------------------------------------\n";

  out.close();
*/

  return true;
}

bool VMCDriver::parse(xmlNodePtr cur)
{
  if(cur==NULL) return false;

  std::string str;
  ParameterSet m_param;
  m_param.add(nBlock,"blocks","int");
  m_param.add(nStep,"steps","int");
  m_param.add(nSubstep,"substeps","int");
  m_param.add(nWalkers,"nWalkers","int");
  m_param.add(nCheckpoint,"checkpoint","int");
  m_param.add(dt,"dt","double");
  m_param.add(dt,"timestep","double");
  m_param.add(diagHam,"diagHam","int");
  m_param.add(hdf_read_tag,"hdf_read_tag","std::string");
  m_param.add(hdf_read_restart,"hdf_read_file","std::string");
  m_param.add(hdf_write_tag,"hdf_write_tag","std::string");
  m_param.add(hdf_write_restart,"hdf_write_file","std::string");
  m_param.put(cur);

  return true;
}


bool VMCDriver::setup(HamPtr h0, WSetPtr w0, PropPtr p0, WfnPtr wf0)
{

  ham0=h0;
  wlkBucket=w0;
  prop0=p0;
  wfn0=wf0;
  restarted=false;

  app_log()<<"\n****************************************************\n"   
           <<"          Beginning VMC Driver initialization.\n" 
           <<"****************************************************\n"   
           <<std::endl;

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

  // hamiltonian
//  if(!ham0->init()) {
//    app_error()<<"Error initializing Hamiltonian in VCDriver::setup" <<std::endl; 
//    return false; 
//  }   

  // wavefunction
  //if(!wfn0->init(read,hdf_read_tag)) {
  //  app_error()<<"Error initializing Wavefunction in VMCDriver::setup" <<std::endl; 
  //  return false; 
  //}   
  //if(!wfn0->setup(ham0)) {
  //  app_error()<<"Error in WavefunctionHandler::setup in VMCDriver::setup" <<std::endl; 
  //  return false; 
  //}   

  // walker set
//  wlkBucket->setup(,ncores_per_TG);
  // A VMC Walker has 2 states, so it is size (4*NMO,NAEA)
  ComplexMatrix HF;
  HF.resize(4*NMO,NAEA);
  for(int i=0; i<NAEA; i++) HF(i,i)=ComplexType(1.0,0.0);
  for(int i=0; i<NAEB; i++) HF(NMO+i,i)=ComplexType(1.0,0.0); 
  for(int i=0; i<NAEA; i++) HF(2*NMO+i,i)=ComplexType(1.0,0.0);
  for(int i=0; i<NAEB; i++) HF(3*NMO+i,i)=ComplexType(1.0,0.0); 
  wlkBucket->setHF(HF);
  if(restarted) { 
    wlkBucket->restartFromHDF5(nWalkers,read,hdf_read_tag,false);
  } else {
    wlkBucket->initWalkers(nWalkers);
  }

  // propagator
//  if(!prop0->setup(core_rank,ncores_per_TG,ham0,wfn0,dt,read,hdf_read_tag)) { 
//    app_error()<<"Error in PropagatorBase::setup in VMCDriver::setup" <<std::endl; 
//    return false; 
//  }   

  app_log()<<"\n****************************************************\n"   
           <<"          Finished Driver initialization.\n" 
           <<"****************************************************\n"   
           <<std::endl;

  return true;
}

// writes checkpoint file
bool VMCDriver::checkpoint(int block, int step) 
{

  hdf_archive dump(myComm,false);
  if(myComm->rank() == 0) {
    std::string file;
    char fileroot[128];
    int nproc = myComm->size();
    int nodeid = myComm->rank();
    int groupid=myComm->getGroupID();
    bool no_gtag= (qmc_common.mpi_groups==1);
    if(no_gtag)
      sprintf(fileroot,"%s.s%03d",project_title.c_str(),m_series);
    else
      sprintf(fileroot,"%s.g%03d.s%03d",project_title.c_str(),groupid,m_series);

    if(hdf_write_restart != std::string("")) 
      file = hdf_write_restart;
    else
      file = std::string(fileroot)+std::string(".chk.h5"); 

    if(!dump.create(file)) {
      app_error()<<" Error opening checkpoint file for write. \n";
      return false;
    }

    std::vector<IndexType> Idata(2);
    Idata[0]=block;
    Idata[1]=step;

    // always write driver data and walkers 
    dump.push("VMCDriver"); 
    if(hdf_write_tag != std::string("")) dump.push(hdf_write_tag);
    dump.write(Idata,"DriverInts");
    //dump.write(Rdata,"DriverReals");
    if(hdf_write_tag != std::string("")) dump.pop();
    dump.pop();
  }

  if(!wlkBucket->dumpToHDF5(dump,hdf_write_tag) ) {
    app_error()<<" Problems writting checkpoint file in Driver/VMCDriver::checkpoint(). \n";
    return false;
  }

  if(myComm->rank() == 0) {
    dump.close();
  }

  return true;
}

// sets up restart archive and reads  
bool VMCDriver::restart(hdf_archive&)
{
  return true;
}

bool VMCDriver::clear()
{
  return true;
}

}
