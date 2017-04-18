
#include<cassert>
#include<random>
#include<cstdlib>
#if defined(HAVE_MPI)
#include<mpi.h>
#endif

#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include "Message/CommOperators.h"
#include"AFQMC/Walkers/SlaterDetWalker.h"
#include"AFQMC/Walkers/DistWalkerHandler.h"

namespace qmcplusplus
{


enum DistWalkerTimers {
  LoadBalance_Setup,
  LoadBalance_Resize,
  LoadBalance_Exchange
};

TimerNameList_t<DistWalkerTimers> DistWalkerTimerNames =
{
  {LoadBalance_Setup, "WalkerHandler::loadBalance::setup"},
  {LoadBalance_Resize, "WalkerHandler::loadBalance::resize"},
  {LoadBalance_Exchange, "WalkerHandler::loadBalance::exchange"},
};

bool DistWalkerHandler::restartFromXML() 
{ 
  return true;
}

bool DistWalkerHandler::dumpToXML() 
{
  return true;
}

// not thinking about efficiency since this is only used once during initialization 
bool DistWalkerHandler::restartFromHDF5(int n, hdf_archive& read, const std::string& tag, bool set_to_target)
{
/*
//  return initWalkers(n);
  int nproc = myComm->size();
  Walker dummy;
  dummy.initWalker(HFMat,true);
  int sz = dummy.sizeForDump();
  std::vector<int> from(nproc);
  int cnt,nWtot = 0;   
  std::vector<char> buffer;
  std::vector<char> bufferall;
  std::vector<int> counts(nproc), displ(nproc);

  if(myComm->rank()==0) {

    std::vector<int> Idata(2);

    std::string path = "/Walkers/DistWalkerHandler_SlaterDetWalker";
    if(tag != std::string("")) path += std::string("/")+tag;
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group  could not find /Walkers/DistWalkerHandler_SlaterDetWalker/{tag} group in file. No restart data for walkers. \n"; 
      return false;
    }

    if(!read.push("Walkers")) return false;
    if(!read.push("DistWalkerHandler_SlaterDetWalker")) return false;
    if(tag != std::string("")) if(!read.push(tag)) return false;
    if(!read.read(Idata,"DistWalkerHandler_dims")) return false;
    nWtot=Idata[0];
    app_log()<<"Found " <<nWtot <<" walkers on restart file." <<std::endl; 
    if(Idata[1] != sz) {
      app_error()<<" ERROR: Size of walker is not consistent in hdf5 file. \n";
      app_error()<<" sz_sim, sz_file: " <<sz <<" " <<Idata[1] <<std::endl; 
      return false;
    } 
    bufferall.resize(nWtot*sz); 
    if(!read.read(bufferall,"DistWalkerHandler_walkers")) return false;
    if(tag != std::string("")) read.pop();
    read.pop();
    read.pop();


    int nWperProc = nWtot/nproc; 
    int nWextra = nWtot%nproc; 
    if( set_to_target && nWperProc >= n ) {
      nWperProc = n;
      nWextra = 0; 
    }

    for(int i=0; i<nWextra; i++) from[i] = nWperProc+1;
    for(int i=nWextra; i<nproc; i++) from[i] = nWperProc;
    myComm->bcast(from); 
    int nW = from[myComm->rank()];

    walkers.clear();
    if(set_to_target) walkers.resize(n);
    else walkers.resize(nW);
    cnt=0;
    for(WalkerIterator it=walkers.begin(); it!=walkers.begin()+nW; it++) {
      it->initWalker(HFMat);
      it->restartFromChar( bufferall.data()+cnt );
      cnt+=sz;
    }
    cnt=0;
    for(WalkerIterator it=walkers.begin()+nW; it!=walkers.end(); it++) {
      it->initWalker(HFMat);
      *it = walkers[cnt];
      cnt++;
      if(cnt == nW) cnt=0;
    }


    displ[0]=0;
    displ[1]=0;
    counts[0]=0;
    for(int i=1; i<nproc-1; i++) {
      counts[i] = from[i]*sz; 
      displ[i+1] = displ[i]+counts[i];
    }
    counts[nproc-1] = from[nproc-1]*sz; 
    myComm->scatterv( bufferall, buffer, counts, displ, 0);

  } else {

    myComm->bcast(from); 
    int nW = from[myComm->rank()];

    buffer.resize(nW*sz);
    myComm->scatterv( bufferall, buffer, counts, displ, 0);

    walkers.clear();
    if(set_to_target) walkers.resize(n);
    else walkers.resize(nW);
    cnt=0;
    for(WalkerIterator it=walkers.begin(); it!=walkers.begin()+nW; it++) {
      it->initWalker(HFMat);
      it->restartFromChar( buffer.data()+cnt );
      cnt+=sz;
    }
    cnt=0;
    for(WalkerIterator it=walkers.begin()+nW; it!=walkers.end(); it++) {
      it->initWalker(HFMat);
      *it = walkers[cnt];
      cnt++;
      if(cnt == nW) cnt=0;
    }    

  }
*/
  return true;

}

bool DistWalkerHandler::dumpToHDF5(hdf_archive& dump, const std::string& tag)
{
/*
  // check that restart data doesnot exist 
  std::string path = "/Walkers/DistWalkerHandler_SlaterDetWalker";
  if(tag != std::string("")) path += std::string("/")+tag;
  if(dump.is_group( path )) {
    app_error()<<" ERROR: H5Group /Walkers/DistWalkerHandler_SlaterDetWalker/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
    return false;
  }

  // doing one big gatherV, not sure if this is the best way right now 
  int nW = numWalkers();
  int sz = walkers[0].sizeForDump(); 
  std::vector<int> to(1,nW); 
  std::vector<int> from(myComm->size());
  int nWtot = 0;

  myComm->allgather(to,from,1);
  for(int i=0; i<myComm->size(); i++) nWtot += from[i];
  
  std::vector<char> buffer(nW*sz);
  std::vector<char> bufferall;
  if(myComm->rank() == 0) 
    bufferall.resize(nWtot*sz); 

  int cnt=0;
  for(WalkerIterator it=walkers.begin(); it!=walkers.end(); it++) 
    if(it->alive) {
      it->dumpToChar( buffer.data()+cnt );
      cnt+=sz;
    }

  std::vector<int> displ(myComm->size());
  for(int i=0; i<myComm->size(); i++) from[i]*=sz;
  cnt=0;
  for(int i=0; i<myComm->size(); i++) {
    displ[i] = cnt; 
    cnt+=from[i];
  }
  myComm->gatherv(buffer,bufferall,from,displ,0);

  if(myComm->rank()==0) {
    // now write to HDF5 file
    std::vector<int> Idata(2);
    Idata[0]=nWtot;
    Idata[1]=sz;

    dump.push("Walkers");
    dump.push("DistWalkerHandler_SlaterDetWalker");
    if(tag != std::string("")) dump.push(tag);
    dump.write(Idata,"DistWalkerHandler_dims");
    dump.write(bufferall,"DistWalkerHandler_walkers");
    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    dump.flush();  

  }
*/
  return true;
}

bool DistWalkerHandler::parse(xmlNodePtr cur)
{

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);
    walkerType = "collinear";

    load_balance_alg = "all";

    ParameterSet m_param;
    m_param.add(reset_weight,"reset_weight","double");
    m_param.add(max_weight,"max_weight","double");
    m_param.add(min_weight,"min_weight","double");
    m_param.add(extra_empty_spaces,"extra_spaces","int");
    m_param.add(walkerType,"walker_type","std::string");
    m_param.add(load_balance_alg,"algo","std::string");
    m_param.put(cur);

    if(load_balance_alg == "all")
      app_log()<<" Using all-to-all (gather/scatter) load balancing algorithm. " <<std::endl;
    else if(load_balance_alg == "iter" || load_balance_alg == "iterative" || load_balance_alg == "seq") {
      load_balance_alg = "iter";
      app_log()<<" Using iterative (1-to-1 communication based) load balancing algorithm. " <<std::endl;
    } else {
      std::cerr<<" Error: Unknown load balancing algorithm: " <<load_balance_alg <<std::endl;
      return false;  
    }  

    std::transform(walkerType.begin(),walkerType.end(),walkerType.begin(),(int (*)(int)) tolower);
    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;
}

void DistWalkerHandler::setHF(const ComplexMatrix& HF) 
{ 
  assert( HF.rows() == HFMat.rows() && HF.cols() == HFMat.cols() );
  HFMat = HF;
}

bool DistWalkerHandler::setup(int cr, int nc, int tgn, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, myTimer* timer)
{
  LocalTimer=timer;
  setup_timers(Timers, DistWalkerTimerNames);

  core_rank=cr;
  ncores_per_TG=nc; 
  head = core_rank==0;
  ncol=NAEA;
  nrow=NMO;
  if(walkerType == "closed_shell")  { 
    type = 1;
  } else if(walkerType == "collinear") {  
    type = 2;
    nrow*=2;
  } else if(walkerType == "non-collinear") {   
    type = 4;
    nrow*=2;
    ncol*=2;
  } else {
    app_error()<<" Error: Incorrect walker type on DistWalkerHandler::setup \n";
    return false; 
  }

  // Stored data [all assumed std::complex numbers]:
  //   - INFO:                 1  (e.g. alive, init, etc)  
  //   - SlaterMatrix:         NCOL*NROW 
  //        type = 1 for closed_shell
  //        type = 2 for non-closed shell collinear
  //        type = 4 for non-collinear
  //   - weight:               1 
  //   - eloc:                 2
  //   - eloc_old:             2
  //   - overlap_alpha:        2
  //   - overlap_beta:         2
  //   - old_overlap_alpha:    2
  //   - old_overlap_beta:     2
  //   Total: 14+NROW*NCOL
  int cnt=0;
  data_displ[0] = cnt; cnt++;
  data_displ[1] = cnt; cnt+=nrow*ncol;
  data_displ[2] = cnt; cnt+=1; // weight 
  data_displ[3] = cnt; cnt+=2; // eloc
  data_displ[4] = cnt; cnt+=2; // eloc_old
  data_displ[5] = cnt; cnt+=2; // overlap_alpha
  data_displ[6] = cnt; cnt+=2; // overlap_beta
  data_displ[7] = cnt; cnt+=2; // old_overlap_alpha
  data_displ[8] = cnt; cnt+=2; // old_overlap_beta
  walker_size = cnt;  
  walker_memory_usage = walker_size*sizeof(ComplexType);

  walkers.setup(head,std::string("DistWalkerHandler_")+std::to_string(tgn),tg_comm);
  tot_num_walkers=maximum_num_walkers=0;

  MPI_COMM_TG_LOCAL = tg_comm;
  MPI_COMM_TG_LOCAL_HEADS = heads_comm;
#if defined(HAVE_MPI)
  MPI_Comm_size(MPI_COMM_TG_LOCAL_HEADS,&nproc_heads);
  MPI_Comm_rank(MPI_COMM_TG_LOCAL_HEADS,&rank_heads);
#else
  nproc_heads=rank_heads=0;
#endif
  empty_spots.reserve(1000);
  outgoing.reserve(nproc_heads);
  incoming.reserve(nproc_heads); 
  counts.resize(nproc_heads);
  displ.resize(nproc_heads);
  nwalk_counts_old.resize(nproc_heads);
  nwalk_counts_new.resize(nproc_heads);
  
  // default setup
  HFMat.resize(nrow,ncol);
  for(int i=0; i<nrow; i++)
   for(int j=0; j<ncol; j++)
     HFMat(i,j) = ComplexType(0.0,0.0);

  for(int j=0; j<ncol; j++)
   HFMat(j,j) = ComplexType(1.0,0.0);

  if(type == 2)
   for(int j=0; j<NAEB; j++)
    HFMat(j+NMO,j) = ComplexType(1.0,0.0);

  // setup buffers for dump and communication 

  min_weight = std::max(std::abs(min_weight),1e-2); 
  return true;
}

bool DistWalkerHandler::clean()
{
  walkers.clear();
  return true;
}

// called at the beginning of each executable section
void DistWalkerHandler::resetNumberWalkers(int n, bool a, ComplexMatrix* S) 
{

  if(head) {

    ComplexMatrix* S0 = S;
    if(S0==NULL) S0 = &HFMat;
    // assuming everything is alive
    int ns = size(); 
    if(n+extra_empty_spaces > ns) {
      std::vector<ComplexType> tmp(tot_num_walkers*walker_size); // store old walker info
      int cnt=0;
      for(int i=0; i<tot_num_walkers; i++)
        if(walkers[walker_size*i+data_displ[INFO]].real() > 0) 
          std::copy(walkers.begin()+walker_size*i,walkers.begin()+walker_size*(i+1),tmp.begin()+(cnt++)*walker_size);
      walkers.resize((n+extra_empty_spaces)*walker_size);
      tot_num_walkers = std::min(n,cnt); 
      std::copy(tmp.begin(),tmp.begin()+walker_size*tot_num_walkers,walkers.begin());
      for(int i=tot_num_walkers; i<n+extra_empty_spaces; i++)
        walkers[walker_size*i+data_displ[INFO]] = ComplexType(-1.0);   
    }


    // now I have enough space, add or remove walkers if needed 
    if(tot_num_walkers < n) {
      // adding walkers
      int cnt=tot_num_walkers;  
      for(int i=0; i<size(); i++) {
        if(walkers[walker_size*i+data_displ[INFO]].real() < 0) {   
          walkers[walker_size*i+data_displ[INFO]] = ComplexType(1.0);
          std::copy(S0->begin(),S0->end(),walkers.begin()+walker_size*i+data_displ[SM]);
          *(walkers.begin()+walker_size*i+data_displ[WEIGHT]) = ComplexType(1.0); 
          *(walkers.begin()+walker_size*i+data_displ[ELOC]) = ComplexType(0); 
          *(walkers.begin()+walker_size*i+data_displ[OVLP_A]) = ComplexType(0); 
          *(walkers.begin()+walker_size*i+data_displ[OVLP_B]) = ComplexType(0); 
          cnt++;
        }
        if(cnt == n) break;
      }
      if(cnt != n) APP_ABORT("Error: Problems in DistWalkerHandler::resetNumberWalkers-add. \n");
    } else if(tot_num_walkers > n) { 
      // removing walkers
      int cnt=tot_num_walkers;
      for(int i=size()-1; i>=0; i--) {
        if(walkers[walker_size*i+data_displ[INFO]].real() > 0) {  
          walkers[walker_size*i+data_displ[INFO]] = ComplexType(-1.0);
          cnt--;
        }
        if(cnt == n) break;
      }
      if(cnt != n) APP_ABORT("Error: Problems in DistWalkerHandler::resetNumberWalkers-remove. \n");
    }


  } else {
    int ns = size(); 
    if(n+extra_empty_spaces > ns) 
      walkers.resize((n+extra_empty_spaces)*walker_size);
  }
  walkers.barrier(); 
  maximum_num_walkers = n+extra_empty_spaces;
  tot_num_walkers = n; 
}

// load balancing algorithm
void DistWalkerHandler::loadBalance() 
{
#if defined(HAVE_MPI)
  
  nwalk_min = nwalk_max= tot_num_walkers;
  if(nproc_heads==1) return;

  int nW = numWalkers();
  int sz = size(), nw_new=0;
  MPI_Request request;
 
  LocalTimer->start("WalkerHandler::loadBalance::setup");
  Timers[LoadBalance_Setup]->start();
  // determine new number of walkers
  if(head) {

    for(int i=0; i<nproc_heads; i++)
      nwalk_counts_old[i] = nwalk_counts_new[i] = 0;
    nwalk_counts_new[0] = nW;
    // in case iallgather is not available
    //myComm->allgather(nwalk_counts_new,nwalk_counts_old,1,MPI_COMM_TG_LOCAL_HEADS);
#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER
#endif
#ifdef HAVE_MPI_IALLGATHER
    MPI_Iallgather(nwalk_counts_new.data(), 1, MPI_INT, nwalk_counts_old.data(), 1, MPI_INT, MPI_COMM_TG_LOCAL_HEADS,&request); 

    // push empty spots to the end of the list 
    push_walkers_to_front();

    // wait for mpi_iallgather 
    MPI_Wait(&request, MPI_STATUS_IGNORE); 
#else
    MPI_Allgather(nwalk_counts_new.data(), 1, MPI_INT, nwalk_counts_old.data(), 1, MPI_INT, MPI_COMM_TG_LOCAL_HEADS); 
    push_walkers_to_front();
#endif
      
    nwalk_global=0;
    for(int i=0; i<nproc_heads; i++) 
      nwalk_global+=nwalk_counts_old[i];
    for(int i=0; i<nproc_heads; i++) 
      nwalk_counts_new[i] = nwalk_global/nproc_heads + ((i<nwalk_global%nproc_heads)?(1):(0));
    auto min_max = std::minmax_element(nwalk_counts_old.begin(),nwalk_counts_old.end());
    nwalk_min = *min_max.first;
    nwalk_max = *min_max.second;

    nw_new = nwalk_counts_new[rank_heads];

//if(rank()==0)
//  for(int i=0; i<nproc_heads; i++)
//    app_log()<<i <<" " <<nwalk_counts_old[i] <<" " <<nwalk_min <<" " <<nwalk_max <<" " <<nw_new <<std::endl; 

//app_log()<<" Before: " <<std::endl;
//cout<<rank() <<" " <<nW <<" " <<nw_new <<" " <<nwalk_global <<std::endl;
//MPI_Barrier(MPI_COMM_TG_LOCAL_HEADS);

  }
  Timers[LoadBalance_Setup]->stop();
  LocalTimer->stop("WalkerHandler::loadBalance::setup");

  Timers[LoadBalance_Setup]->start();
  LocalTimer->start("WalkerHandler::loadBalance::resize");
  // resize arrays if necessary. This requires all cores in a TG 
  walkers.share(&nw_new,1,head);
  if(nw_new > sz) walkers.resize(walker_size*(nw_new+extra_empty_spaces));
  LocalTimer->stop("WalkerHandler::loadBalance::resize");
  Timers[LoadBalance_Resize]->stop();

  Timers[LoadBalance_Exchange]->start();
  LocalTimer->start("WalkerHandler::loadBalance::exchange");
  if(load_balance_alg == "all") {

    if(head) {

      int ncomm = nW - nwalk_counts_new[rank_heads];
      int pos = 0, cnt=0;
      char *buffer; 
      if(ncomm > 0) 
        buffer = reinterpret_cast<char*>( &( *(walkers.begin() + walker_size*nwalk_counts_new[rank_heads]) ) );
       else 
        buffer = reinterpret_cast<char*>( &( *(walkers.begin() + walker_size*nW) ) );
      

      if(rank_heads==0) {

        // setup gatherv call 
        displ[0]=0;
        for(int i=0; i<nproc_heads-1; i++)
          if(nwalk_counts_old[i] > nwalk_counts_new[i]) {
            counts[i] = (nwalk_counts_old[i]-nwalk_counts_new[i])*walker_size;
            displ[i+1] = displ[i]+counts[i];
          } else {
            counts[i]=0;
            displ[i+1] = displ[i];
          }
        if(nwalk_counts_old[nproc_heads-1] > nwalk_counts_new[nproc_heads-1]) 
          counts[nproc_heads-1] = (nwalk_counts_old[nproc_heads-1]-nwalk_counts_new[nproc_heads-1])*walker_size;
        else
          counts[nproc_heads-1]=0;
        bufferall.clear();
        bufferall.resize(displ[nproc_heads-1]+counts[nproc_heads-1]);
        int nn = walker_size*((ncomm>0)?ncomm:0); 
        myComm->gatherv( buffer, bufferall.data(), nn , counts, displ, 0, MPI_COMM_TG_LOCAL_HEADS); 
  
        // setup scatterv call
        displ[0]=0;
        for(int i=0; i<nproc_heads-1; i++)
          if(nwalk_counts_old[i] < nwalk_counts_new[i]) {
            counts[i] = (nwalk_counts_new[i]-nwalk_counts_old[i])*walker_size;
            displ[i+1] = displ[i]+counts[i];
          } else {
            counts[i]=0;
            displ[i+1] = displ[i];
          }
        if(nwalk_counts_old[nproc_heads-1] < nwalk_counts_new[nproc_heads-1])
          counts[nproc_heads-1] = (nwalk_counts_new[nproc_heads-1]-nwalk_counts_old[nproc_heads-1])*walker_size;
        else
          counts[nproc_heads-1]=0;    

        nn = walker_size*((ncomm<0)?(std::abs(ncomm)):0); 
        myComm->scatterv( bufferall.data(), buffer, nn, counts, displ, 0, MPI_COMM_TG_LOCAL_HEADS); 

        if(ncomm > 0) {
          register ComplexType zero = ComplexType(-1.0,0.0);
          for(ComplexSMVector::iterator it=walkers.begin()+nwalk_counts_new[rank_heads]*walker_size; it<walkers.end(); it+=walker_size)
            *(it+data_displ[INFO]) = zero; 
        }

      } else {

        int nn = walker_size*((ncomm>0)?ncomm:0);
        myComm->gatherv( buffer, bufferall.data(), nn, counts, displ, 0, MPI_COMM_TG_LOCAL_HEADS);
        nn = walker_size*((ncomm<0)?(std::abs(ncomm)):0);
        myComm->scatterv( bufferall.data(), buffer, nn, counts, displ, 0, MPI_COMM_TG_LOCAL_HEADS); 

        if(ncomm > 0) {
          register ComplexType zero = ComplexType(-1.0,0.0);
          for(ComplexSMVector::iterator it=walkers.begin()+nwalk_counts_new[rank_heads]*walker_size; it<walkers.end(); it+=walker_size)
            *(it+data_displ[INFO]) = zero; 
        }

      }
    } // head 

  } else if(load_balance_alg == "iter") {
  // 1. (only heads) exchange data directly into SM without need for temporary
  //  - Implement both 1. stupid gather/scatter algorithm
  //                   2. clever algorithm based on 1-1 communications

    

  }
  LocalTimer->stop("WalkerHandler::loadBalance::exchange");
  Timers[LoadBalance_Exchange]->stop();

  walkers.barrier();
  tot_num_walkers=0;
  for(ComplexSMVector::iterator it=walkers.begin()+data_displ[INFO]; it<walkers.end(); it+=walker_size)
    if(it->real() > 0) tot_num_walkers++;

//if(head) {
//app_log()<<" After: " <<std::endl;
//cout<<rank() <<" " <<tot_num_walkers <<std::endl;
//MPI_Barrier(MPI_COMM_TG_LOCAL_HEADS);
//}
#endif
} 

// moves all walkers to the front and empty spots to the end 
void DistWalkerHandler::push_walkers_to_front()
{
  if(!head) return;
  ComplexSMVector::iterator empty = walkers.begin(), wlk = walkers.end()-walker_size;

  // 0. while wlk != empty
  // 1. find next walker from the end
  // 2. find next empty from front
  // 3. swap 

  while( wlk < empty ) {

    // 1. find next walker 
    while( (wlk+data_displ[INFO])->real() < 0 && wlk > empty ) 
      wlk -= walker_size;
    if(wlk <= empty)
      return;

    // 2. find next walker
    while( (empty+data_displ[INFO])->real() > 0 && empty < wlk ) 
      empty += walker_size;
    if(wlk <= empty)
      return;

    // 3. swap
    std::copy(wlk,wlk+walker_size,empty); 
    *(wlk+data_displ[INFO]) = ComplexType(-1.0,0.0); 
  }

}

// population control algorithm
void DistWalkerHandler::popControl()
{
  ComplexType minus = ComplexType(-1.0,0.0);

  // do I need to sync here?
  walkers.barrier();
  int buff_reset = 0; 
  if( head ) { 
    int max=size(), nw=0, cnt=0;   
    empty_spots.clear();
    for(ComplexSMVector::iterator it=walkers.begin(); it<walkers.end(); it+=walker_size, cnt++) {
      if( (it+data_displ[INFO])->real() < 0 || std::abs(*(it+data_displ[WEIGHT])) <= 1e-6) {
        // walker is not alive 
        empty_spots.push_back(cnt);  
        *(it+data_displ[INFO]) = minus;
      } else {
        ComplexType w0 = *(it+data_displ[WEIGHT]); 
        if( std::abs(w0) < std::abs(min_weight)) {
          if( (int)(distribution(generator) + std::abs(w0)/reset_weight) == 0 ) { 
            *(it+data_displ[INFO]) = minus; 
            empty_spots.push_back(cnt);
          } else {
            *(it+data_displ[WEIGHT]) = reset_weight;
          } 
        }
      }
    }

    // number of walkers after killing small weights 
    nw = max - empty_spots.size();

    // tentative algorithm to avoid multiple memory reallocations
    //  1. go through the list and decide how many copies of a walker are needed. Use upper bound 
    //  2. Expand memory allocation if necessary
    //  3. Make copies of walkers 
    int num_new=0;
    for(ComplexSMVector::iterator it=walkers.begin(); it<walkers.end(); it+=walker_size) 
      if( (it+data_displ[INFO])->real() > 0 && std::abs(*(it+data_displ[WEIGHT])) > std::abs(max_weight)) 
        num_new += std::floor( std::abs( *(it+data_displ[WEIGHT]) ) / std::abs(reset_weight) );  // since I only need n-1 copies 

    if( num_new > empty_spots.size() ) {
      buff_reset = 1; // reset SHM pointers when done
      int ntot = nw + num_new + extra_empty_spaces; 
      std::vector<ComplexType> tmp(nw*walker_size); // store old walker info
      cnt=0;
      for(ComplexSMVector::iterator it=walkers.begin(); it<walkers.end(); it+=walker_size) 
        if( (it+data_displ[INFO])->real() > 0 )
          std::copy( it, it+walker_size, tmp.begin()+(cnt++)*walker_size);
      int nn = ntot*walker_size;
      walkers.share(&nn,1,head);
      walkers.resize(ntot*walker_size,false);
      std::copy(tmp.begin(),tmp.begin()+walker_size*nw,walkers.begin());
      for(ComplexSMVector::iterator it=walkers.begin()+nw*walker_size; it<walkers.end(); it+=walker_size) 
        *(it+data_displ[INFO]) = minus; 

      empty_spots.clear();
      int newsz = size();
      empty_spots.reserve( newsz-nw );
      for(int i=nw; i<newsz; i++) empty_spots.push_back(i);
    } else {
      int nn=0;
      walkers.share<int>(&nn,1,head);
    }

    // insert new walkers in empty spots
    cnt=0; 
    for(ComplexSMVector::iterator it=walkers.begin(); it<walkers.end(); it+=walker_size) 
      if( (it+data_displ[INFO])->real() > 0 && std::abs(*(it+data_displ[WEIGHT])) > std::abs(max_weight)) {   
        RealType w = std::abs(*(it+data_displ[WEIGHT]));
        int n = (int) (w/std::abs(reset_weight));
        RealType rem = w-n*std::abs(reset_weight);
        if( ( (int)(distribution(generator) + std::abs(rem/reset_weight) ) ) != 0 ) n++;      
        *(it+data_displ[WEIGHT]) *= reset_weight/w; 
        for(int i=0; i<n-1; i++,cnt++) 
          std::copy( it, it+walker_size, walkers.begin()+walker_size*empty_spots[cnt] ); 
      }
  
    empty_spots.clear();
     
  } else {
    int nn; 
    walkers.share<int>(&nn,1,head);
    if(nn > 0)
      walkers.resize(nn,false);
  }
  walkers.barrier();
  tot_num_walkers=0;
  for(ComplexSMVector::iterator it=walkers.begin()+data_displ[INFO]; it<walkers.end(); it+=walker_size)
    if(it->real() > 0) tot_num_walkers++;
}

}

