
#include<cassert>
#include<random>
#include<cstdlib>

#include "Configuration.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include "Message/CommOperators.h"
#include"AFQMC/Walkers/SlaterDetWalker.h"
#include"AFQMC/Walkers/LocalWalkerHandler.h"

namespace qmcplusplus
{

bool LocalWalkerHandler::restartFromXML() 
{ 
  return true;
}

bool LocalWalkerHandler::dumpToXML() 
{
  return true;
}

// not thinking about efficiency since this is only used once during initialization 
bool LocalWalkerHandler::restartFromHDF5(int n, hdf_archive& read, const std::string& tag, bool set_to_target)
{

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

    std::string path = "/Walkers/LocalWalkerHandler";
    if(tag != std::string("")) path += std::string("/")+tag;
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group  could not find /Walkers/LocalWalkerHandler/{tag} group in file. No restart data for walkers. \n"; 
      return false;
    }

    if(!read.push("Walkers")) return false;
    if(!read.push("LocalWalkerHandler")) return false;
    if(tag != std::string("")) if(!read.push(tag)) return false;
    if(!read.read(Idata,"LocalWalkerHandler_dims")) return false;
    nWtot=Idata[0];
    app_log()<<"Found " <<nWtot <<" walkers on restart file." <<std::endl; 
    if(Idata[1] != sz) {
      app_error()<<" ERROR: Size of walker is not consistent in hdf5 file. \n";
      app_error()<<" sz_sim, sz_file: " <<sz <<" " <<Idata[1] <<std::endl; 
      return false;
    } 
    bufferall.resize(nWtot*sz); 
    if(!read.read(bufferall,"LocalWalkerHandler_walkers")) return false;
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

  return true;

}

bool LocalWalkerHandler::dumpToHDF5(hdf_archive& dump, const std::string& tag)
{

  // check that restart data doesnot exist 
  std::string path = "/Walkers/LocalWalkerHandler";
  if(tag != std::string("")) path += std::string("/")+tag;
  if(dump.is_group( path )) {
    app_error()<<" ERROR: H5Group /Walkers/LocalWalkerHandler/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
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
    dump.push("LocalWalkerHandler");
    if(tag != std::string("")) dump.push(tag);
    dump.write(Idata,"LocalWalkerHandler_dims");
    dump.write(bufferall,"LocalWalkerHandler_walkers");
    if(tag != std::string("")) dump.pop();
    dump.pop();
    dump.pop();

    dump.flush();  

  }

  return true;
}

bool LocalWalkerHandler::parse(xmlNodePtr cur)
{

    if(cur == NULL)
      return false;

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);

    ParameterSet m_param;
    m_param.add(reset_weight,"reset_weight","double");
    m_param.add(max_weight,"max_weight","double");
    m_param.add(min_weight,"min_weight","double");
    m_param.add(extra_empty_spaces,"extra_spaces","int");
    m_param.put(cur);

    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }

    return true;
}

void LocalWalkerHandler::setHF(const ComplexMatrix& HF) 
{ 
  //assert( (HF.rows() == 2*NMO || HF.rows() == 4*NMO) && HF.cols()==NAEA );
  HFMat.resize(HF.rows(),HF.cols());
  HFMat = HF;
}

bool LocalWalkerHandler::setup(int cr, int nc, int tgn, MPI_Comm heads_comm, MPI_Comm tg_comm, MPI_Comm node_comm, myTimer* timer)
{

  LocalTimer=timer;
  if(cr != 0) {
    app_error()<<" Error: Found core_rank != 0 in LocalWalkerHandler::setup. \n";  
    APP_ABORT(" Error: Found ncores_per_TG != 1 in LocalWalkerHandler::setup. \n");
    return false;
  }    
  if(nc != 1) {
    app_error()<<" Error: Found ncores_per_TG != 1 in LocalWalkerHandler::setup. \n";  
    APP_ABORT(" Error: Found ncores_per_TG != 1 in LocalWalkerHandler::setup. \n");
    return false;
  }    
  ncores_per_TG=1;
  core_rank=0;
  // default setup
  HFMat.resize(2*NMO,NAEA);
  for(int i=0; i<2*NMO; i++)
   for(int j=0; j<NAEA; j++)
     HFMat(i,j) = ComplexType(0.0,0.0);

  for(int j=0; j<NAEA; j++)
   HFMat(j,j) = ComplexType(1.0,0.0);

  for(int j=0; j<NAEB; j++)
   HFMat(j+NMO,j) = ComplexType(1.0,0.0);
  
  // walkerSizeForCommunication

  // walkerSizeForDump 

  // setup buffers for dump and communication 

  min_weight = std::max(std::abs(min_weight),1e-2); 
  return true;
}

bool LocalWalkerHandler::clean()
{
  for(int i=0; i<walkers.size(); i++) 
    walkers[i].clear();
  return true;
}

// called at the beginning of each executable section
void LocalWalkerHandler::resetNumberWalkers(int n, bool a, ComplexMatrix* S) 
{
  ComplexMatrix* S0 = S;
  if(S0==NULL) S0 = &HFMat;
  // assuming everything is alive
  int nold = walkers.size();
  if(nold > n ) {
     walkers.erase(walkers.begin()+n,walkers.end());
  } else {
     walkers.resize(n);
     for(int i=nold; i<n; i++)
       walkers[i].initWalker(*S0,true);;
  }  


}

/*
void LocalWalkerHandler::resetNumberWalkers(int n, bool a, ComplexMatrix* S) 
{
  ComplexMatrix* S0 = S;
  if(S0==NULL) S0 = &HFMat;
  int nold = numWalkers(); 
  int old_size = walkers.size();
  assert( walkers.size() == nold+emptyWalkers.size() );
  if(nold == n) {
    // still resize to add empty spaces if needed
    if(emptyWalkers.size() >= extra_empty_spaces) return;   
    int add = extra_empty_spaces - emptyWalkers.size();
    emptyWalkers.reserve(2*extra_empty_spaces);
// careful with c++11 move semantics
    walkers.resize(old_size+add); 
    for(int i=old_size; i<walkers.size(); i++) { 
      walkers[i].initWalker(*S0,false); 
      emptyWalkers.push_back(i);
    }
    // sort in decreasing order to keep list "tight" 
    std::sort(emptyWalkers.begin(),emptyWalkers.end(),sortDecreasing);
    return;
  } 
  if(nold == 0) {
    int add = n+extra_empty_spaces - old_size;
    emptyWalkers.reserve(2*extra_empty_spaces);
    if(add > 0) {
      walkers.resize(old_size+add); 
      for(int i=old_size; i<walkers.size(); i++) { 
        walkers[i].initWalker(*S0,false); 
        emptyWalkers.push_back(i);
      }
    }
    // sort in decreasing order to keep list "tight" 
    std::sort(emptyWalkers.begin(),emptyWalkers.end(),sortDecreasing);
    for(int i=0; i<n; i++) {
      int pos = emptyWalkers.back();
      emptyWalkers.pop_back();
      walkers[pos].initWalker(*S0,a);
    }
    return;
  }
  if( nold > n ) {
    int res=0;
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive)
        if(res < n) res++;
        else {
          res++;
          it->weight=RealType(0.0);
          //it->weight=ComplexType(0.0);
          it->alive=false; 
          emptyWalkers.push_back(std::distance(walkers.begin(),it)); 
        } 
    int add = n+extra_empty_spaces - old_size;
    emptyWalkers.reserve(2*extra_empty_spaces);
    if(add > 0) {
      walkers.resize(old_size+add);
      for(int i=old_size; i<walkers.size(); i++) {
        walkers[i].initWalker(*S0,false);  
        emptyWalkers.push_back(i);
      }
    }
    // sort in decreasing order to keep list "tight" 
    std::sort(emptyWalkers.begin(),emptyWalkers.end(),sortDecreasing);
    return;
  } 
  if(nold < n) {
    int add = n+extra_empty_spaces - old_size;
    emptyWalkers.reserve(2*extra_empty_spaces);
    if(S == NULL) {
      for(int i=0; i<walkers.size(); i++)
        if(walkers[i].alive) {
          S0 = &(walkers[i].SlaterMat);
          break;
        } 
    }
    if(add > 0) {
      walkers.resize(old_size+add);
      for(int i=old_size; i<walkers.size(); i++) {
        walkers[i].initWalker(*S0,false);
        emptyWalkers.push_back(i);
      }
    }
    // sort in decreasing order to keep list "tight" 
    std::sort(emptyWalkers.begin(),emptyWalkers.end(),sortDecreasing);
    for(int i=0; i<n-nold; i++) {
      int pos = emptyWalkers.back();
      emptyWalkers.pop_back();
      walkers[pos].initWalker(*S0,a);
    }
    return;
  } 
}  
*/
// load balancing algorithm
void LocalWalkerHandler::loadBalance(MPI_Comm comm) 
{

  int rank = myComm->rank(); 
  int nproc = myComm->size(); 
  if(nproc==1) return;
  int nW = numWalkers();
  int sz = walkers[0].sizeForDump(); 
  nwalk_counts_old.resize(nproc);
  nwalk_counts_new.resize(nproc);
  nwalk_global=0;
  std::vector<int> counts(nproc),displ(nproc);

  for(int i=0; i<nproc; i++)
    nwalk_counts_old[i] = nwalk_counts_new[i] = 0;
  nwalk_counts_new[0] = nW;
  myComm->allgather(nwalk_counts_new,nwalk_counts_old,1); 
  for(int i=0; i<nproc; i++) 
    nwalk_global+=nwalk_counts_old[i];
  for(int i=0; i<nproc; i++) 
    nwalk_counts_new[i] = nwalk_global/nproc + ((i<nwalk_global%nproc)?(1):(0));
  auto min_max = std::minmax_element(nwalk_counts_old.begin(),nwalk_counts_old.end());
  nwalk_min = *min_max.first;
  nwalk_max = *min_max.second;
/*
if(rank==0) {
app_log()<<"min,max: " <<nwalk_min <<" " <<nwalk_max <<std::endl;
for(int i=0; i<nproc; i++)
app_log()<<nwalk_counts_old[i] <<" ";
app_log()<<std::endl;  
for(int i=0; i<nproc; i++)
app_log()<<nwalk_counts_new[i] <<" ";
app_log()<<std::endl;  
}
*/
  std::vector<char> buffer;
  std::vector<char> bufferall;
  int ncomm = nW - nwalk_counts_new[rank];
  int pos = 0, cnt=0;
  buffer.reserve(std::abs(sz*ncomm));
  if(ncomm > 0) {
    buffer.resize(sz*ncomm);
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it) {
      if(it->alive) {
        cnt++;
        if(cnt > nwalk_counts_new[rank]) {
          it->dumpToChar(buffer.data()+pos);
          pos+=sz;
          it->alive = false;
          it->weight = ComplexType(0,0);
        } 
      } 
    }
    // until I code the use of it->alive properly
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ) { 
      if (!it->alive ) it=walkers.erase(it);
      else ++it;
    }
  }

  if(rank==0) {

    // setup gatherv call 
    displ[0]=0;
    for(int i=0; i<nproc-1; i++)
      if(nwalk_counts_old[i] > nwalk_counts_new[i]) {
        counts[i] = (nwalk_counts_old[i]-nwalk_counts_new[i])*sz;
        displ[i+1] = displ[i]+counts[i];
      } else {
        counts[i]=0;
        displ[i+1] = displ[i];
      }
    if(nwalk_counts_old[nproc-1] > nwalk_counts_new[nproc-1]) 
      counts[nproc-1] = (nwalk_counts_old[nproc-1]-nwalk_counts_new[nproc-1])*sz;
    else
      counts[nproc-1]=0;
    bufferall.resize(displ[nproc-1]+counts[nproc-1]);
    myComm->gatherv( buffer, bufferall, counts, displ, 0); 

    // setup scatterv call
    displ[0]=0;
    for(int i=0; i<nproc-1; i++)
      if(nwalk_counts_old[i] < nwalk_counts_new[i]) {
        counts[i] = (nwalk_counts_new[i]-nwalk_counts_old[i])*sz;
        displ[i+1] = displ[i]+counts[i];
      } else {
        counts[i]=0;
        displ[i+1] = displ[i];
      }
    if(nwalk_counts_old[nproc-1] < nwalk_counts_new[nproc-1])
      counts[nproc-1] = (nwalk_counts_new[nproc-1]-nwalk_counts_old[nproc-1])*sz;
    else
      counts[nproc-1]=0;    

    if(ncomm < 0) 
      buffer.resize(sz*std::abs(ncomm));    
    else 
      buffer.resize(0);    

    myComm->scatterv( bufferall, buffer, counts, displ, 0); 

  } else {
    myComm->gatherv( buffer, bufferall, counts, displ, 0); 
    if(ncomm < 0) 
      buffer.resize(sz*std::abs(ncomm));    
    else 
      buffer.resize(0);    
    myComm->scatterv( bufferall, buffer, counts, displ, 0); 
  }

  if(ncomm >= 0) return;
  ncomm = std::abs(ncomm);

  int oldsz = walkers.size();
  walkers.resize(oldsz+ncomm);

  pos=0;
  for(int i=0; i<ncomm; i++) {
    walkers[oldsz+i].initWalker(HFMat,true);
    walkers[oldsz+i].restartFromChar(buffer.data()+pos);
    pos+=sz;
  }

/*

  if(walkers.size() < nwalk_counts_new[rank]) {
     int oldSz = walkers.size();
     walkers.resize(nwalk_counts_new[rank]+extra_empty_spaces);
     for(int i=oldSz; i<walkers.size(); i++)
       walkers[i].initWalker(HFMat,false);
  } 

  cnt=0;
  pos=0;
  for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it) {
    if(!it->alive) {
      it->restartFromChar(buffer.data()+pos);
      pos+=sz;
      cnt++;
      it->alive=true;
    }
    if(cnt==ncomm) return;
  }
*/

} 

// population control algorithm
void LocalWalkerHandler::popControl(MPI_Comm comm, std::vector<ComplexType>& curData)
{

  for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ) 
    if (!it->alive) 
     it=walkers.erase(it);
    else
     ++it;

  // handle small weights first 
  for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ) 
    // if below cutoff, flip a coin and kill walker if necessary
    // check if this is correct, do I reset the weight to min_weight if successful??? 
    if( it->alive && std::abs(it->weight) < std::abs(min_weight) ) { 
      if( (int)((*rng)() + std::abs(it->weight)/std::abs(reset_weight)) == 0  ) {
         it=walkers.erase(it);
      } else {
        (it++)->weight = reset_weight;
      } 
    } else {
      ++it;
    }

  std::vector<Walker> w_;

  // now handle large weights
  for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)  
    if( it->alive && std::abs(it->weight) > std::abs(max_weight) ) {
  
      RealType w = std::abs(it->weight);
      int n = (int) (w/std::abs(reset_weight));
      RealType rem = w-n*std::abs(reset_weight);
      if( ( (int)((*rng)() + std::abs(rem/reset_weight) ) ) != 0 ) n++;      
      it->weight *= reset_weight/w; 
      for(int i=0; i<n-1; i++) {
        w_.push_back(SlaterDetWalker(*it));
      }

    } 

  int sz = walkers.size();
  walkers.resize(sz+w_.size());

  for(int i=0; i<w_.size(); i++)
    walkers[sz+i] = w_[i];

  // call load balance at the end of every population control event
  loadBalance(comm);
}

}

