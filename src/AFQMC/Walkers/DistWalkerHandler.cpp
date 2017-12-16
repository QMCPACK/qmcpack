
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
#include"AFQMC/Walkers/WalkerControl.hpp"
#include"AFQMC/Walkers/WalkerUtilities.hpp"


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

bool DistWalkerHandler::dumpSamplesHDF5(hdf_archive& dump, int nW_to_file) 
{
  if(nW_to_file==0) return true;
  if(head) { 

    int nW = numWalkers();

    push_walkers_to_front();

    std::vector<int> from(nproc_heads);
    MPI_Allgather(&nW,1,MPI_INT,from.data(),1,MPI_INT,MPI_COMM_TG_LOCAL_HEADS);

    int nWtot = std::accumulate(from.begin(),from.end(),int(0));
    int w0 = std::accumulate(from.begin(),from.begin()+rank_heads,int(0));

    if(nW_to_file<0) nW_to_file=nWtot;
    nW_to_file = std::min(nW_to_file,nWtot);

    // careful here, avoid sending extra information (e.g. B mats for back propg)
    int wlk_nterms = (9 + nrow*ncol);
    int wlk_sz = wlk_nterms*sizeof(ComplexType);
    int nwlk_per_block = std::min(std::max(1,hdf_block_size/wlk_sz),nW_to_file);
    int nblks = nW_to_file/nwlk_per_block + ((nW_to_file%nwlk_per_block>0)?1:0);
    std::vector<int> wlk_per_blk;
  
    if(rank_heads==0) {

      // check that restart data doesnot exist 
      counts.resize(nproc_heads);
      displ.resize(nproc_heads);
      wlk_per_blk.reserve(nblks);

      commBuff.reserve(wlk_nterms*nwlk_per_block);

      std::vector<int> Idata(6);
      Idata[0]=nW_to_file;
      Idata[1]=nblks;
      Idata[2]=wlk_nterms;
      Idata[3]=wlk_sz;
      Idata[4]=nrow;
      Idata[5]=ncol;

      dump.push("Walkers");
      dump.write(Idata,"dims");
      
    }

    // to extract a single walker from walker storage 
    MPI_Datatype rtype, stype;
    {
      MPI_Datatype stype_;
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &stype_);
      MPI_Type_commit(&stype_);
      MPI_Aint loc0, loc1;
      MPI_Get_address(walkers.values(), &loc0);
      MPI_Get_address(walkers.values()+walker_size, &loc1);
      MPI_Aint dis = loc1-loc0;
      MPI_Type_create_resized(stype_,0,dis,&stype); 
      MPI_Type_commit(&stype);
      MPI_Type_free(&stype_);

      // to deposit a single walker on contiguous memory
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &rtype);
      MPI_Type_commit(&rtype);
    }

    int nsent=0;
    // ready to send walkers to head in blocks
    for(int i=0, ndone=0; i<nblks; i++, ndone+=nwlk_per_block) {

      int nwlk_tot = std::min(nwlk_per_block,nW_to_file-ndone);  
      int nw_to_send=0; 
      if( w0+nsent >= ndone && w0+nsent < ndone + nwlk_tot) 
        nw_to_send = std::min(nW-nsent,(ndone+nwlk_tot)-(w0+nsent)); 
 
      if(rank_heads==0) {

        for(int p=0, nt=0; p<nproc_heads; p++) {
   
          int n_ = 0;
          int nn = nt + nW;   
          if( ndone+nwlk_tot > nt && ndone < nt+nW ) {
            if(ndone <= nt)
              n_ = std::min(nW,(ndone+nwlk_tot)-nt);    
            else    
              n_ = std::min(nt+nW-ndone,nwlk_tot);    
          }    

          counts[p]=n_;
          nt+=from[p];

        }    
        displ[0]=0;
        for(int p=1, nt=0; p<nproc_heads; p++) {
          nt += counts[p-1];
          displ[p]=nt;
        }

        commBuff.resize(wlk_nterms*nwlk_tot);
      }  

      // assumes current storage structure
      MPI_Gatherv( getSM(nsent), nw_to_send, stype,
                    commBuff.data(), counts.data(), displ.data(), rtype, 
                    0, MPI_COMM_TG_LOCAL_HEADS);  
      nsent += nw_to_send;
      
      if(rank_heads==0) {  
        dump.write(commBuff,std::string("walkers_")+std::to_string(i));
        wlk_per_blk.push_back(nwlk_tot);
      }  
  
    }

    if(rank_heads==0) {
      dump.write(wlk_per_blk,"wlk_per_blk");
      dump.pop();
    }

    MPI_Type_free(&rtype);
    MPI_Type_free(&stype);

  }

  myComm->barrier();

  return true;
}

// not thinking about efficiency since this is only used once during initialization 
bool DistWalkerHandler::restartFromHDF5(int nW_per_tg, hdf_archive& read, const std::string& tag, bool set_to_target)
{
  targetN_per_TG = nW_per_tg;
  // if working with a fixed number of walkers, set extra_empty_spaces to 2*targetN_per_TG 
  if(pop_control_type.find("simple")==std::string::npos)
    extra_empty_spaces = nW_per_tg;

  std::vector<int> Idata(6);
  if(head && rank_heads==0) {
    std::string path = "/Walkers/DistWalkerHandler";
    if(tag != std::string("")) path += std::string("/")+tag;
    if(!read.is_group( path )) {
      app_error()<<" ERROR: H5Group  could not find /Walkers/DistWalkerHandler/{tag} group in file. No restart data for walkers. \n";
      return false;
    }

    if(!read.push("Walkers")) return false;
    if(!read.push("DistWalkerHandler")) return false;
    if(tag != std::string("")) if(!read.push(tag)) return false;

    if(!read.read(Idata,"dims")) return false;
  }

  MPI_Bcast(Idata.data(),6,MPI_INT,0,myComm->getMPI());

  int nWtot = Idata[0];
  int wlk_nterms = (9 + nrow*ncol);
  int wlk_sz = wlk_nterms*sizeof(ComplexType);
  std::vector<int> wlk_per_blk;

  int nW_needed = set_to_target?nW_per_tg*nproc_heads:nWtot;

  int n0, nn;
  std::tie(n0,nn) = FairDivideBoundary(rank_heads,std::min(nW_needed,nWtot),nproc_heads);  
  
  int nw_local = set_to_target?nW_per_tg:(nn-n0); 

  if(nw_local+extra_empty_spaces > walkers.size()/walker_size)
    walkers.resize((nw_local+extra_empty_spaces)*walker_size);
  maximum_num_walkers=walkers.size()/walker_size;
  tot_num_walkers=0;

  // single reader for now
  if(head && rank_heads==0) {

    std::vector<int> from(nproc_heads);

    int nblks = Idata[1];
    if(Idata[4] != nrow || Idata[5] != ncol) {
      app_error()<<" Error reading walker restart file: Walker dimensions do not agree: " 
                 <<Idata[4] <<" "
                 <<Idata[5] <<" "
                 <<nrow <<" "
                 <<ncol <<std::endl;
    }    
    if(Idata[2] != wlk_nterms) {
      app_error()<<" ERROR: Size of walker data is not consistent in hdf5 file. \n";
      app_error()<<" sz_sim, sz_file: " <<wlk_nterms <<" " <<Idata[2] <<std::endl; 
      return false;
    } 
    if(Idata[3] != wlk_sz) {
      app_error()<<" ERROR: Memory usage of walker data is not consistent in hdf5 file. \n";
      app_error()<<" sz_sim, sz_file: " <<wlk_sz <<" " <<Idata[3] <<std::endl;
      return false;
    }

    app_log()<<"Found " <<nWtot <<" walkers on restart file." <<std::endl; 

    wlk_per_blk.resize(nblks);
    if(!read.read(wlk_per_blk,"wlk_per_blk")) return false;
    MPI_Bcast(wlk_per_blk.data(),wlk_per_blk.size(),MPI_INT,0,MPI_COMM_TG_LOCAL_HEADS);

    int maxnW = *std::max_element(wlk_per_blk.begin(),wlk_per_blk.end());
    commBuff.reserve(maxnW*wlk_nterms);

    // 1. read/bcast blocks of walkers
    // 2. fair share of walkers among procs
    // 3. if too many on file, stop reading
    // 4. if too few, replicate walkers at the end
    int nread=0;
    for(int i=0, nt=0; i<nblks; i++) {

      commBuff.resize(wlk_per_blk[i]*wlk_nterms); 
      if(!read.read(commBuff,std::string("walkers_")+std::to_string(i))) return false;

      MPI_Bcast(commBuff.data(),commBuff.size()*2,MPI_DOUBLE,0,MPI_COMM_TG_LOCAL_HEADS);

      if( nt+wlk_per_blk[i] > n0 && nn >= nt) {

        int from_i=0, ntake=0;
        if(nt <= n0) {
          if(nread != 0) {
            app_error()<<" Error in walker restart algorithm. nread!=0. " <<std::endl;
            return false;
          }    
          from_i = n0-nt;
          ntake = std::min(nn-n0,nt+wlk_per_blk[i]-n0);  
        } else {
          if(nread == 0) {
            app_error()<<" Error in walker restart algorithm. nread==0. " <<std::endl;
            return false;
          }    
          from_i = 0;
          ntake = std::min(nn-nt,wlk_per_blk[i]);  
        }
        nread += ntake;

        for(int k=0; k<ntake; k++) {
          std::copy(commBuff.data()+(from_i+k)*wlk_nterms,
                    commBuff.data()+(from_i+k+1)*wlk_nterms,
                    getSM(tot_num_walkers));
          walkers[walker_size*tot_num_walkers+data_displ[INFO]] = ComplexType(1.0);   
          tot_num_walkers++;  
        }

      } 

      nt += wlk_per_blk[i];
      if( nt >= nW_needed )
        break;
    }

    if(tot_num_walkers < nw_local) {
      int n0 = tot_num_walkers; 
      for(int to=n0+1, from=0; to<nw_local; to++,from++) {
          std::copy(getSM(from),
                    getSM(from)+walker_size, 
                    getSM(to));
          walkers[walker_size*to+data_displ[INFO]] = ComplexType(1.0);
          tot_num_walkers++;        
      }    
    }

    if(tag != std::string("")) read.pop();
    read.pop();
    read.pop();

  } else if(head) {

    std::vector<int> from(nproc_heads);

    int nblks = Idata[1];

    wlk_per_blk.resize(nblks);
    MPI_Bcast(wlk_per_blk.data(),wlk_per_blk.size(),MPI_INT,0,MPI_COMM_TG_LOCAL_HEADS);

    int maxnW = *std::max_element(wlk_per_blk.begin(),wlk_per_blk.end());
    commBuff.reserve(maxnW*wlk_nterms);

    // 1. read/bcast blocks of walkers
    // 2. fair share of walkers among procs
    // 3. if too many on file, stop reading
    // 4. if too few, replicate walkers at the end
    int nread=0;
    int n0, nn;
    std::tie(n0,nn) = FairDivideBoundary(rank_heads,std::min(nW_needed,nWtot),nproc_heads);  
    for(int i=0, nt=0; i<nblks; i++) {

      commBuff.resize(wlk_per_blk[i]*wlk_nterms); 
      MPI_Bcast(commBuff.data(),commBuff.size()*2,MPI_DOUBLE,0,MPI_COMM_TG_LOCAL_HEADS);

      if( nt+wlk_per_blk[i] > n0 && nn >= nt) {

        int from_i=0, ntake=0;
        if(nt <= n0) {
          if(nread != 0) {
            app_error()<<" Error in walker restart algorithm. nread!=0. " <<std::endl;
            return false;
          }    
          from_i = n0-nt;
          ntake = std::min(nn-n0,nt+wlk_per_blk[i]-n0);  
        } else {
          if(nread == 0) {
            app_error()<<" Error in walker restart algorithm. nread==0. " <<std::endl;
            return false;
          }    
          from_i = 0;
          ntake = std::min(nn-nt,wlk_per_blk[i]);  
        }
        nread += ntake;

        for(int k=0; k<ntake; k++) {
          std::copy(commBuff.data()+(from_i+k)*wlk_nterms,
                    commBuff.data()+(from_i+k+1)*wlk_nterms,
                    getSM(tot_num_walkers));
          walkers[walker_size*tot_num_walkers+data_displ[INFO]] = ComplexType(1.0);   
          tot_num_walkers++;  
        }

      } 

      nt += wlk_per_blk[i];
      if( nt >= nW_needed )
        break;
    }

    if(tot_num_walkers < nw_local) {
      int n0 = tot_num_walkers; 
      for(int to=n0+1, from=0; to<nw_local; to++,from++) {
          std::copy(getSM(from),
                    getSM(from)+walker_size, 
                    getSM(to));
          walkers[walker_size*to+data_displ[INFO]] = ComplexType(1.0);
          tot_num_walkers++;        
      }    
    }
  }

  myComm->barrier();
  reset_walker_count();  
  targetN = GlobalPopulation();

  return true;
}

bool DistWalkerHandler::dumpToHDF5(hdf_archive& dump, const std::string& tag)
{

  if(head) { 

    int nW = numWalkers();

    push_walkers_to_front();

    std::vector<int> from(nproc_heads);
    MPI_Allgather(&nW,1,MPI_INT,from.data(),1,MPI_INT,MPI_COMM_TG_LOCAL_HEADS);

    int nWtot = std::accumulate(from.begin(),from.end(),int(0));
    int w0 = std::accumulate(from.begin(),from.begin()+rank_heads,int(0));

    // careful here, avoid sending extra information (e.g. B mats for back propg)
    int wlk_nterms = (9 + nrow*ncol);
    int wlk_sz = wlk_nterms*sizeof(ComplexType);
    int nwlk_per_block = std::min(std::max(1,hdf_block_size/wlk_sz),nWtot);
    int nblks = nWtot/nwlk_per_block + ((nWtot%nwlk_per_block>0)?1:0);
    std::vector<int> wlk_per_blk;
  
    if(rank_heads==0) {

      // check that restart data doesnot exist 
      std::string path = "/Walkers/DistWalkerHandler";
      if(tag != std::string("")) path += std::string("/")+tag;
      if(dump.is_group( path )) {
        app_error()<<" ERROR: H5Group /Walkers/DistWalkerHandler/{tag} already exists in restart file. This is a bug and should not happen. Contact a developer.\n";
        return false;
      }

      counts.resize(nproc_heads);
      displ.resize(nproc_heads);
      wlk_per_blk.reserve(nblks);

      commBuff.reserve(wlk_nterms*nwlk_per_block);

      std::vector<int> Idata(6);
      Idata[0]=nWtot;
      Idata[1]=nblks;
      Idata[2]=wlk_nterms;
      Idata[3]=wlk_sz;
      Idata[4]=nrow;
      Idata[5]=ncol;

      dump.push("Walkers");
      dump.push("DistWalkerHandler");
      if(tag != std::string("")) dump.push(tag);       
      dump.write(Idata,"dims");
      
    }

    // to extract a single walker from walker storage 
    MPI_Datatype rtype, stype;
    {
      MPI_Datatype stype_;
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &stype_);
      MPI_Type_commit(&stype_);
      MPI_Aint loc0, loc1;
      MPI_Get_address(walkers.values(), &loc0);
      MPI_Get_address(walkers.values()+walker_size, &loc1);
      MPI_Aint dis = loc1-loc0;
      MPI_Type_create_resized(stype_,0,dis,&stype); 
      MPI_Type_commit(&stype);
      MPI_Type_free(&stype_);

      // to deposit a single walker on contiguous memory
      MPI_Type_contiguous(2*wlk_nterms, MPI_DOUBLE, &rtype);
      MPI_Type_commit(&rtype);
    }

    int nsent=0;
    // ready to send walkers to head in blocks
    for(int i=0, ndone=0; i<nblks; i++, ndone+=nwlk_per_block) {

      int nwlk_tot = std::min(nwlk_per_block,nWtot-ndone);  
      int nw_to_send=0; 
      if( w0+nsent >= ndone && w0+nsent < ndone + nwlk_tot) 
        nw_to_send = std::min(nW-nsent,(ndone+nwlk_tot)-(w0+nsent)); 
 
      if(rank_heads==0) {

        for(int p=0, nt=0; p<nproc_heads; p++) {
   
          int n_ = 0;
          int nn = nt + nW;   
          if( ndone+nwlk_tot > nt && ndone < nt+nW ) {
            if(ndone <= nt)
              n_ = std::min(nW,(ndone+nwlk_tot)-nt);    
            else    
              n_ = std::min(nt+nW-ndone,nwlk_tot);    
          }    

          counts[p]=n_;
          nt+=from[p];

        }    
        displ[0]=0;
        for(int p=1, nt=0; p<nproc_heads; p++) {
          nt += counts[p-1];
          displ[p]=nt;
        }

        commBuff.resize(wlk_nterms*nwlk_tot);
      }  

      // assumes current storage structure
      MPI_Gatherv( getSM(nsent), nw_to_send, stype,
                    commBuff.data(), counts.data(), displ.data(), rtype, 
                    0, MPI_COMM_TG_LOCAL_HEADS);  
      nsent += nw_to_send;

      if(rank_heads==0) {  
        dump.write(commBuff,std::string("walkers_")+std::to_string(i));
        wlk_per_blk.push_back(nwlk_tot);
      }  

      // not sure if necessary, but avoids avalanche of messages on head node
      MPI_Barrier(MPI_COMM_TG_LOCAL_HEADS);  
  
    }

    if(rank_heads==0) {
      dump.write(wlk_per_blk,"wlk_per_blk");
      if(tag != std::string("")) dump.pop();
      dump.pop();
      dump.pop();
    }

    MPI_Type_free(&rtype);
    MPI_Type_free(&stype);

  }

  myComm->barrier();
  
  return true;
}

bool DistWalkerHandler::parse(xmlNodePtr cur)
{

    if(cur == NULL)
      return false;

    app_log()<<"\n\n --------------- Parsing DistWalkerHandler input ------------------ \n\n";

    xmlNodePtr curRoot=cur;
    OhmmsAttributeSet oAttrib;
    oAttrib.add(name,"name");
    oAttrib.put(cur);
    walkerType = "collinear";

    load_balance_alg = "async";
    pop_control_type = "pair";

    ParameterSet m_param;
    m_param.add(reset_weight,"reset_weight","double");
    m_param.add(max_weight,"max_weight","double");
    m_param.add(min_weight,"min_weight","double");
    m_param.add(extra_empty_spaces,"extra_spaces","int");
    m_param.add(walkerType,"walker_type","std::string");
    m_param.add(load_balance_alg,"load_balance","std::string");
    m_param.add(pop_control_type,"pop_control","std::string");
    m_param.put(cur);

    std::for_each(load_balance_alg.begin(), load_balance_alg.end(), [](char & c){
        c = ::tolower(c);
    });
    if(load_balance_alg.find("all")!=std::string::npos) {
      app_log()<<" Using all-to-all (gather/scatter) load balancing algorithm.  \n";
      app_log()<<"    This algorithm is slow, does not scale, and is not recommended. Use async if possible. \n";  
    } else if(load_balance_alg.find("simple")!=std::string::npos) {
      app_log()<<" Using blocking (1-1) swap load balancing algorithm. " <<"\n";
    } else if(load_balance_alg.find("async")!=std::string::npos) {
      app_log()<<" Using asynchronous non-blocking swap load balancing algorithm. " <<"\n";
    } else {
      std::cerr<<" Error: Unknown load balancing algorithm: " <<load_balance_alg <<" \n";
      return false;  
    }  

    std::for_each(pop_control_type.begin(), pop_control_type.end(), [](char & c){
        c = ::tolower(c);
    });
    if(pop_control_type.find("mod_simple")!=std::string::npos) {
      app_log()<<" Using simple population control algorithm with fluctuating polulation size. \n";
      app_log()<<" Measurements are performed after population control is applied. \n";
    } else if(pop_control_type.find("simple")!=std::string::npos) {
      app_log()<<" Using simple population control algorithm with fluctuating polulation size. \n";
      app_log()<<" Measurements are performed before population control is applied. \n";
    } else if(pop_control_type.find("pair")!=std::string::npos) {
      app_log()<<" Using population control algorithm based on paired walker branching ( a la QWalk). \n"; 
    } else if(pop_control_type.find("comb")!=std::string::npos) {
      app_log()<<" Using population control algorithm based on comb method (See Booth, Gubernatis, PRE 2009). \n"; 
    } else if(pop_control_type.find("min")!=std::string::npos) {
      app_log()<<" Using population control algorithm based on minimum reconfiguration (Caffarel et al., 2000). \n"; 
    } else {
      std::cerr<<" Error: Unknown population control algorithm: " <<pop_control_type <<"\n"; 
      return false;
    }

    popcontrol_on_master = true;
    if(pop_control_type.find("dist")!=std::string::npos) popcontrol_on_master = false;

    std::transform(walkerType.begin(),walkerType.end(),walkerType.begin(),(int (*)(int)) tolower);
    cur = curRoot->children;
    while (cur != NULL) {
      std::string cname((const char*)(cur->name));
      if(cname =="something") {
      }
      cur = cur->next;
    }
    app_log()<<std::endl;
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
  int dummy=0;
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
  maximum_num_walkers=tot_num_walkers=0;
}

// called at the beginning of each executable section
void DistWalkerHandler::resetNumberWalkers(int n, bool a, ComplexMatrix* S) 
{

  if(head) {

    ComplexMatrix* S0 = S;
    if(S0==NULL) S0 = &HFMat;
    // assuming everything is alive
    int ns = maximum_num_walkers; //size(); 
    if(n+extra_empty_spaces > ns) {
      std::vector<ComplexType> tmp(tot_num_walkers*walker_size); // store old walker info
      int cnt=0;
      for(int i=0; i<maximum_num_walkers; i++)
        if(walkers[walker_size*i+data_displ[INFO]].real() > 0) 
          std::copy(walkers.begin()+walker_size*i,walkers.begin()+walker_size*(i+1),tmp.begin()+(cnt++)*walker_size);
      walkers.resize((n+extra_empty_spaces)*walker_size);
      maximum_num_walkers = n+extra_empty_spaces;
      tot_num_walkers = std::min(n,cnt); 
      std::copy(tmp.begin(),tmp.begin()+walker_size*tot_num_walkers,walkers.begin());
      for(int i=tot_num_walkers; i<maximum_num_walkers; i++)
        walkers[walker_size*i+data_displ[INFO]] = ComplexType(-1.0);   
    }


    // now I have enough space, add or remove walkers if needed 
    if(tot_num_walkers < n) {
      // adding walkers
      int cnt=tot_num_walkers;  
      for(int i=0; i<maximum_num_walkers; i++) {
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
      for(int i=maximum_num_walkers-1; i>=0; i--) {
        if(walkers[walker_size*i+data_displ[INFO]].real() > 0) {  
          walkers[walker_size*i+data_displ[INFO]] = ComplexType(-1.0);
          cnt--;
        }
        if(cnt == n) break;
      }
      if(cnt != n) APP_ABORT("Error: Problems in DistWalkerHandler::resetNumberWalkers-remove. \n");
    }


  } else {
    int ns = maximum_num_walkers; //size(); 
    if(n+extra_empty_spaces > ns) { 
      walkers.resize((n+extra_empty_spaces)*walker_size);
      maximum_num_walkers = n+extra_empty_spaces;
    }
  }
  walkers.barrier(); 
  tot_num_walkers = n; 
  nwalk_min = nwalk_max= tot_num_walkers;
  if(pop_control_type.find("pair")!=std::string::npos)
    nwalk_min = nwalk_max= 0; // in this case, this counts number of walkers left unpaired
}

// remove n walkers from the set and return their data
void DistWalkerHandler::pop_walkers(int n, ComplexType* data) {

  if(!head) {
    tot_num_walkers -= n;
    return;
  }  
  assert(tot_num_walkers >= n);
  ComplexSMVector::iterator wlk = walkers.end()-walker_size;

  int cnt = 0;
  while( cnt < n ) {

    while( (wlk+data_displ[INFO])->real() < 0 && wlk > walkers.begin() ) 
      wlk -= walker_size;
    if(wlk == walkers.begin() && ( (wlk+data_displ[INFO])->real() < 0 || (cnt+1 < n)  ))
      APP_ABORT(" Error in DistWalkerHandler::pop_walkers. Could not find walkers. \n\n\n"); 

    std::copy(wlk,wlk+walker_size,data); 
    *(wlk+data_displ[INFO]) = ComplexType(-1.0,0.0); 
    data += walker_size;
    tot_num_walkers--;
    cnt++; 

  } 

}

// add n walkers to the set. There must be enough space since this is called in serial
void DistWalkerHandler::push_walkers(int n, ComplexType* data) {

  if(!head) {
    tot_num_walkers += n;
    return;
  }  
  assert(tot_num_walkers+n <= maximum_num_walkers);
  ComplexSMVector::iterator wlk = walkers.begin();

  int cnt = 0;
  while( cnt < n ) {

    while( (wlk+data_displ[INFO])->real() > 0 && wlk != walkers.end() ) 
      wlk += walker_size;
    if(wlk == walkers.end())
      APP_ABORT(" Error in DistWalkerHandler::push_walkers. Could not find space. \n\n\n"); 

    std::copy(data,data+walker_size,wlk);
    data += walker_size;
    tot_num_walkers++;
    cnt++; 

  } 

}

// load balancing algorithm
//  CAREFUL!!!! This assumes a communicator of heads_of_TG!!!!
//  Will give unexpected results for anything else!!!
void DistWalkerHandler::loadBalance(MPI_Comm comm) 
{
  
#if defined(HAVE_MPI)
  nwalk_min = nwalk_max= tot_num_walkers;
  if(nproc_heads==1) return;

  int nW = numWalkers();
  int sz = size(), nw_new=0;
  MPI_Request request;

  int npr,comm_rank; 
  if(head) {
    MPI_Comm_size(comm,&npr);
    MPI_Comm_rank(comm,&comm_rank);
    assert(nwalk_counts_new.size() == npr);
    assert(nwalk_counts_old.size() == npr);
    nw_new = nwalk_counts_new[comm_rank];
  }

  Timers[LoadBalance_Setup]->start();
  LocalTimer->start("WalkerHandler::loadBalance::resize");
  // resize arrays if necessary. This requires all cores in a TG 
  walkers.share(&nw_new,1,head);
  if(nw_new > sz) {
    walkers.resize(walker_size*(nw_new+extra_empty_spaces));
    maximum_num_walkers = nw_new+extra_empty_spaces; 
  }
  LocalTimer->stop("WalkerHandler::loadBalance::resize");
  Timers[LoadBalance_Setup]->stop();

// trying to understand the actual cost of exchanging walkers,
// remove this when algorithm works as expected
myComm->barrier();

  Timers[LoadBalance_Exchange]->start();
  LocalTimer->start("WalkerHandler::loadBalance::exchange");
  if(load_balance_alg.find("all")!=std::string::npos) {

    if(head) {

      int ncomm = nW - nwalk_counts_new[comm_rank];
      int pos = 0, cnt=0;
      ComplexType *buffer=NULL; 
      if(ncomm > 0) 
        buffer = walkers.values() + walker_size*nwalk_counts_new[comm_rank];
       else 
        buffer = walkers.values() + walker_size*nW;

      if(comm_rank==0) {

        counts.resize(npr);
        displ.resize(npr);

        // setup gatherv call 
        displ[0]=0;
        for(int i=0; i<npr-1; i++)
          if(nwalk_counts_old[i] > nwalk_counts_new[i]) {
            counts[i] = (nwalk_counts_old[i]-nwalk_counts_new[i])*walker_size;
            displ[i+1] = displ[i]+counts[i];
          } else {
            counts[i]=0;
            displ[i+1] = displ[i];
          }
        if(nwalk_counts_old[npr-1] > nwalk_counts_new[npr-1]) 
          counts[npr-1] = (nwalk_counts_old[npr-1]-nwalk_counts_new[npr-1])*walker_size;
        else
          counts[npr-1]=0;
        commBuff.clear();
        commBuff.resize(displ[npr-1]+counts[npr-1]);
        int nn = walker_size*((ncomm>0)?ncomm:0); 
        myComm->gatherv( buffer, commBuff.data(), nn , counts, displ, 0, comm); 

        // setup scatterv call
        displ[0]=0;
        for(int i=0; i<npr-1; i++)
          if(nwalk_counts_old[i] < nwalk_counts_new[i]) {
            counts[i] = (nwalk_counts_new[i]-nwalk_counts_old[i])*walker_size;
            displ[i+1] = displ[i]+counts[i];
          } else {
            counts[i]=0;
            displ[i+1] = displ[i];
          }
        if(nwalk_counts_old[npr-1] < nwalk_counts_new[npr-1])
          counts[npr-1] = (nwalk_counts_new[npr-1]-nwalk_counts_old[npr-1])*walker_size;
        else
          counts[npr-1]=0;    

        nn = walker_size*((ncomm<0)?(std::abs(ncomm)):0); 
        myComm->scatterv( commBuff.data(), buffer, nn, counts, displ, 0, comm); 

        if(ncomm > 0) {
          register ComplexType minus = ComplexType(-1.0,0.0);
          for(ComplexSMVector::iterator it=walkers.begin()+nwalk_counts_new[comm_rank]*walker_size; it<walkers.end(); it+=walker_size)
            *(it+data_displ[INFO]) = minus; 
        }

      } else {

        int nn = walker_size*((ncomm>0)?ncomm:0);
        myComm->gatherv( buffer, commBuff.data(), nn, counts, displ, 0, comm);
        nn = walker_size*((ncomm<0)?(std::abs(ncomm)):0);
        myComm->scatterv( commBuff.data(), buffer, nn, counts, displ, 0, comm); 

        if(ncomm > 0) {
          register ComplexType minus = ComplexType(-1.0,0.0);
          for(ComplexSMVector::iterator it=walkers.begin()+nwalk_counts_new[comm_rank]*walker_size; it<walkers.end(); it+=walker_size)
            *(it+data_displ[INFO]) = minus; 
        }

      }
    } // head 

  } else if(load_balance_alg.find("simple")!=std::string::npos) {

    if(head) 
      afqmc::swapWalkersSimple(*this,nwalk_counts_old,nwalk_counts_new,comm);
    
  } else if(load_balance_alg.find("async")!=std::string::npos) {

    if(head) 
      afqmc::swapWalkersAsync(*this,nwalk_counts_old,nwalk_counts_new,comm);

  }

// trying to understand the actual cost of exchanging walkers,
// remove this when algorithm works as expected
myComm->barrier();

  LocalTimer->stop("WalkerHandler::loadBalance::exchange");
  Timers[LoadBalance_Exchange]->stop();

  walkers.barrier();
  reset_walker_count();

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

  while( wlk > empty ) {

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

// performs pair branches in current list of walkers. The location of
// unpaired walkers are returned in windx and the return value defines
// the type of walker in windx (<0 small weight, >0 large weight)
int DistWalkerHandler::local_pair_branching(std::vector<int>& windx)
{
  windx.clear();
  push_walkers_to_front();
  reset_walker_count();
  int type = 0;

  int np=-1;
  int nm=-1;
  while(np < tot_num_walkers && nm < tot_num_walkers) {

    if(np < tot_num_walkers) {
      while(++np < tot_num_walkers) {
        if( std::abs(getWeight(np)) > std::abs(max_weight))
          break;
      }
    }
    if(nm < tot_num_walkers) {
      while(++nm < tot_num_walkers) {
        if( std::abs(getWeight(nm)) < std::abs(min_weight))
          break;
      }
    }

    if( np < tot_num_walkers ) {
    
      if( nm < tot_num_walkers ) {

        double w12 = getWeight(np).real() + getWeight(nm).real();
        if((*rng)() < getWeight(np).real()/w12) 
          pair_branch(w12*0.5,np,nm);
        else 
          pair_branch(w12*0.5,nm,np);

      } else {
        windx.push_back(np);
        type=1;
      }

    } else {

      if( nm < tot_num_walkers ) {
        windx.push_back(nm);
        type=-1;
      }

    }

  }
  return type;

}

// copies walker from branch_from to branch_to and resets the weights of both walkers to w
void DistWalkerHandler::pair_branch(RealType w, int branch_from, int branch_to)
{
  if(!isAlive(branch_from)) {
    app_error()<<"  Error: Calling pair_branch on dead walker. \n\n\n";
    APP_ABORT("  Error: Calling pair_branch on dead walker. \n\n\n");
  }
  setWeight(branch_from,w);
  std::copy(walkers.begin()+branch_from*walker_size,  
           walkers.begin()+(branch_from+1)*walker_size, 
           walkers.begin()+branch_to*walker_size);
}

// sets the weight of walker "n" to "w" and adds "num" copies to the list
// list must have sufficient space to add the requested number of walkers
// list must be compact 
void DistWalkerHandler::branch(RealType w, int n, int num)
{
  assert(maximum_num_walkers >= tot_num_walkers+num); 
  setWeight(n,ComplexType(w,0.0));
  ComplexSMVector::iterator itn=walkers.begin()+walker_size*n;
  ComplexSMVector::iterator it=walkers.begin();
  assert( (itn+data_displ[INFO])->real() > 0 );
  for(int i=0; i<num; i++) {
    while( it < walkers.end() ) {
      if( (it+data_displ[INFO])->real() < 0 ) break;
      it+=walker_size;
    }
    assert( it < walkers.end() );
    std::copy(itn,itn+walker_size,it);
  }
  tot_num_walkers += num;
}

// population control algorithm
//  CAREFUL!!!! This assumes a communicator of heads_of_TG!!!!
//  Will give unexpected results for anything else!!!
//
//  curData:
//  0: factor used to rescale the weights
//  1: sum_i w_i * Eloc_i   (where w_i is the unnormalized weight)
//  2: sum_i w_i            (where w_i is the unnormalized weight)
//  3: sum_i abs(w_i)       (where w_i is the unnormalized weight)
//  4: sum_i abs(<psi_T|phi_i>)
//  5: total number of walkers  
//  6: total number of "healthy" walkers (those with weight > 1e-6, ovlp>1e-8, etc) 
void DistWalkerHandler::popControl(MPI_Comm comm, std::vector<ComplexType>& curData)
{
  ComplexType minus = ComplexType(-1.0,0.0);
  bool needsLoadBalance = true;

  // temporarily
  comm = MPI_COMM_TG_LOCAL_HEADS;

  curData.resize(7);
  std::fill(curData.begin(),curData.begin()+7,ComplexType(0));

// trying to understand time spend idle
  LocalTimer->start("WalkerHandler::popControl::idle");
  walkers.barrier();   // DO NOT REMOVE THIS BARRIER
  if(head)
   MPI_Barrier(comm);
  LocalTimer->stop("WalkerHandler::popControl::idle");

  // the "simple" algorithm is inherently local, so check regardless of popcontrol_on_master 
  LocalTimer->start("WalkerHandler::popControl");
  if(pop_control_type.find("simple")!=std::string::npos) {

    bool before = pop_control_type.find("mod_simple")==std::string::npos;
    // first rescale weights
    ComplexType wfact=ComplexType(1.0,0.0);
    if(head) {

      // calculate basic averages and rescale weights 
      afqmc::BasicWalkerData(*this,curData,MPI_COMM_TG_LOCAL_HEADS);
      RealType scl = 1.0/curData[0].real();
      scaleWeight(scl);
      // not sure why I'm doing this, doesn't change anything at the end
      curData[1] *= scl; 
      curData[2] *= scl;
      curData[3] *= scl;

      // branching step 
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
            if( (int)((*rng)() + std::abs(w0)/reset_weight) == 0 ) { 
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
          num_new += std::floor( std::abs( *(it+data_displ[WEIGHT]) ) / std::abs(reset_weight) )-1;  

      if( num_new > empty_spots.size() ) {
        int ntot = nw + num_new + extra_empty_spaces; 
        std::vector<ComplexType> tmp(nw*walker_size); // store old walker info
        cnt=0;
        for(ComplexSMVector::iterator it=walkers.begin(); it<walkers.end(); it+=walker_size) 
          if( (it+data_displ[INFO])->real() > 0 )
            std::copy( it, it+walker_size, tmp.begin()+(cnt++)*walker_size);
        assert(nw==cnt);
        int nn = ntot*walker_size;
        walkers.share(&ntot,1,head);
        walkers.resize(nn,false);
        maximum_num_walkers = ntot;
        std::copy(tmp.begin(),tmp.begin()+walker_size*nw,walkers.begin());
        for(ComplexSMVector::iterator it=walkers.begin()+nw*walker_size; it<walkers.end(); it+=walker_size) 
          *(it+data_displ[INFO]) = minus; 
 
        empty_spots.clear();
        empty_spots.reserve( maximum_num_walkers-nw );
        for(int i=nw; i<maximum_num_walkers; i++) empty_spots.push_back(i);
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
          RealType w0 = w/static_cast<RealType>(n);
          *(it+data_displ[WEIGHT]) *= w0/w; 
          for(int i=0; i<n-1; i++,cnt++) 
            std::copy( it, it+walker_size, walkers.begin()+walker_size*empty_spots[cnt] ); 
          // keep count of number of walkers
          nw += (n-1);
        }
  
      empty_spots.clear();

      tot_num_walkers = nw;
     
    } else {
      int ntot; 
      walkers.share<int>(&ntot,1,head);
      if(ntot > 0) {
        walkers.resize(ntot*walker_size,false);
        maximum_num_walkers = ntot;
      }
    }
    // end of branching

    if(!before && head) {
      RealType scl = curData[0].real();
      // recalculate curData 
      afqmc::BasicWalkerData(*this,curData,MPI_COMM_TG_LOCAL_HEADS);
      curData[0] = ComplexType(scl,0.0);
    } 

    // setup walker counts
    if(head) {
      int npr; 
      MPI_Comm_size(comm,&npr);
      nwalk_counts_new.resize(npr);
      nwalk_counts_old.resize(npr);

      LocalTimer->start("WalkerHandler::loadBalance::setup");
      Timers[LoadBalance_Setup]->start();
      // determine new number of walkers

      for(int i=0; i<npr; i++)
        nwalk_counts_old[i] = nwalk_counts_new[i] = 0;
      nwalk_counts_new[0] = tot_num_walkers;
      // in case iallgather is not available
#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER 
#define HAVE_MPI_IGATHER 
#define HAVE_MPI_IBCAST 
#endif
#ifdef HAVE_MPI_IALLGATHER
      MPI_Request request;
      MPI_Iallgather(nwalk_counts_new.data(), 1, MPI_INT, nwalk_counts_old.data(), 1, MPI_INT, comm,&request); 

      // push empty spots to the end of the list 
      push_walkers_to_front();

      // wait for mpi_iallgather 
      MPI_Wait(&request, MPI_STATUS_IGNORE); 
#else
      MPI_Allgather(nwalk_counts_new.data(), 1, MPI_INT, nwalk_counts_old.data(), 1, MPI_INT, comm);
      push_walkers_to_front();
#endif 
      nwalk_global=0;
      for(int i=0; i<npr; i++) 
        nwalk_global+=nwalk_counts_old[i];
      // Fair Share distribution
      for(int i=0; i<npr; i++) 
        nwalk_counts_new[i] = nwalk_global/npr + ((i<nwalk_global%npr)?(1):(0));
      auto min_max = std::minmax_element(nwalk_counts_old.begin(),nwalk_counts_old.end());
      nwalk_min = *min_max.first;
      nwalk_max = *min_max.second;

      LocalTimer->stop("WalkerHandler::loadBalance::setup");
      Timers[LoadBalance_Setup]->stop();
    }

  // are we doing popcontol on the master node???
  } else if(popcontrol_on_master) {

    // gather data and walker information
    if(head) {
      afqmc::BasicWalkerData(*this,curData,MPI_COMM_TG_LOCAL_HEADS);
      RealType scl = 1.0/curData[0].real();
      scaleWeight(scl);
      // not sure why I'm doing this, doesn't change anything at the end
      curData[1] *= scl; 
      curData[2] *= scl;
      curData[3] *= scl;
    } 

    int npr,rk;
    if(head) {
      MPI_Comm_size(comm,&npr);
      MPI_Comm_rank(comm,&rk);
      nwalk_counts_old.resize(npr);
      nwalk_counts_new.resize(npr);

      int ione = 1;
      assert(tot_num_walkers==targetN_per_TG);
      if(rk==0) 
        bufferall.resize(targetN*(sizeof(double)+sizeof(int)));  
      else
        bufferall.resize(1);
      bufferlocal.resize(targetN_per_TG*(sizeof(double)+sizeof(int)));  
      push_walkers_to_front();

      for(int i=0,sz=0; i<tot_num_walkers; i++) {
        double e = static_cast<double>(getWeight(i).real());
        std::memcpy(bufferlocal.data()+sz,&(e),sizeof(double));
        sz += sizeof(double);
        std::memcpy(bufferlocal.data()+sz,&(ione),sizeof(int));
        sz += sizeof(int);
      } 
      MPI_Gather(bufferlocal.data(),bufferlocal.size(),MPI_CHAR,
                 bufferall.data(),bufferlocal.size(),MPI_CHAR,
                  0,comm);

      if(rk==0) 
        if(pop_control_type.find("pair")!=std::string::npos) { 
          afqmc::SerialPairBranching(bufferall,*rng,max_weight,min_weight); 
        } else if(pop_control_type.find("comb")!=std::string::npos) { 
          APP_ABORT(" Error: Not implemented yet. \n\n\n");
        } else if(pop_control_type.find("comb")!=std::string::npos) { 
          APP_ABORT(" Error: Not implemented yet. \n\n\n");
        }

#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER 
#define HAVE_MPI_IGATHER 
#define HAVE_MPI_IBCAST 
#endif
#ifdef HAVE_MPI_IGATHER
      MPI_Request req1;
      MPI_Iscatter(bufferall.data(),bufferlocal.size(),MPI_CHAR,
                 bufferlocal.data(),bufferlocal.size(),MPI_CHAR,
                  0,comm,&req1);
#else
      MPI_Scatter(bufferall.data(),bufferlocal.size(),MPI_CHAR,
                 bufferlocal.data(),bufferlocal.size(),MPI_CHAR,
                  0,comm);
#endif

      if(rk==0) {
        int nbranch=0;
        for(int i=0, sz=0; i<npr; i++) {  
          int cnt=0;
          for(int k=0; k<targetN_per_TG; k++) { 
            int n;
            double e;
            std::memcpy(&e,bufferall.data()+sz,sizeof(double));
            sz += sizeof(double);
            std::memcpy(&(n),bufferall.data()+sz,sizeof(int));
            sz += sizeof(int);
            cnt += n;
            if(n==0) nbranch++; 
          }
          nwalk_counts_old[i]=cnt;
          max_nexch = std::max(max_nexch,std::abs(targetN_per_TG-cnt)); 
        }
        max_nbranch = std::max(max_nbranch,nbranch);
        if(targetN!=std::accumulate(nwalk_counts_old.begin(),nwalk_counts_old.end(),0)) {
          app_log()<<" Error: targetN != nwold: " <<targetN <<" " 
                   <<std::accumulate(nwalk_counts_old.begin(),nwalk_counts_old.end(),0) 
                   <<std::endl;
        } 
        assert(targetN==std::accumulate(nwalk_counts_old.begin(),nwalk_counts_old.end(),0));
      }
      std::fill(nwalk_counts_new.begin(),nwalk_counts_new.end(),targetN_per_TG);

#ifdef HAVE_MPI_IBCAST
      MPI_Request req2;
      MPI_Ibcast(nwalk_counts_old.data(),nwalk_counts_old.size(),MPI_INT,0,comm,&req2);
#else
      MPI_Bcast(nwalk_counts_old.data(),nwalk_counts_old.size(),MPI_INT,0,comm);
#endif

      // birth/death step 
#ifdef HAVE_MPI_IGATHER
      MPI_Status st1;
      MPI_Wait(&req1,&st1);
#endif      

      // missing!!! need to count the number of new walkers and make sure they
      // fit, otherwise need to resize SM array in coordination with TG
      // Easy to do, but finish later!!!
      int nnew = 0;
      int nempty = size() - tot_num_walkers;
      for(int i=0,sz=0,ni=0; i<targetN_per_TG; i++) {
        sz += sizeof(double);
        std::memcpy(&ni,bufferlocal.data()+sz,sizeof(int));
        sz += sizeof(int);
        assert(ni >= 0);
        if(ni > 1) nnew+=ni-1; 
      }
      assert(nempty >= nnew);

      // from now until load balance, targetN_per_TG != tot_num_walkers
      double e;
      for(int i=0,sz=0,ni=0; i<targetN_per_TG; i++) {
        std::memcpy(&e,bufferlocal.data()+sz,sizeof(double));
        sz += sizeof(double);
        std::memcpy(&ni,bufferlocal.data()+sz,sizeof(int));
        sz += sizeof(int);
        if(ni == 0) {
          kill(i); 
        } else if(ni>1) {
          branch(e,i,ni-1);
        }
      }
      push_walkers_to_front();

#ifdef HAVE_MPI_IBCAST
      MPI_Status st2;
      MPI_Wait(&req2,&st2);
#endif

    } else {

      // if resizing is needed, to it here  
      assert(tot_num_walkers==targetN_per_TG);

    }

  // otherwise call distributed routine
  } else if(pop_control_type.find("pair")!=std::string::npos) { 

// right now all the work is done at the master node (myComm->rank() == 0), 
// minimize global communication later

    needsLoadBalance = false;

    // gather data and walker information
    if(head) {
      afqmc::BasicWalkerData(*this,curData,MPI_COMM_TG_LOCAL_HEADS);
      RealType scl = 1.0/curData[0].real();
      scaleWeight(scl);
      // not sure why I'm doing this, doesn't change anything at the end
      curData[1] *= scl; 
      curData[2] *= scl;
      curData[3] *= scl;
    } 
  
    // perform branching step
    // 1. local pair branching 
    if(head) {
      int nleft = afqmc::NaivePairBranching(*this,*rng,comm); 
      if(nleft < 0) nwalk_min = std::min(nwalk_min,nleft); 
      if(nleft > 0) nwalk_max = std::max(nwalk_max,nleft); 
    }

  } else if(pop_control_type.find("comb")!=std::string::npos) { 
    APP_ABORT(" Error: Not implemented yet. \n\n\n");
  } else if(pop_control_type.find("min")!=std::string::npos) { 
    APP_ABORT(" Error: Not implemented yet. \n\n\n");
  }
  LocalTimer->stop("WalkerHandler::popControl");

  walkers.barrier();
  reset_walker_count();

  LocalTimer->start("WalkerHandler::popControl::bcast");
#if MPI_VERSION >= 3
#define HAVE_MPI_IALLGATHER 
#define HAVE_MPI_IGATHER 
#define HAVE_MPI_IBCAST 
#endif
#ifdef HAVE_MPI_IBCAST
  MPI_Request req3;
  MPI_Ibcast(curData.data(),curData.size(),MPI_INT,0,MPI_COMM_TG_LOCAL,&req3);
#else
  MPI_Bcast(curData.data(),curData.size(),MPI_INT,0,MPI_COMM_TG_LOCAL);
#endif
  LocalTimer->stop("WalkerHandler::popControl::bcast");

  // load balance after population control events
  if(needsLoadBalance) loadBalance(comm);

  if(pop_control_type.find("simple")==std::string::npos) {
    assert(tot_num_walkers == targetN_per_TG);
  }

  LocalTimer->start("WalkerHandler::popControl::bcast");
#ifdef HAVE_MPI_IBCAST
  MPI_Status st3;
  MPI_Wait(&req3,&st3);
#endif 
  LocalTimer->stop("WalkerHandler::popControl::bcast");

}

void DistWalkerHandler::benchmark(std::string& blist,int maxnW,int delnW,int repeat)
{

  if(blist.find("comm")!=std::string::npos) {

    app_log()<<" Testing communication times in WalkerHandler. This should be done using a single TG per node, to avoid timing communication between cores on the same node. \n";
    std::ofstream out;
    if(myComm->rank() == 0)
      out.open(myComm->getName()+".icomm.dat");  

    MPI_Comm_size(MPI_COMM_TG_LOCAL_HEADS,&nproc_heads);
    MPI_Comm_rank(MPI_COMM_TG_LOCAL_HEADS,&rank_heads);
    assert(nproc_heads > 1);

    std::vector<std::string> tags(3);
    tags[0]="M1";
    tags[1]="M2";
    tags[2]="M3";

    for( std::string& str: tags) Timer.reset(str);

    int nw=1;
    while(nw <= maxnW) {

      if(head && (rank_heads==0 || rank_heads==1)) {
        int sz = nw*walker_size; 
        std::vector<ComplexType> Cbuff(sz);
        MPI_Request req;
        MPI_Status st;    
        MPI_Barrier(MPI_COMM_TG_LOCAL_HEADS);
        for(int i=0; i<repeat; i++) {

          if(rank_heads==0) {
            Timer.start("M1");
            MPI_Isend(Cbuff.data(),2*Cbuff.size(),MPI_DOUBLE,1,999,MPI_COMM_TG_LOCAL_HEADS,&req);  
            MPI_Wait(&req,&st);
            Timer.stop("M1");
          } else {
            MPI_Irecv(Cbuff.data(),2*Cbuff.size(),MPI_DOUBLE,0,999,MPI_COMM_TG_LOCAL_HEADS,&req);  
            MPI_Wait(&req,&st);
          }

        }

        if(rank_heads == 0) {
          out<<nw <<" " ;
          for( std::string& str: tags) out<<Timer.total(str)/double(repeat) <<" ";
          out<<std::endl;
        }
      } else if(head) {
        MPI_Barrier(MPI_COMM_TG_LOCAL_HEADS);
      }

      if(delnW <=0) nw *= 2;
      else nw += delnW;
    }


  } else if(blist.find("comm")!=std::string::npos) {
    std::ofstream out;
    if(myComm->rank() == 0)
      out.open(myComm->getName()+".comm.dat"); 

  }
 
}

}

