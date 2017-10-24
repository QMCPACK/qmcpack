
// -*- C++ -*-
// /**@file WalkerHandler 
//  * @brief Virtual Class for walker handlers. 
//   */

#ifndef QMCPLUSPLUS_AFQMC_DISTWALKERHANDLER_H
#define QMCPLUSPLUS_AFQMC_DISTWALKERHANDLER_H

#include<random>

#include "OhmmsData/libxmldefs.h"
#include "io/hdf_archive.h"

#include"AFQMC/config.h"
#include"AFQMC/Walkers/WalkerHandlerBase.h"
#include<Message/MPIObjectBase.h>
#include "Message/CommOperators.h"
#include "Utilities/NewTimer.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"


namespace qmcplusplus
{

/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class DistWalkerHandler: public WalkerHandlerBase 
{

  public:

  /// constructor
  DistWalkerHandler(Communicate* c,  RandomGenerator_t* r): WalkerHandlerBase(c,r,std::string("DistWalkerHandler")),
      min_weight(0.05),max_weight(4.0)
      ,reset_weight(1.0),extra_empty_spaces(10)
      ,head(false),walker_memory_usage(0),tot_num_walkers(0),maximum_num_walkers(0)
      ,popcontrol_on_master(true)
  { }

  /// destructor
  ~DistWalkerHandler() {}

  inline int size() { 
    assert( maximum_num_walkers == walkers.size()/walker_size );
    return maximum_num_walkers; 
  } 

  bool restartFromXML(); 
  bool dumpToXML();

  bool restartFromHDF5(int n, hdf_archive&, const std::string&, bool set_to_target); 
  bool dumpToHDF5(hdf_archive&, const std::string&); 

  bool dumpSamplesHDF5(hdf_archive& dump, int nW); 

  // reads xml and performs setup
  bool parse(xmlNodePtr cur); 

  // performs setup
  bool setup(int,int,int,MPI_Comm,MPI_Comm,MPI_Comm,myTimer*);

  // cleans state of object. 
  //   -erases allocated memory 
  bool clean(); 

  // called at the beginning of each executable section
  void resetNumberWalkers(int n, bool a=true, ComplexMatrix* S=NULL); 

  // perform and report tests/timings
  void benchmark(std::string& blist,int maxnW,int delnW,int repeat);

  inline bool initWalkers(int n) {
   targetN_per_TG = n;
   // if working with a fixed number of walkers, set extra_empty_spaces to 2*targetN_per_TG 
   if(pop_control_type.find("simple")==std::string::npos) 
     extra_empty_spaces = n;
   resetNumberWalkers(n,true,&HFMat);
   targetN = GlobalPopulation();
   return true;   
  } 

  inline void set_TG_target_population(int N) { targetN_per_TG = N; }

  inline int get_TG_target_population() { return targetN_per_TG; }
  inline int get_global_target_population() { return targetN; }

  inline int numWalkers(bool dummy=false) {
// checking
    int cnt=0;
    for(ComplexSMVector::iterator it=walkers.begin()+data_displ[INFO]; 
       it<walkers.end(); it+=walker_size)
      if(it->real() > 0) cnt++; 
    if(cnt != tot_num_walkers)
      APP_ABORT(" Error in DistWalkerHandler::numWalkers(): Incorrect number of walkers. \n"); 

    if(dummy)
      return size();
    else
      return tot_num_walkers;
  } 

  inline int GlobalPopulation() {
    std::vector<int> res(1);
    res[0]=0;
    if(head)
      res[0] += tot_num_walkers;
    myComm->gsum(res);
    return res[0];
  }

  inline RealType GlobalWeight() {
    std::vector<RealType> res(1);
    res[0]=0;
    if(head)
      for(int i=0; i<maximum_num_walkers; i++)
       if(isAlive(i)) 
        res[0] += std::abs(getWeight(i)); 
    myComm->gsum(res);
    return res[0];
  }

  // load balancing algorithm
  void loadBalance(MPI_Comm comm); 

  // population control algorithm
  void popControl(MPI_Comm comm, std::vector<ComplexType>& curData); 

  void setHF(const ComplexMatrix& HF);

  inline void Orthogonalize(int i) {

    if(walkerType == "closed_shell")  {
      DenseMatrixOperators::GeneralizedGramSchmidt(&(walkers[walker_size*i+data_displ[SM]]),NAEA,NMO,NAEA);
    } else if(walkerType == "collinear") {
      DenseMatrixOperators::GeneralizedGramSchmidt(&(walkers[walker_size*i+data_displ[SM]]),NAEA,NMO,NAEA);
      DenseMatrixOperators::GeneralizedGramSchmidt(&((walkers)[walker_size*i+data_displ[SM]])+NAEA*NMO,NAEA,NMO,NAEB);
    } else if(walkerType == "non-collinear") {
      APP_ABORT("ERROR: non-collinear not implemented in Orthogonalize. \n\n\n");
    }
  }

  inline void Orthogonalize() {
    for(int i=0; i<maximum_num_walkers; i++)
      if( i%ncores_per_TG == core_rank && isAlive(i))
        Orthogonalize(i);
  } 

  //private:
 
  enum walker_data { INFO=0, SM=1, WEIGHT=2, ELOC=3, ELOC_OLD=4, OVLP_A=5, OVLP_B=6, OLD_OVLP_A=7, OLD_OVLP_B=8};

  // n is zero-based
  ComplexType* getSM(int n) { if(n>=maximum_num_walkers) {return NULL;} else {return &((walkers)[walker_size*n+data_displ[SM]]);} }
  ComplexType getWeight(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[WEIGHT]];} }
  ComplexType getEloc(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC]];} }
  ComplexType getOldEloc(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC_OLD]];} }
  ComplexType getOvlpAlpha(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_A]];} }
  ComplexType getOvlpBeta(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_B]];} }
  ComplexType getOldOvlpAlpha(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_A]];} }
  ComplexType getOldOvlpBeta(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_B]];} }
  void setWeight(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else {(walkers)[walker_size*n+data_displ[WEIGHT]]=Q;} }
  void setEloc(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC]] = Q;} }
  void setOldEloc(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC_OLD]] = Q;} }
  void setOvlp(int n, ComplexType Q1, ComplexType Q2) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[OVLP_A]] = Q1; 
      (walkers)[walker_size*n+data_displ[OVLP_B]] = Q2;
    } 
  }
  void setOldOvlp(int n, ComplexType Q1, ComplexType Q2) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]] = Q1; 
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]] = Q2;
    } 
  }
  void setCurrToOld(int n) { 
    if(n>=maximum_num_walkers) {return;} 
    else { 
      (walkers)[walker_size*n+data_displ[ELOC_OLD]] = (walkers)[walker_size*n+data_displ[ELOC]];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]] = (walkers)[walker_size*n+data_displ[OVLP_A]];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]] = (walkers)[walker_size*n+data_displ[OVLP_B]];
    } 
  }
  ComplexType getEloc2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC]+1];} }
  ComplexType getOldEloc2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[ELOC_OLD]+1];} }
  ComplexType getOvlpAlpha2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_A]+1];} }
  ComplexType getOvlpBeta2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OVLP_B]+1];} }
  ComplexType getOldOvlpAlpha2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1];} }
  ComplexType getOldOvlpBeta2(int n) { if(n>=maximum_num_walkers) {return zero;} else {return (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1];} }
  void setEloc2(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC]+1] = Q;} }
  void setOldEloc2(int n, ComplexType Q) { if(n>=maximum_num_walkers) {return;} else { (walkers)[walker_size*n+data_displ[ELOC_OLD]+1] = Q;} }
  void setOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[OVLP_A]+1] = Q1;
      (walkers)[walker_size*n+data_displ[OVLP_B]+1] = Q2;
    }
  }
  void setOldOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1] = Q1;
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1] = Q2;
    }
  }
  void setCurrToOld2(int n) {
    if(n>=maximum_num_walkers) {return;}
    else {
      (walkers)[walker_size*n+data_displ[ELOC_OLD]+1] = (walkers)[walker_size*n+data_displ[ELOC]+1];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_A]+1] = (walkers)[walker_size*n+data_displ[OVLP_A]+1];
      (walkers)[walker_size*n+data_displ[OLD_OVLP_B]+1] = (walkers)[walker_size*n+data_displ[OVLP_B]+1];
    }
  }
  void setWalker(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]] = ob;
  }
  void setWalker(int n, ComplexType w0 ,ComplexType eloc) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[WEIGHT]] = w0; 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
  }
  void setWalker(int n, ComplexType w0 ,ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[WEIGHT]] = w0; 
    (walkers)[walker_size*n+data_displ[ELOC]] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]] = ob;
  }
  void setWalker2(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    if(n>=maximum_num_walkers) {return;} 
    (walkers)[walker_size*n+data_displ[ELOC]+1] = eloc; 
    (walkers)[walker_size*n+data_displ[OVLP_A]+1] = oa;
    (walkers)[walker_size*n+data_displ[OVLP_B]+1] = ob;
  }
  void getOldWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return ;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC_OLD]]; 
    oa = (walkers)[walker_size*n+data_displ[OLD_OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OLD_OVLP_B]];
  } 
  ComplexType* getWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC]]; 
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]];
    return &((walkers)[walker_size*n+data_displ[SM]]); 
  }
  ComplexType* getWalker2(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;} 
    eloc = (walkers)[walker_size*n+data_displ[ELOC]+1]; 
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]+1];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]+1];
    return &((walkers)[walker_size*n+data_displ[SM]+1]); 
  }
  ComplexType* getWalker(int n, ComplexType& w, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    if(n>=maximum_num_walkers) {return NULL;}
    w = (walkers)[walker_size*n+data_displ[WEIGHT]];
    eloc = (walkers)[walker_size*n+data_displ[ELOC]];
    oa = (walkers)[walker_size*n+data_displ[OVLP_A]];
    ob = (walkers)[walker_size*n+data_displ[OVLP_B]];
    return &((walkers)[walker_size*n+data_displ[SM]]);
  }

  void kill(int n) {
    walkers[walker_size*n+data_displ[INFO]] = ComplexType(-1.0);
    tot_num_walkers--;
  }

  bool isAlive(int n) {
    return ((n>=maximum_num_walkers)?false:(walkers[walker_size*n+data_displ[INFO]].real()>0)); // for now
  }
  
  void scaleWeight(RealType w0) {
    if(!head) return;
    for(int i=0; i<maximum_num_walkers; i++) 
      (walkers)[walker_size*i+data_displ[WEIGHT]] *= w0;
  } 

  void push_walkers_to_front();

  // useful when consistency breaks in algorithm
  void reset_walker_count() {
    maximum_num_walkers == walkers.size()/walker_size;
    tot_num_walkers=0;
    for(ComplexSMVector::iterator it=walkers.begin()+data_displ[INFO]; it<walkers.end(); it+=walker_size)
      if(it->real() > 0) tot_num_walkers++;
  }

  int single_walker_memory_usage() { return walker_memory_usage; } 
  int single_walker_size() { return walker_size; } 

  // copies a single walker at position "n" to/from buffer "data"
  void copy_to_buffer(int n, ComplexType* data) { 
    std::copy(walkers.values()+n*walker_size,walkers.values()+(n+1)*walker_size,data); 
  } 
  void copy_from_buffer(int n, ComplexType* data) { 
    std::copy(data,data+walker_size,walkers.values()+n*walker_size);
  } 

  // adds/removes n walkers from data. Puts/takes from the end of the list
  void pop_walkers(int n, ComplexType* data);
  void push_walkers(int n, ComplexType* data);

  void pop_walkers(int n, std::vector<ComplexType>& data) { pop_walkers(n,data.data()); }
  void push_walkers(int n, std::vector<ComplexType>& data) { push_walkers(n,data.data()); }

  int local_pair_branching(std::vector<int>& windx);
  // modifies the weight of walker at position branch_from and copies the resulting walker to branch_to
  void pair_branch(RealType w, int branch_from, int branch_to);
  // modifies the weight of walker at position "n" to "w" 
  // and adds "num" copies of the walker to the list 
  void branch(RealType w, int n, int num);


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
  int type, nrow, ncol; 
  int walker_size, data_displ[9], walker_memory_usage; 
  bool popcontrol_on_master;

  int targetN_per_TG;

  int tot_num_walkers;
  int maximum_num_walkers;

  ComplexType zero = ComplexType(0.0,0.0);
  RealType min_weight, max_weight, reset_weight;

  int extra_empty_spaces; 
  
  ComplexMatrix HFMat; 

  std::vector<int> empty_spots;

  bool head;

  // container with walker data 
  ComplexSMVector walkers; 

  MPI_Comm MPI_COMM_TG_LOCAL;
  MPI_Comm MPI_COMM_TG_LOCAL_HEADS; 
  int nproc_heads, rank_heads;
  
  std::vector<std::tuple<int,int>> outgoing, incoming; 
  std::vector<int> counts,displ;
  std::vector<char> bufferall;
  std::vector<char> bufferlocal;
  std::vector<ComplexType> commBuff;

  myTimer* LocalTimer;
  TimerList_t Timers;

};
}

#endif
