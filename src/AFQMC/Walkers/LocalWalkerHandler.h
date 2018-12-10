
// -*- C++ -*-
// /**@file WalkerHandler 
//  * @brief Virtual Class for walker handlers. 
//   */

#ifndef QMCPLUSPLUS_AFQMC_WALKERHANDLER_H
#define QMCPLUSPLUS_AFQMC_WALKERHANDLER_H

#include<random>

#include "OhmmsData/libxmldefs.h"
#include "io/hdf_archive.h"

#include"AFQMC/config.h"
#include"AFQMC/Walkers/WalkerHandlerBase.h"
#include<Message/MPIObjectBase.h>
#include "Message/CommOperators.h"

#include "AFQMC/Numerics/DenseMatrixOperations.h"


namespace qmcplusplus
{

/*
 * Class that contains and handles walkers.
 * Implements communication, load balancing, and I/O operations.   
 * Walkers are always accessed through the handler.
 */
class LocalWalkerHandler: public WalkerHandlerBase 
{

  typedef SlaterDetWalker Walker;
  typedef SlaterDetWalker* WalkerPtr;
  typedef std::vector<SlaterDetWalker>::iterator WalkerIterator;

  public:

  /// constructor
  LocalWalkerHandler(Communicate* c,  RandomGenerator_t* r): WalkerHandlerBase(c,r,std::string("LocalWalkerHandler")), 
      walkerSizeForCommunication(0),walkerSizeForDump(0),
      min_weight(0.05),max_weight(4.0)
      ,reset_weight(1.0),extra_empty_spaces(10)
  { }

  /// destructor
  ~LocalWalkerHandler() {}


  inline int size() { return walkers.size(); } 

  bool restartFromXML(); 
  bool dumpToXML();

  bool restartFromHDF5(int n, hdf_archive&, const std::string&, bool set_to_target); 
  bool dumpToHDF5(hdf_archive&, const std::string&); 

  bool dumpSamplesHDF5(hdf_archive& dump, int nW) {
    std::cerr<<" LocalWalkerHandler:dumpSamplesHDF5() not implemented. Not writing anything. \n";
    return true;
  }

  // reads xml and performs setup
  bool parse(xmlNodePtr cur); 

  // performs setup
  bool setup(int a,int b, int c, MPI_Comm, MPI_Comm , MPI_Comm,myTimer*);

  // cleans state of object. 
  //   -erases allocated memory 
  bool clean(); 

  // called at the beginning of each executable section
  void resetNumberWalkers(int n, bool a=true, ComplexMatrix* S=NULL); 

  inline bool initWalkers(int n) {
   resetNumberWalkers(n,true,&HFMat);
   // walkers.resize(n);
   // for(int i=0; i<n; i++) walkers[i].initWalker(HFMat);
   return true;   
  } 

  inline int numWalkers(bool includeDead=false) {
    int res=0;
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive || includeDead) res++; 
    return res;
  } 

  inline int numWalkers2() {
    int res=0;
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive && std::abs(it->weight) > 1e-6 ) res++; 
    return res;
  } 

  inline int GlobalPopulation() {
    std::vector<int> res(1);
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive) res[0]++; 
    myComm->gsum(res);
    return res[0];
  } 

  inline RealType GlobalWeight() {
    std::vector<RealType> res(1);
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive) res[0] += std::abs(it->weight);
    myComm->gsum(res);
    return res[0];
  }

  // load balancing algorithm
  void loadBalance(MPI_Comm comm); 

  // population control algorithm
  void popControl(MPI_Comm comm, std::vector<ComplexType>& curData); 

  void benchmark(std::string& blist,int maxnW,int delnW,int repeat) {}

  void setHF(const ComplexMatrix& HF);

  inline void Orthogonalize(int i) {
    DenseMatrixOperators::GeneralizedGramSchmidt(walkers[i].SlaterMat.data(),NAEA,NMO,NAEA);
    DenseMatrixOperators::GeneralizedGramSchmidt(walkers[i].SlaterMat.data()+NAEA*NMO,NAEA,NMO,NAEB);
  }

  inline void Orthogonalize() {
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive) { 
        DenseMatrixOperators::GeneralizedGramSchmidt(it->SlaterMat.data(),NAEA,NMO,NAEA);
        DenseMatrixOperators::GeneralizedGramSchmidt(it->SlaterMat.data()+NAEA*NMO,NAEA,NMO,NAEB);
      }
  }

  ComplexType* getSM(int n) { return (walkers[n].SlaterMat).data(); }
  ComplexType getWeight(int n) { return walkers[n].weight; }
  ComplexType getEloc(int n) { return std::get<0>(walkers[n].eloc); } 
  ComplexType getEloc2(int n) { return std::get<1>(walkers[n].eloc); } 
  ComplexType getOldEloc(int n) { return std::get<0>(walkers[n].eloc_old); }
  ComplexType getOldEloc2(int n) { return std::get<1>(walkers[n].eloc_old); }
  ComplexType getOvlpAlpha(int n) { return std::get<0>(walkers[n].overlap_alpha); }
  ComplexType getOvlpAlpha2(int n) { return std::get<1>(walkers[n].overlap_alpha); }
  ComplexType getOvlpBeta(int n) { return std::get<0>(walkers[n].overlap_beta); }
  ComplexType getOvlpBeta2(int n) { return std::get<1>(walkers[n].overlap_beta); }
  ComplexType getOldOvlpAlpha(int n) { return std::get<0>(walkers[n].old_overlap_alpha); }
  ComplexType getOldOvlpAlpha2(int n) { return std::get<1>(walkers[n].old_overlap_alpha); }
  ComplexType getOldOvlpBeta(int n) { return std::get<0>(walkers[n].old_overlap_beta); }
  ComplexType getOldOvlpBeta2(int n) { return std::get<1>(walkers[n].old_overlap_beta); }
  void setWeight(int n, ComplexType Q) { walkers[n].weight=Q; }
  void setEloc(int n, ComplexType Q) { std::get<0>(walkers[n].eloc) = Q; }
  void setEloc2(int n, ComplexType Q) { std::get<1>(walkers[n].eloc) = Q; }
  void setOldEloc(int n, ComplexType Q) { std::get<0>(walkers[n].eloc_old)=Q; }
  void setOldEloc2(int n, ComplexType Q) { std::get<1>(walkers[n].eloc_old)=Q; }
  void setOvlp(int n, ComplexType Q1, ComplexType Q2) {
    std::get<0>(walkers[n].overlap_alpha) = Q1;
    std::get<0>(walkers[n].overlap_beta) = Q2;
  }
  void setOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    std::get<1>(walkers[n].overlap_alpha) = Q1;
    std::get<1>(walkers[n].overlap_beta) = Q2;
  }
  void setOldOvlp(int n, ComplexType Q1, ComplexType Q2) {
    std::get<0>(walkers[n].old_overlap_alpha) = Q1;
    std::get<0>(walkers[n].old_overlap_beta) = Q2;
  }
  void setOldOvlp2(int n, ComplexType Q1, ComplexType Q2) {
    std::get<1>(walkers[n].old_overlap_alpha) = Q1;
    std::get<1>(walkers[n].old_overlap_beta) = Q2;
  }
  void setCurrToOld(int n) {
    std::get<0>(walkers[n].eloc_old) = std::get<0>(walkers[n].eloc);
    std::get<0>(walkers[n].old_overlap_alpha) = std::get<0>(walkers[n].overlap_alpha);
    std::get<0>(walkers[n].old_overlap_beta) = std::get<0>(walkers[n].overlap_beta);
  }
  void setCurrToOld2(int n) {
    std::get<1>(walkers[n].eloc_old) = std::get<1>(walkers[n].eloc);
    std::get<1>(walkers[n].old_overlap_alpha) = std::get<1>(walkers[n].overlap_alpha);
    std::get<1>(walkers[n].old_overlap_beta) = std::get<1>(walkers[n].overlap_beta);
  }
  void setWalker(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    std::get<0>(walkers[n].eloc) = eloc;
    std::get<0>(walkers[n].overlap_alpha) = oa;
    std::get<0>(walkers[n].overlap_beta) = ob;
  }
  void setWalker(int n, ComplexType w0, ComplexType eloc, ComplexType oa, ComplexType ob) {
    walkers[n].weight = w0;
    std::get<0>(walkers[n].eloc) = eloc;
    std::get<0>(walkers[n].overlap_alpha) = oa;
    std::get<0>(walkers[n].overlap_beta) = ob;
  }
  void setWalker(int n, ComplexType w0, ComplexType eloc) {
    walkers[n].weight = w0;
    std::get<0>(walkers[n].eloc) = eloc;
  }
  void setWalker2(int n, ComplexType eloc, ComplexType oa, ComplexType ob) {
    std::get<1>(walkers[n].eloc) = eloc;
    std::get<1>(walkers[n].overlap_alpha) = oa;
    std::get<1>(walkers[n].overlap_beta) = ob;
  }
  ComplexType* getWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    eloc = std::get<0>(walkers[n].eloc);
    oa = std::get<0>(walkers[n].overlap_alpha);
    ob = std::get<0>(walkers[n].overlap_beta);
    return (walkers[n].SlaterMat).data();
  }
  ComplexType* getWalker2(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    eloc = std::get<1>(walkers[n].eloc);
    oa = std::get<1>(walkers[n].overlap_alpha);
    ob = std::get<1>(walkers[n].overlap_beta);
    return (walkers[n].SlaterMat).data();
  }
  void getOldWalker(int n, ComplexType& eloc, ComplexType& oa, ComplexType& ob) {
    eloc = std::get<0>(walkers[n].eloc_old);
    oa = std::get<0>(walkers[n].old_overlap_alpha);
    ob = std::get<0>(walkers[n].old_overlap_beta);
  }
  ComplexType* getWalker(int n, ComplexType& w, ComplexType& eloc, ComplexType& oa, ComplexType& ob) { 
    w = walkers[n].weight;
    eloc = std::get<0>(walkers[n].eloc);
    oa = std::get<0>(walkers[n].overlap_alpha);
    ob = std::get<0>(walkers[n].overlap_beta);
    return (walkers[n].SlaterMat).data(); 
  }

  bool isAlive(int n) {
    return walkers[n].alive; 
  }

  void scaleWeight(RealType w0) {
    for(WalkerIterator it=walkers.begin(); it!=walkers.end(); ++it)
      if(it->alive) it->weight *= w0; 
  }

  //private:
 
  ComplexType min_weight, max_weight, reset_weight;

  int extra_empty_spaces; 
  
  ComplexMatrix HFMat; 

  // memory footprint (in bytes) of a walker for exchange between processors 
  int walkerSizeForCommunication;

  // memory footprint (in bytes) of a walker for dump 
  int walkerSizeForDump;

  // locations of empty spots in list 
  std::vector<int> emptyWalkers; 

  // container of walker pointers  
  std::vector<Walker> walkers; 

  private:

  inline WalkerIterator begin(int i) { return walkers.begin()+i; } 

  inline WalkerIterator begin() { return walkers.begin(); } 

  inline WalkerIterator end() { return walkers.end(); } 

  inline WalkerPtr getWalkerPtr(size_t i) { return &(walkers[i]); } 

  myTimer* LocalTimer;
     

};
}

#endif
