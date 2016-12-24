
// -*- C++ -*-
// /**@file WalkerHandlerBase 
//  * @brief Virtual Class for walker handlers. 
//   */

#ifndef QMCPLUSPLUS_AFQMC_SLATERDETWALKER_H
#define QMCPLUSPLUS_AFQMC_SLATERDETWALKER_H

#include<iostream>
#include <cstdlib>
#include <cstring>
#include"AFQMC/config.h"

namespace qmcplusplus
{

/*
 * Base (virtual) class for a walker.  
 */
class SlaterDetWalker: public AFQMCInfo
{

  typedef SlaterDetWalker WalkerPtr;
  typedef std::vector<WalkerPtr>::iterator WalkerIterator;

  public:

  /// constructor
  SlaterDetWalker():alive(false),init(false) {}

  SlaterDetWalker(const SlaterDetWalker& w)
  { 
    init=w.init;
    alive=w.alive;
    SlaterMat.resize(w.SlaterMat.rows(),w.SlaterMat.cols());
    SlaterMat=w.SlaterMat;
    weight=w.weight;
    energy_full=w.energy_full;
    eloc=w.eloc;
    eloc_old=w.eloc_old;
    overlap_alpha=w.overlap_alpha;
    overlap_beta=w.overlap_beta;     
    old_overlap_alpha=w.old_overlap_alpha;
    old_overlap_beta=w.old_overlap_beta;     
  }

  SlaterDetWalker& operator=(const SlaterDetWalker& w) 
  {
    if(this == &w) 
      return *this;

    init=w.init;
    alive=w.alive;
    SlaterMat.resize(w.SlaterMat.rows(),w.SlaterMat.cols());
    SlaterMat=w.SlaterMat;
    weight=w.weight;
    energy_full=w.energy_full;
    eloc=w.eloc;
    eloc_old=w.eloc_old;
    overlap_alpha=w.overlap_alpha;
    overlap_beta=w.overlap_beta;
    old_overlap_alpha=w.old_overlap_alpha;
    old_overlap_beta=w.old_overlap_beta; 
    return *this;
  }

  /// destructor
  ~SlaterDetWalker() {}

  int sizeForDump() { return sizeForComm();}
  // do I need to save *_old?
  int sizeForComm() {
    if(SlaterMat.size() == 0) return -1;
    return sizeof(ComplexType)*(14+SlaterMat.size());
  }

  //bool restartFromXML() {}
  //bool dumpToXML() {} 

  void clear() {}

  void unpackFromChar(char *arr, ComplexMatrix& A, ComplexType& w, ComplexType& el, ComplexType& ov )
  {
    int n=0;
    ComplexType o1;
    A.resize(SlaterMat.rows(),SlaterMat.cols());
    memcpy(&(w), arr ,sizeof(ComplexType)); n += sizeof(ComplexType)*2;
    memcpy(&(el), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType)*4;
    memcpy(&(o1), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType)*4;
    memcpy(&(ov), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType)*4;
    memcpy(A.data(), arr+n, A.size()*sizeof(ComplexType));    
    ov *= o1; 
  }  

  void restartFromChar(char *arr) {
    int n=0;
    memcpy(&(weight), arr ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(energy_full), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(eloc)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<1>(eloc)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<2>(eloc)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(eloc_old)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(overlap_alpha)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<1>(overlap_alpha)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<2>(overlap_alpha)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(old_overlap_alpha)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(overlap_beta)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<1>(overlap_beta)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<2>(overlap_beta)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(&(std::get<0>(old_overlap_beta)), arr+n ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(SlaterMat.data(), arr+n, SlaterMat.size()*sizeof(ComplexType)); 
  }

  void dumpToChar(char* arr ) {

    int n=0;
    memcpy(arr, &(weight) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(energy_full) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(eloc)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<1>(eloc)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<2>(eloc)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(eloc_old)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(overlap_alpha)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<1>(overlap_alpha)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<2>(overlap_alpha)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(old_overlap_alpha)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(overlap_beta)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<1>(overlap_beta)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<2>(overlap_beta)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, &(std::get<0>(old_overlap_beta)) ,sizeof(ComplexType)); n += sizeof(ComplexType);
    memcpy(arr+n, SlaterMat.data(), SlaterMat.size()*sizeof(ComplexType)); 

  }

  void initWalker(ComplexMatrix& S, bool a=true) 
  {
    init=true;
    alive=a;
    SlaterMat.resize(S.rows(),S.cols());
    SlaterMat=S;
    //weight=RealType(1.0);
    weight=ComplexType(1.0);
    eloc=std::forward_as_tuple(0.0,0.0,0.0);
    eloc_old=std::forward_as_tuple(0.0,0.0,0.0);
    overlap_alpha=std::forward_as_tuple(0.0,0.0,0.0);
    overlap_beta=std::forward_as_tuple(0.0,0.0,0.0);
    old_overlap_alpha=std::forward_as_tuple(0.0,0.0,0.0);
    old_overlap_beta=std::forward_as_tuple(0.0,0.0,0.0);
  }

  bool alive;
  bool init;

  // Storage for alpha and beta Sater matrix
  // Alpha from [0,NMO*NMO-1], beta from [NMO*NMO,2*NMO*NMO-1] 
  ComplexMatrix SlaterMat;

  // weight 
  //RealType weight; 
  ComplexType weight; 

  // <phi|H|psi> / <psi|psi> 
  ComplexType energy_full;

  // local energy: (importance sampling,  phaseless,  estimator) 
  std::tuple<ComplexType,ComplexType,ComplexType> eloc;
  std::tuple<ComplexType,ComplexType,ComplexType> eloc_old;

  // overlaps
  std::tuple<ComplexType,ComplexType,ComplexType> overlap_alpha;
  std::tuple<ComplexType,ComplexType,ComplexType> overlap_beta;
  std::tuple<ComplexType,ComplexType,ComplexType> old_overlap_alpha;
  std::tuple<ComplexType,ComplexType,ComplexType> old_overlap_beta;

};
}

#endif
