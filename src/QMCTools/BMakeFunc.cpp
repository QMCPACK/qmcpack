//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file BMakeFunc.cpp
 * @brief specializations of BMakeFunc and builder function
 *
 * Consult makefun.f90 by M. Casula. 
 */
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;
#include "OhmmsData/libxmldefs.h"
#include "QMCTools/BMakeFunc.h"

vector<double> BMakeFuncBase::YlmNorm;
void BMakeFuncBase::init() {
  YlmNorm.resize(10);
  for(int l=0; l<10; l++) {
    YlmNorm[l] = std::sqrt(4.0*M_PI/(2*static_cast<double>(l)+1));
  }
}

/** general basisGroup with one radial function */
template<int FLAG>
struct BMakeFunc: public BMakeFuncBase { 
  void put(vector<string>& words) {
    addRadFunc(atof(words[1].c_str()),1.0,1);
  }
}; 

/** case 34: c(exp(-dd1*x)+dd1*x*exp(-dd1*x))
 */
template<>
struct BMakeFunc<34>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=1; L=0; RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),0.377964473009227,1);
    addRadFunc(atof(words[1].c_str()),0.654653670707977,1);
  }
}; 

/** case 10: r**2*exp(-z1*x)
 */
template<>
struct BMakeFunc<10>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=2; L=0; RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,3);
  }
}; 

/** case 12: r**3*exp(-z1*x)
 */
template<>
struct BMakeFunc<12>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=3; L=0; RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,4);
  }
}; 

/** case 20: exp(-z1*x)
 */
template<>
struct BMakeFunc<20>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=1; L=1; RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,2);
  }
}; 

/** case 22: r*exp(-z1*x)
 */
template<>
struct BMakeFunc<22>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=2; L=1; RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,3);
  }
}; 

/** case 100: s Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<100>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=1; L=0; RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
}; 

/** case 103: p Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<103>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=1; L=1; RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
}; 

/** case 147: d Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<147>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=1; L=2; RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
}; 
/** case 200: constant
 */
template<>
struct BMakeFunc<200>: public BMakeFuncBase { 
  void put(vector<string>& words) {
    N=0; L=0; RadFuncType=GAUSSIANTYPE;
    addRadFunc(0.0,YlmNorm[L],1);
  }
}; 

BMakeFuncBase* createBMakeFunc(int iflag) {
  BMakeFuncBase* b;
  switch(iflag) {
    case(10): b=new BMakeFunc<10>; break;
    case(12): b=new BMakeFunc<12>; break;
    case(20): b=new BMakeFunc<20>; break;
    case(22): b=new BMakeFunc<22>; break;
    case(34): b=new BMakeFunc<34>; break;
    case(100): b=new BMakeFunc<100>; break;
    case(103): b=new BMakeFunc<103>; break;
    case(147): b=new BMakeFunc<147>; break;
    case(200): b=new BMakeFunc<200>; break;
    default: b=new BMakeFunc<0>;
  }
  return b;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
