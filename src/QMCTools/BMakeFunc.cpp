//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file BMakeFunc.cpp
 * @brief specializations of BMakeFunc and builder function
 *
 * Consult makefun.f90 by M. Casula.
 */
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include "OhmmsData/libxmldefs.h"
#include "QMCTools/BMakeFunc.h"

std::vector<double> BMakeFuncBase::YlmNorm;
void BMakeFuncBase::init()
{
  YlmNorm.resize(10);
  for(int l=0; l<10; l++)
  {
    YlmNorm[l] = std::sqrt(4.0*M_PI/(2*static_cast<double>(l)+1));
  }
}

/** general basisGroup with one radial function */
template<int FLAG>
struct BMakeFunc: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    addRadFunc(atof(words[1].c_str()),1.0,1);
  }
};

/** case 34: c(exp(-dd1*x)+dd1*x*exp(-dd1*x))
 */
template<>
struct BMakeFunc<34>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=0;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),0.377964473009227,1);
    addRadFunc(atof(words[1].c_str()),0.654653670707977,1);
  }
};

/** case 10: r**2*exp(-z1*x)
 */
template<>
struct BMakeFunc<10>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=2;
    L=0;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,3);
  }
};

/** case 12: r**3*exp(-z1*x)
 */
template<>
struct BMakeFunc<12>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=3;
    L=0;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,4);
  }
};

/** case 20: exp(-z1*x)
 */
template<>
struct BMakeFunc<20>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=1;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,2);
  }
};

/** case 22: r*exp(-z1*x)
 */
template<>
struct BMakeFunc<22>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=2;
    L=1;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,3);
  }
};

/** case 30: r*exp(-z1*x)
 */
template<>
struct BMakeFunc<30>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=2;
    L=2;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),1.0,3);
  }
};

/** case 100: s Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<100>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=0;
    RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
};

/** case 103: p Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<103>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=1;
    RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
};

/** case 127: d Slater for J3 exp(-z*r)
 */
template<>
struct BMakeFunc<127>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=2;
    RadFuncType=SLATERTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
};

/** case 147: d Gaussian for J3 exp(-z*r*r)
 */
template<>
struct BMakeFunc<147>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=2;
    RadFuncType=GAUSSIANTYPE;
    addRadFunc(atof(words[1].c_str()),YlmNorm[L],1);
  }
};
/** case 200: constant
 */
template<>
struct BMakeFunc<200>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=0;
    L=0;
    RadFuncType=GAUSSIANTYPE;
    addRadFunc(0.0,YlmNorm[L],1);
  }
};

/** case 300: s with contracted GTO
 */
template<>
struct BMakeFunc<300>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=0;
    RadFuncType=GAUSSIANTYPE;
    int nc=(words.size()-1)/2;
    std::cout << " 300 = " << nc << std::endl;
    for(int ic=1; ic<=nc; ic++)
    {
      double e=atof(words[ic].c_str());
      double c=atof(words[ic+nc].c_str());
      std::cout << " 300:" << e << " " << c  << std::endl;
      addRadFunc(e,c,1);
    }
  }
};

/** case 400: p with contracted GTO
 */
template<>
struct BMakeFunc<400>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=1;
    RadFuncType=GAUSSIANTYPE;
    int nc=(words.size()-1)/2;
    std::cout << " 400 = " << nc << std::endl;
    for(int ic=1; ic<=nc; ic++)
    {
      double e=atof(words[ic].c_str());
      double c=atof(words[ic+nc].c_str());
      std::cout << " 400:" << e << " " << c  << std::endl;
      addRadFunc(e,c,1);
    }
  }
};

/** case 3000: s with contracted STO
 */
template<>
struct BMakeFunc<3000>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=1;
    L=0;
    RadFuncType=SLATERTYPE;
    int nc=(words.size()-1)/2;
    for(int ic=1; ic<=nc; ic++)
    {
      double e=atof(words[ic].c_str());
      double c=atof(words[ic+nc].c_str());
      int p=static_cast<int>(e/1000.);
      e -= static_cast<double>(p)*1000.;
      addRadFunc(e,c,p+1);
    }
  }
};

/** case 3100: p with  contracted STO
 */
template<>
struct BMakeFunc<3100>: public BMakeFuncBase
{
  void put(std::vector<std::string>& words)
  {
    N=2;
    L=1;
    RadFuncType=SLATERTYPE;
    int nc=(words.size()-1)/2;
    for(int ic=1; ic<=nc; ic++)
    {
      double e=atof(words[ic].c_str());
      double c=atof(words[ic+nc].c_str());
      int p=static_cast<int>(e/1000.);
      e -= static_cast<double>(p)*1000.;
      addRadFunc(e,c,p+2);
    }
  }
};

BMakeFuncBase* createBMakeFunc(int iflag)
{
  BMakeFuncBase* b;
  switch(iflag)
  {
  case(10):
    b=new BMakeFunc<10>;
    break;
  case(12):
    b=new BMakeFunc<12>;
    break;
  case(20):
    b=new BMakeFunc<20>;
    break;
  case(22):
    b=new BMakeFunc<22>;
    break;
  case(30):
    b=new BMakeFunc<30>;
    break;
  case(34):
    b=new BMakeFunc<34>;
    break;
  case(100):
    b=new BMakeFunc<100>;
    break;
  case(103):
    b=new BMakeFunc<103>;
    break;
  case(127):
    b=new BMakeFunc<127>;
    break;
  case(147):
    b=new BMakeFunc<147>;
    break;
  case(200):
    b=new BMakeFunc<200>;
    break;
  case(300):
    b=new BMakeFunc<300>;
    break;
  case(400):
    b=new BMakeFunc<400>;
    break;
  case(3000):
    b=new BMakeFunc<3000>;
    break;
  case(3100):
    b=new BMakeFunc<3100>;
    break;
  default:
    b=new BMakeFunc<0>;
  }
  return b;
}

