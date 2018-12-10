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
    
    



/**@file BMakeFunc.h
 * @brief Declaration of BMakeFunc
 */
#ifndef QMCPLUSPLUS_TOOLS_BMAKEFUNC_H
#define QMCPLUSPLUS_TOOLS_BMAKEFUNC_H

/** base class to build a basisGroup
 */
struct BMakeFuncBase
{
  enum {GAUSSIANTYPE, SLATERTYPE, UNKNOWNTYPE};
  int N;
  int L;
  int M;
  int RadFuncType;
  std::vector<double> exponent;
  std::vector<double> contraction;
  std::vector<int> node;
  std::string BasisID;

  static std::vector<double> YlmNorm;
  static void init();

  BMakeFuncBase():
    N(1),L(0),M(0),RadFuncType(GAUSSIANTYPE),BasisID("none") {}

  xmlNodePtr createBasisGroup(bool shortform=false)
  {
    xmlNodePtr bptr = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
    std::ostringstream n,l,m;
    n<<N;
    l<<L;
    m<<M;
    xmlNewProp(bptr,(const xmlChar*)"rid",(const xmlChar*)BasisID.c_str());
    xmlNewProp(bptr,(const xmlChar*)"n",(const xmlChar*)n.str().c_str());
    xmlNewProp(bptr,(const xmlChar*)"l",(const xmlChar*)l.str().c_str());
    xmlNewProp(bptr,(const xmlChar*)"m",(const xmlChar*)m.str().c_str());
    switch(RadFuncType)
    {
    case(GAUSSIANTYPE):
      xmlNewProp(bptr,(const xmlChar*)"type",(const xmlChar*)"Gaussian");
      break;
    case(SLATERTYPE):
      xmlNewProp(bptr,(const xmlChar*)"type",(const xmlChar*)"Slater");
      break;
    default:
      xmlNewProp(bptr,(const xmlChar*)"type",(const xmlChar*)"Gaussian");
      break;
    }
    if(shortform)
      return bptr;
    for(int i=0; i<exponent.size(); i++)
    {
      std::ostringstream a,b,c;
      a.setf(std::ios::scientific, std::ios::floatfield);
      a.precision(8);
      b.setf(std::ios::scientific, std::ios::floatfield);
      b.precision(8);
      a << exponent[i];
      b << contraction[i];
      c << node[i];
      xmlNodePtr anode = xmlNewNode(NULL,(const xmlChar*)"radfunc");
      xmlNewProp(anode,(const xmlChar*)"exponent",(const xmlChar*)a.str().c_str());
      xmlNewProp(anode,(const xmlChar*)"contraction",(const xmlChar*)b.str().c_str());
      xmlNewProp(anode,(const xmlChar*)"node",(const xmlChar*)c.str().c_str());
      xmlAddChild(bptr,anode);
    }
    return bptr;
  }

  void addRadFunc(double e, double c, int n)
  {
    exponent.push_back(e);
    contraction.push_back(c);
    node.push_back(n);
  }

  void get(std::ostream& os)
  {
    os << BasisID << " n,l,m " << N << " " << L << " " << M << std::endl;
    for(int i=0; i<exponent.size(); i++)
    {
      std::cout << " " << exponent[i] << " " << contraction[i]
           << " " << node[i] << std::endl;
    }
  }
  virtual void put(std::vector<std::string>& words) = 0;
};


BMakeFuncBase* createBMakeFunc(int iflag);


#endif
