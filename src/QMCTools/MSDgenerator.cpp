//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include "OhmmsData/OhmmsElementBase.h"
#include "QMCTools/CasinoParser.h"
#include "QMCTools/GaussianFCHKParser.h"
#include "QMCTools/GamesXmlParser.h"
#include "QMCTools/GamesAsciiParser.h"
#include "QMCTools/BParser.h"
#include "Message/Communicate.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGenerator.h"


int main(int argc, char* argv[])
{
  int nea=16,neb=16,nca=0,ncb=0,nstates=32;
  double wgt=0.00,thr=0.0;
  bool useCSF=true;
  std::string name = "dummy.xml";
  OHMMS::Controller->initialize(argc,argv);
  if (OHMMS::Controller->rank() != 0) {
    outputManager.shutOff();
  }
  Random.init(0,1,-1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
  QMCGaussianParserBase::init();
  QMCGaussianParserBase *parser=0;
  parser = new GamesAsciiParser(argc,argv);
  parser->multideterminant = true;
  parser->SpinRestricted=true;
  parser->zeroCI=false;
  parser->orderByExcitation=false;
  parser->NumberOfEls=nca+ncb+nea+neb;
  parser->SpinMultiplicity=(parser->NumberOfEls-2*neb)+1;
  parser->NumberOfAlpha=nca+nea;
  parser->NumberOfBeta=ncb+neb;
  parser->ci_nca=nca;
  parser->ci_ncb=ncb;
  parser->ci_nea=nea;
  parser->ci_neb=neb;
  parser->ci_nstates=nstates;
  parser->ci_threshold=thr;
  parser->usingCSF=useCSF;
  std::vector<std::string> &CIalpha=parser->CIalpha;
  std::vector<std::string> &CIbeta=parser->CIbeta;
  std::vector<std::string> &CSFocc=parser->CSFocc;
  std::vector<std::vector<std::string> > &CSFalpha=parser->CSFalpha;
  std::vector<std::vector<std::string> > &CSFbeta=parser->CSFbeta;
  std::vector<std::vector<double> > &CSFexpansion=parser->CSFexpansion;
  std::vector<double> &CIcoeff=parser->CIcoeff;
  std::vector<std::pair<int,double> > &coeff2csf=parser->coeff2csf;
  if(useCSF)
  {
// dummy for now
    std::string alpha0(nstates,'0'),beta0(nstates,'0'),occ0(nstates,'0');
    std::vector<std::string> dummy;
    for(int i=0; i<nea; i++)
      alpha0[i]='1';
    for(int i=0; i<neb; i++)
      beta0[i]='1';
    for(int i=0; i<neb; i++)
      occ0[i]='2';
    for(int i=neb; i<nea; i++)
      occ0[i]='1';
    std::vector<double> vec;
    std::pair<int,double> dum(1,0.95);
    // HF
    coeff2csf.push_back(dum);
    CSFocc.push_back(occ0);
    CSFexpansion.push_back(vec);
    CSFexpansion.back().push_back(1.0);
    CSFalpha.push_back(dummy);
    CSFalpha.back().push_back(alpha0);
    CSFbeta.push_back(dummy);
    CSFbeta.back().push_back(beta0);
    // singles, assume multiplicity of 1
    for(int i=0; i<nea; i++)
    {
      for(int v=0; v<nstates-nea; v++)
      {
        coeff2csf.push_back(dum);
        coeff2csf.back().second = wgt;
        CSFocc.push_back(occ0);
        CSFocc.back().at(i) = '1';
        CSFocc.back().at(nea+v) = '1';
        CSFexpansion.push_back(vec);
        CSFexpansion.back().push_back(1.0/std::sqrt(2.0));
        CSFexpansion.back().push_back(1.0/std::sqrt(2.0));
        CSFalpha.push_back(dummy);
        CSFalpha.back().push_back(alpha0);
        CSFalpha.back().back().at(i)='0';
        CSFalpha.back().back().at(nea+v)='1';
        CSFalpha.back().push_back(alpha0);
        CSFalpha.back().back().at(i)='1';
        CSFalpha.back().back().at(nea+v)='0';
        CSFbeta.push_back(dummy);
        CSFbeta.back().push_back(beta0);
        CSFbeta.back().back().at(i)='1';
        CSFbeta.back().back().at(nea+v)='0';
        CSFbeta.back().push_back(beta0);
        CSFbeta.back().back().at(i)='0';
        CSFbeta.back().back().at(nea+v)='1';
      }
    }
  }
  else
  {
  }
  if(useCSF)
    parser->ci_size = CSFocc.size();
  else
    parser->ci_size = CIcoeff.size();
  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qm_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  xmlNodePtr multislaterdetPtr=NULL;
  multislaterdetPtr = parser->createMultiDeterminantSet();
  xmlAddChild(qm_root,multislaterdetPtr);
  xmlDocSetRootElement(doc, qm_root);
  xmlSaveFormatFile(name.c_str(),doc,1);
  xmlFreeDoc(doc);
  return 0;
}



