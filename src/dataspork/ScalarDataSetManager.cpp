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
    
    


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include "ScalarDataSet.h"
#include "ScalarDataSetManager.h"

ScalarDataSetManager::ScalarDataSetManager() { }
ScalarDataSetManager::~ScalarDataSetManager()
{
  std::map<std::string,ScalarDataSet*>::iterator
  it(Files.begin()),it_end(Files.end());
  while(it != it_end)
  {
    delete (*it).second;
    ++it;
  }
}

void
ScalarDataSetManager::registerObservables(std::map<std::string, int>& olist)
{
  Observables.resize(olist.size());
  std::map<std::string,int>::iterator it(olist.begin()), it_end(olist.end());
  while(it != it_end)
  {
    Observables[(*it).second]=(*it).first;
    ++it;
  }
}

void
ScalarDataSetManager::registerCollectables(std::map<std::string, int>& olist)
{
  Collectables.resize(olist.size());
  std::map<std::string,int>::iterator it(olist.begin()), it_end(olist.end());
  while(it != it_end)
  {
    Collectables[(*it).second]=(*it).first;
    ++it;
  }
}

void ScalarDataSetManager::addDataSet(const std::string& fname, int first, int last)
{
  Files[fname] = new ScalarDataSet(fname,first,last);
}

/** Add multiple data to the list and create a merged data
 */
void ScalarDataSetManager::addDataSet(const std::vector<std::string>& flist,
                                      const std::string& gname, int first, int last)
{
  int nrow=100000000;
  for(int i=0; i<flist.size(); i++)
  {
    ScalarDataSet* a=new ScalarDataSet(flist[i],first,last);
    nrow = std::min(nrow,a->rows());
    Files[flist[i]]=a;
  }
  std::map<std::string,ScalarDataSet*>::iterator dit(Files.begin());
  ScalarDataSet* g = (*dit).second->clone(nrow);
  g->group(Files);
  Files[gname]=g;
}

void ScalarDataSetManager::write(std::ostream& os)
{
  std::map<std::string,ScalarDataSet*>::iterator it(Files.begin()),it_end(Files.end());
  while(it != it_end)
  {
    std::cout << "============================================\n"
         << "\t source = " << (*it).first
         << "\n============================================\n";
    os << "  <dataset src=\""<<(*it).first<<"\" grouped=\"";
    if((*it).second->Grouped)
      os<<"yes\">\n";
    else
      os << "no\">\n";
    ScalarDataSet& sdf(*((*it).second));
    os << "    <coldesc "
       <<"first=\""<<sdf.FirstRow << "\" "
       <<"last=\"" <<sdf.LastRow << "\" "
       <<"size=\"" <<sdf.NumRows << "\"/>\n";
    std::vector<std::string>::iterator oit(Collectables.begin()), oit_end(Collectables.end());
    while(oit != oit_end)
    {
      sdf.writeCollect(os,*oit++);
    }
    std::cout << std::endl;
    oit=Observables.begin();
    oit_end=Observables.end();
    while(oit != oit_end)
    {
      sdf.writeObservable(os,*oit++);
    }
    os << "  </dataset>" << std::endl;
    //write(os, (*it).first, *((*it).second));
    ++it;
  }
}
