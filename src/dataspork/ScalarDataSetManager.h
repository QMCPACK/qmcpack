//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_SCALARDATAFILEMAGANER_H
#define QMCPLUSPLUS_SCALARDATAFILEMAGANER_H
#include <map>
///forward declation of ScalarObservable
class ScalarDataSet;

class ScalarDataSetManager
{

  std::vector<std::string> Observables;
  std::vector<std::string> Collectables;


  std::map<std::string,ScalarDataSet*> Files;

public:

  ScalarDataSetManager();
  ~ScalarDataSetManager();

  void addDataSet(const std::vector<std::string>& flist,
                  const std::string& gname, int first, int last);
  void addDataSet(const std::string& fname, int first, int last);

  void registerObservables(std::map<std::string, int>& olist);
  void registerCollectables(std::map<std::string, int>& olist);

  void write(std::ostream& os);

};
#endif
