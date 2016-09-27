//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_POSTPROCESSOR_H
#define QMCPLUSPLUS_POSTPROCESSOR_H

#include<Estimators/PostProcessorBase.h>

namespace qmcplusplus
{

class ParticleSetPool;
class WaveFunctionPool;
class HamiltonianPool;

class PostProcessor 
{
 public:
  
  PostProcessor(const std::string& id,int ss,int se);

  ~PostProcessor() { }

  void set(const std::string& id,int ss,int se)
  {
    project_id   = id;
    series_start = ss;
    series_end   = se;
  }

  void put(xmlNodePtr cur,ParticleSetPool& ppool,
           WaveFunctionPool& wpool,HamiltonianPool& hpool);

  void postprocess();

 private:
  std::string project_id;
  int series_start;
  int series_end;

  std::vector<PostProcessorBase*> postprocessors;

  inline void add(PostProcessorBase* pp)
  {
    postprocessors.push_back(pp);
  }

};

}

#endif
