//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

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
  
  PostProcessor(const string& id,int ss,int se);

  ~PostProcessor() { }

  void set(const string& id,int ss,int se)
  {
    project_id   = id;
    series_start = ss;
    series_end   = se;
  }

  void put(xmlNodePtr cur,ParticleSetPool& ppool,
           WaveFunctionPool& wpool,HamiltonianPool& hpool);

  void postprocess();

 private:
  string project_id;
  int series_start;
  int series_end;

  vector<PostProcessorBase*> postprocessors;

  inline void add(PostProcessorBase* pp)
  {
    postprocessors.push_back(pp);
  }

};

}

#endif
