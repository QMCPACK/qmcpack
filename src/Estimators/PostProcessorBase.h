//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_POSTPROCESSOR_BASE_H
#define QMCPLUSPLUS_POSTPROCESSOR_BASE_H

#include<Message/Communicate.h>
#include <QMCHamiltonians/QMCHamiltonianBase.h>


namespace qmcplusplus
{

class PostProcessorBase
{
 public:
  typedef ParticleSet::ParticleLayout_t Lattice_t;
  typedef ParticleSet::PosType PosType;
  typedef ParticleSet::RealType RealType;

  enum{DIM=OHMMS_DIM};

  string type;
  string project_id;
  int series_start;
  int series_end;
  
  PostProcessorBase()
  {
    set("","",-1,-1);
  }

  ~PostProcessorBase() { }

  void set(const string& t,const string& id,int ss,int se)
  {
    type = t;
    project_id   = id;
    series_start = ss;
    series_end   = se;
  }

  virtual void put(xmlNodePtr cur)
  {
    APP_ABORT("PostProcessorBase::put is not implemented");
  }

  virtual void postprocess()
  {
    APP_ABORT("PostProcessorBase::postprocess is not implemented");
  }
};

}

#endif
