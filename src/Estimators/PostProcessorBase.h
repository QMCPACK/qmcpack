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

  std::string type;
  std::string project_id;
  int series_start;
  int series_end;
  
  PostProcessorBase()
  {
    set("","",-1,-1);
  }

  ~PostProcessorBase() { }

  void set(const std::string& t,const std::string& id,int ss,int se)
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
