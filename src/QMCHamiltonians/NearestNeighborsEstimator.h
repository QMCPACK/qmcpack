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
    
    
#ifndef QMCPLUSPLUS_NEAREST_NEIGHBORS_ESTIMATOR_H
#define QMCPLUSPLUS_NEAREST_NEIGHBORS_ESTIMATOR_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>


namespace qmcplusplus
{



struct NeighborsTrace
{
  typedef DistanceTableData::ripair ripair;

  std::string ind_name;
  std::string dist_name;
  int count;
  ParticleSet* neighbors;
  ParticleSet* centers;
  ParticleSet* dtable_owner;
  int dtable_id;
  bool transposed;
  std::vector<ripair> ri;

  bool streaming_particles;
  TraceRequest& request;
  Array<TraceInt,2>*  index_sample;
  Array<TraceReal,2>* distance_sample;


  NeighborsTrace(int cnt, ParticleSet* p, ParticleSet* n,TraceRequest& req)
    : neighbors(p), centers(n), index_sample(0), distance_sample(0),request(req)
  {
    count=cnt;
    if(count>neighbors->getTotalNum())
      APP_ABORT("NeighborsTrace::NeighborsTrace  requested number of nearest neighbors is greater than the total number of neighbors");
    ind_name  = neighbors->parentName()+"_near_"+centers->parentName()+"_indices";
    dist_name = neighbors->parentName()+"_near_"+centers->parentName()+"_distances";
    dtable_id = neighbors->getTable(*centers);
    transposed = dtable_id<0;
    if(transposed)
    {
      dtable_id = centers->getTable(*neighbors);
      if(dtable_id<0)
        APP_ABORT("NeighborsTrace::NeighborsTrace  cannot find distance table for particlesets "+neighbors->getName()+" "+centers->getName());
      dtable_owner = centers;
    }
    else
    {
      dtable_owner = neighbors;
    }
    ri.resize(neighbors->getTotalNum());
  }


  inline void set_dynamic(ParticleSet& P)
  {
    bool match_centers   = centers->parentName()==P.parentName();
    bool match_neighbors = neighbors->parentName()==P.parentName();
    if(!match_centers && !match_neighbors)
    {
      APP_ABORT("NeighborsTrace::set_dynamic  particleset "+P.getName()+" is not the same as neighbors ("+neighbors->getName()+") or centers ("+centers->getName()+")");
    }
    if(match_centers)
    {
      if(dtable_owner==centers)
        dtable_owner = &P;
      centers = &P;
    }
    if(match_neighbors)
    {
      if(dtable_owner==neighbors)
        dtable_owner = &P;
      neighbors = &P;
    }
  }


  inline void contribute_arrays()
  {
    request.contribute_array(ind_name);
    request.contribute_array(dist_name);
  }

  inline void checkout_arrays(TraceManager& tm)
  {
    streaming_particles = request.streaming_array(ind_name) || request.streaming_array(dist_name);
    if( streaming_particles)
    {
      index_sample    = tm.checkout_int<2>( ind_name, *centers,count);
      distance_sample = tm.checkout_real<2>(dist_name,*centers,count);
    }
  }


  inline void delete_arrays()
  {
    if( streaming_particles)
    {
      delete index_sample;
      delete distance_sample;
      index_sample = 0;
      distance_sample = 0;
    }
  }


  inline void sample()
  {
    if( streaming_particles)
    {
      DistanceTableData& dtable = *dtable_owner->DistTables[dtable_id];
      dtable.check_neighbor_size(ri,transposed);
      for(int i=0; i<centers->getTotalNum(); ++i)
      {
        dtable.nearest_neighbors(i,count,ri,transposed);
        for(int n=0; n<count; ++n)
          (*distance_sample)(i,n) = ri[n].first;
        for(int n=0; n<count; ++n)
          (*index_sample)(i,n)    = ri[n].second;
      }
    }
  }
};




class NearestNeighborsEstimator : public QMCHamiltonianBase
{
private:
  typedef std::map<std::string,ParticleSet*> PSPool;

  std::vector<NeighborsTrace*> ntraces;
  PSPool& psetpool;

public:

  NearestNeighborsEstimator(PSPool& psp)
    : psetpool(psp)
  {
    myName = "NearestNeighbors";
  }


  ~NearestNeighborsEstimator()
  {
    delete_iter(ntraces.begin(),ntraces.end());
  }


  bool put(xmlNodePtr cur)
  {
    xmlNodePtr element = cur->children;
    while(element!=NULL)
    {
      std::string name((const char*)element->name);
      if(name=="neighbor_trace")
      {
        OhmmsAttributeSet eattrib;
        int count = -1;
        std::string neighbors_name = "none";
        std::string centers_name   = "none";
        eattrib.add(count,"count");
        eattrib.add(neighbors_name,"neighbors");
        eattrib.add(centers_name,"centers");
        eattrib.put(element);
        check_attribute_presence(count,neighbors_name,centers_name);
        ParticleSet* neighbors = get_particleset(neighbors_name);
        ParticleSet* centers   = get_particleset(centers_name);
        check_attribute_values(count,neighbors,centers,neighbors_name,centers_name);
        ntraces.push_back(new NeighborsTrace(count,neighbors,centers,request));
      }
      else if(name=="text")
      {
      }
      else
      {
        APP_ABORT("NearestNeighborsEstimator::put  "+name+" is not a valid sub-element of NearestNeighbors\n  valid options are: nearest");
      }
      element=element->next;
    }
    return true;
  }


  inline void check_attribute_presence(int count,const std::string& neighbors_name,const std::string& centers_name)
  {
    if(count==-1)
      APP_ABORT("NearestNeighborsEstimator::put  element <neighbor_trace/> must have attribute 'count' (int)");
    if(neighbors_name=="none")
      APP_ABORT("NearestNeighborsEstimator::put  element <neighbor_trace/> must have attribute 'neighbors' ( std::string)");
    if(centers_name=="none")
      APP_ABORT("NearestNeighborsEstimator::put  element <neighbor_trace/> must have attribute 'centers' ( std::string)");
  }


  inline void check_attribute_values(int count,const ParticleSet* neighbors,const ParticleSet* centers,const std::string& neighbors_name,const std::string& centers_name)
  {
    if(!neighbors)
      APP_ABORT("NearestNeighborsEstimator::put  <neighbor_trace/> attribute neighbors="+neighbors_name+" does not correspond to a known particleset");
    if(!centers)
      APP_ABORT("NearestNeighborsEstimator::put  <neighbor_trace/> attribute centers="+centers_name+" does not correspond to a known particleset");
    if(count<0)
      APP_ABORT("NearestNeighborsEstimator::put  element <neighbor_trace/> attribute 'count' (int) must be greater than zero");
    if(count>neighbors->getTotalNum()){
      app_log()<<"centers   = "<<centers->parentName()<< std::endl;
      app_log()<<"neighbors = "<<neighbors->parentName()<< std::endl;
      app_log()<<"        max # of neighbors = "<<neighbors->getTotalNum()<< std::endl;
      app_log()<<"  requested # of neighbors = "<<count<< std::endl;
      APP_ABORT("NearestNeighborsEstimator::put  element <neighbor_trace/> attribute 'count' (int) must not be greater than total number of particles in neighbor particleset");
    }
  }


  ParticleSet* get_particleset( std::string& psname)
  {
    if(psetpool.find(psname)==psetpool.end())
    {
      app_log()<<"  ParticleSet "<<psname<<" does not exist"<< std::endl;
      APP_ABORT("NearestNeighborsEstimator::put");
    }
    return psetpool[psname];
  }


  virtual QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    NearestNeighborsEstimator* clone = new NearestNeighborsEstimator(psetpool);
    for(int n=0; n<ntraces.size(); ++n)
    {
      const NeighborsTrace& nt = *ntraces[n];
      NeighborsTrace* ntrace = new NeighborsTrace(nt.count,nt.neighbors,nt.centers,clone->request);
      ntrace->set_dynamic(qp);
      clone->ntraces.push_back(ntrace);
    }
    return clone;
  }


  virtual void contribute_particle_quantities()
  {
    for(int n=0; n<ntraces.size(); ++n)
      ntraces[n]->contribute_arrays();
  }


  virtual void checkout_particle_quantities(TraceManager& tm)
  {
    for(int n=0; n<ntraces.size(); ++n)
      ntraces[n]->checkout_arrays(tm);
    streaming_particles = false;
    for(int n=0; n<ntraces.size(); ++n)
      streaming_particles |= ntraces[n]->streaming_particles;
  }


  virtual void delete_particle_quantities()
  {
    for(int n=0; n<ntraces.size(); ++n)
      ntraces[n]->delete_arrays();
  }


  virtual Return_t evaluate(ParticleSet& P)
  {
    if( streaming_particles){
      for(int n=0; n<ntraces.size(); ++n)
        ntraces[n]->sample();
    }
    return 0.0;
  }

  virtual Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  // does not produce a scalar value
  virtual void contribute_scalar_quantities()               { }
  virtual void checkout_scalar_quantities(TraceManager& tm) { }
  virtual void collect_scalar_quantities()                  { }
  virtual void delete_scalar_quantities()                   { }

  // is not a scalar observable or collectable (produces no value, only a trace stream)
  virtual bool get(std::ostream& os) const
  {
    return true;
  }
  virtual void resetTargetParticleSet(ParticleSet& P) {}
  virtual void addObservables(PropertySetType& plist, BufferType& collectables) {}
  virtual void registerObservables(std::vector<observable_helper*>& h5desc,hid_t gid) const {}
  virtual void setObservables(PropertySetType& plist) {}
  virtual void setParticlePropertyList(PropertySetType& plist, int offset) {}
  virtual void setHistories(Walker_t& ThisWalker) {}
};

}



#endif
