//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Norbert Podhorszki, pnorbert@ornl.gov, Oak Ridge National Laboratory
//                    Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_TRACEMANAGER_H
#define QMCPLUSPLUS_TRACEMANAGER_H


#if !defined(DISABLE_TRACEMANAGER)


#include <Configuration.h>
#include <OhmmsData/OhmmsElementBase.h>
#include <OhmmsData/AttributeSet.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Particle/ParticleSet.h>
#include <Utilities/IteratorUtility.h>
#include <Message/Communicate.h>
#include <io/hdf_archive.h>
#include <map>
#include <set>
#include <algorithm>

#ifdef HAVE_ADIOS
#include "ADIOS/ADIOSTrace.h"
#include <std::string.h>
//#include "adios.h"
//#include "adios_read.h"
//#include "ADIOS/ADIOS_config.h"
#include <boost/lexical_cast.hpp>
//#ifdef IO_PROFILE
//#include "ADIOS/ADIOS_profile.h"
//#endif
//#ifdef ADIOS_VERIFY
//#include "ADIOS/ADIOS_verify.h"
//#endif
#endif


namespace qmcplusplus
{

//#define TRACE_CHECK

const unsigned int DMAX=4;
typedef long   TraceInt;
typedef OHMMS_PRECISION TraceReal;
typedef std::complex<TraceReal> TraceComp;


struct TraceQuantity
{
  std::string name;
  bool default_quantity       ;
  bool combined_quantity      ;
  bool scalar_available       ;
  bool array_available        ;
  bool scalar_stream_requested;
  bool array_stream_requested ;
  bool scalar_write_requested ;
  bool array_write_requested  ;
  bool stream_scalar          ;
  bool stream_array           ;
  bool write_scalar           ;
  bool write_array            ;

  inline TraceQuantity()
  {
    default_quantity        = false;
    combined_quantity       = false;
    scalar_available        = false;
    array_available         = false;
    scalar_stream_requested = false;
    array_stream_requested  = false;
    scalar_write_requested  = false;
    array_write_requested   = false;
    stream_scalar           = false;
    stream_array            = false;
    write_scalar            = false;
    write_array             = false;
  }

  inline void incorporate(const TraceQuantity& other)
  {
    if(name!=other.name)
    {
      APP_ABORT("TraceQuantity::incorporate\n  cannot merge quantities with differing names\n  names: "+name+" "+other.name);
    }
    default_quantity        |= other.default_quantity       ;
    combined_quantity       |= other.combined_quantity      ;
    scalar_available        |= other.scalar_available       ;
    array_available         |= other.array_available        ;
    scalar_stream_requested |= other.scalar_stream_requested;
    array_stream_requested  |= other.array_stream_requested ;
  }
};


//means of control of traces
//  which quantities are available, streaming, writing, etc
//  handles global (user + all observable) trace request 
//medium of exchange of trace information between trace manager and observables
//  handles individual trace requests from observables
struct TraceRequest
{
  //switch to allow or disallow streams/writes of requested quantities
  bool allow_streams;
  bool allow_writes;

  //switch to allow or disallow scalar/array streaming/writing
  bool scalars_on;
  bool arrays_on;

  //all scalar/array quantities should be provided for streaming
  bool stream_all_scalars;
  bool stream_all_arrays;

  //all scalar/array quantities should be provided for writing
  bool write_all_scalars;
  bool write_all_arrays;

  //records whether default scalars/arrays should stream/write
  bool streaming_default_scalars;
  bool streaming_default_arrays ;
  bool writing_default_scalars  ;
  bool writing_default_arrays   ;


  //quantities made available or requested for streaming or writing
  std::map<std::string,TraceQuantity> quantities;

  //dependency lists for combined quantities
  std::map<std::string,std::set<std::string> > combined_dependencies;

  //used to screen checked out quantities for writing
  std::string scalar_domain;

  inline TraceRequest()
  {
    reset();
  }

  inline void reset()
  {
    allow_streams             = false;
    allow_writes              = false;
    scalars_on                = false;
    arrays_on                 = false;
    stream_all_scalars        = false;
    stream_all_arrays         = false;
    write_all_scalars         = false;
    write_all_arrays          = false;
    streaming_default_scalars = false;
    streaming_default_arrays  = false;
    writing_default_scalars   = false;
    writing_default_arrays    = false;

    quantities.clear();
    combined_dependencies.clear();
  }

  inline void set_scalar_domain(const std::string& domain)
  {
    scalar_domain = domain;
  }

  inline bool screen_sample(const std::string& domain,const std::string& name,bool& write)
  {
    bool scalar  = domain == scalar_domain;
    bool present = quantities.find(name)!=quantities.end();
    bool stream  = false;
    write   = false;
    if(present)
    {
      TraceQuantity& q = quantities[name];
      if(scalar)
      {
        stream = q.stream_scalar;
        write  = q.write_scalar;
      }
      else
      {
        stream = q.stream_array;
        write  = q.write_array;
      }
    }
    return stream;
  }

//Contributor API (QMCHamiltonianBase and others)
  //declare that scalars are available for a quantity
  inline void contribute_scalar(const std::string& name,bool default_quantity=false)
  {
    guarantee_presence(name);
    quantities[name].scalar_available = true;
    if(default_quantity)
      quantities[name].default_quantity = true;
  }

  //declare that arrays are available for a quantity
  inline void contribute_array(const std::string& name,bool default_quantity=false)
  {
    guarantee_presence(name);
    quantities[name].array_available = true;
    if(default_quantity)
      quantities[name].default_quantity = true;
  }

  //declare that a combined quantity could be made available
  inline void contribute_combined(const std::string& name,std::vector<std::string>&deps,bool scalar=false,bool array=false,bool default_quantity=false)
  {
    guarantee_presence(name,true);
    TraceQuantity& q = quantities[name];
    q.combined_quantity = true;
    q.scalar_available  = scalar;
    q.array_available   = array;
    if(default_quantity)
      q.default_quantity = true;
    //combined_dependencies[name] = deps;
  }

  //declare that a scalar quantity is desired for streaming
  inline void request_scalar(const std::string& name,bool write=false)
  {
    guarantee_presence(name);
    quantities[name].scalar_stream_requested = true;
    if(write)
      quantities[name].scalar_write_requested = true;
  }

  //declare that a array quantity is desired for streaming
  inline void request_array(const std::string& name,bool write=false)
  {
    guarantee_presence(name);
    quantities[name].array_stream_requested = true;
    if(write)
      quantities[name].array_write_requested = true;
  }

  //query whether a scalar quantity should/will be streaming
  inline bool streaming_scalar(const std::string& name)
  {
    check_presence(name);
    return quantities[name].stream_scalar;
  }

  //query whether a array quantity should/will be streaming
  inline bool streaming_array(const std::string& name)
  {
    check_presence(name);
    return quantities[name].stream_array;
  }

  //query whether a quantity will be streaming in any fashion
  inline bool streaming(const std::string& name)
  {
    check_presence(name);
    TraceQuantity& q = quantities[name];
    return q.stream_scalar || q.stream_array;
  }


//TraceManager API
  //declare that scalar quantities are desired for streaming
  inline void request_scalar(const std::set<std::string>& names,bool write=false)
  {
    std::set<std::string>::iterator name;
    for(name=names.begin();name!=names.end();++name)
      request_scalar(*name,write);
  }

  //declare that array quantities are desired for streaming
  inline void request_array(const std::set<std::string>& names,bool write=false)
  {
    std::set<std::string>::iterator name;
    for(name=names.begin();name!=names.end();++name)
      request_array(*name,write);
  }

  //merge in all quantities from a contributor request
  inline void incorporate(TraceRequest& other)
  {
    std::map<std::string,TraceQuantity>::iterator it;
    for(it=other.quantities.begin();it!=other.quantities.end();++it)
    {
      const TraceQuantity& q = it->second;
      if(quantity_present(q.name))
        quantities[q.name].incorporate(q);
      else
        quantities[q.name] = q;
    }
    std::map<std::string,std::set<std::string> >::iterator d;
    for(d=other.combined_dependencies.begin();d!=other.combined_dependencies.end();++d)
    {
      const std::string& name = d->first;
      std::set<std::string>& deps  = d->second;
      if(combined_dependencies.find(name)!=combined_dependencies.end())
        combined_dependencies[name].insert(deps.begin(),deps.end());
      else
        combined_dependencies[name] = deps;
    }
  }

  //balance requests with availability and turn streaming/writing on/off
  inline void determine_stream_write()
  {
    std::map<std::string,TraceQuantity>::iterator it;
    for(it=quantities.begin();it!=quantities.end();++it)
    {
      TraceQuantity& q = it->second;

      q.stream_scalar = allow_streams && scalars_on && q.scalar_available 
                        && (q.scalar_stream_requested || stream_all_scalars);

      q.write_scalar  = q.stream_scalar
                        && (q.scalar_write_requested || write_all_scalars);

      q.stream_array  = allow_streams && arrays_on  && q.array_available 
                        && (q.array_stream_requested || stream_all_arrays);

      q.write_array   = q.stream_array
                        && (q.array_write_requested  || write_all_arrays);
    }
    // default quantities stream and write if any others do
    streaming_default_scalars = false;
    writing_default_scalars   = false;
    streaming_default_arrays  = false;
    writing_default_arrays    = false;
    for(it=quantities.begin();it!=quantities.end();++it)
    {
      TraceQuantity& q = it->second;
      streaming_default_scalars |= q.stream_scalar;
      writing_default_scalars   |= q.write_scalar;
      streaming_default_arrays  |= q.stream_array;
      writing_default_arrays    |= q.write_array;
    }
    for(it=quantities.begin();it!=quantities.end();++it)
    {
      TraceQuantity& q = it->second;
      if(q.default_quantity)
      {
        q.stream_scalar = streaming_default_scalars;
        q.write_scalar  = writing_default_scalars  ;
        q.stream_array  = streaming_default_arrays ;
        q.write_array   = writing_default_arrays   ;
      }
    }
    // if any combined quantities are streaming, their dependencies must also
    for(it=quantities.begin();it!=quantities.end();++it)
    {
      TraceQuantity& q = it->second;
      if(q.combined_quantity)
      {
        std::set<std::string>& deps = combined_dependencies[q.name];
        std::set<std::string>::iterator name;
        if(q.stream_scalar || q.stream_array)
        {
          for(name=deps.begin();name!=deps.end();++name)
          {
            check_presence(*name);
            TraceQuantity& qd = quantities[*name];
            qd.stream_scalar |= (qd.scalar_available && q.stream_scalar);
            qd.stream_array  |= (qd.array_available  && q.stream_array );
          }
        }
      }
    }
    combined_dependencies.clear();
  }

  //relay updated streaming information to contributor
  inline void relay_stream_info(TraceRequest& other)
  {
    other.allow_streams             = allow_streams     ;
    other.allow_writes              = allow_writes      ;
    other.scalars_on                = scalars_on        ;
    other.arrays_on                 = arrays_on         ;
    other.stream_all_scalars        = stream_all_scalars;
    other.stream_all_arrays         = stream_all_arrays ;
    other.write_all_scalars         = write_all_scalars ;
    other.write_all_arrays          = write_all_arrays  ;
    other.streaming_default_scalars = streaming_default_scalars;
    other.streaming_default_arrays  = streaming_default_arrays ;
    other.writing_default_scalars   = writing_default_scalars  ;
    other.writing_default_arrays    = writing_default_arrays   ;
    std::map<std::string,TraceQuantity>::iterator it;
    for(it=other.quantities.begin();it!=other.quantities.end();++it)
    {
      TraceQuantity& q = it->second;
      check_presence(q.name);
      q = quantities[q.name];
    }
    other.combined_dependencies.clear();
  }

  inline void report()
  {
    //app_log()<<"\n  TraceRequest"<< std::endl;
    app_log()<<"    allow_streams             = "<< allow_streams             << std::endl;
    app_log()<<"    allow_writes              = "<< allow_writes              << std::endl;
    app_log()<<"    scalars_on                = "<< scalars_on                << std::endl;
    app_log()<<"    arrays_on                 = "<< arrays_on                 << std::endl;
    app_log()<<"    stream_all_scalars        = "<< stream_all_scalars        << std::endl;
    app_log()<<"    stream_all_arrays         = "<< stream_all_arrays         << std::endl;
    app_log()<<"    write_all_scalars         = "<< write_all_scalars         << std::endl;
    app_log()<<"    write_all_arrays          = "<< write_all_arrays          << std::endl;
    app_log()<<"    streaming_default_scalars = "<< streaming_default_scalars << std::endl;
    app_log()<<"    streaming_default_arrays  = "<< streaming_default_arrays  << std::endl;
    app_log()<<"    writing_default_scalars   = "<< writing_default_scalars   << std::endl;
    app_log()<<"    writing_default_arrays    = "<< writing_default_arrays    << std::endl;

    write_selected("scalars available","scalar_available");
    write_selected("arrays available","array_available");

    write_selected("scalar streams requested","scalar_stream_requested");
    write_selected("array streams requested","array_stream_requested");

    write_selected("scalar writes requested","scalar_write_requested");
    write_selected("array writes requested","array_write_requested");

    write_selected("scalar streams occurring","stream_scalar");
    write_selected("array streams occurring","stream_array");

    write_selected("scalar writes occurring","write_scalar");
    write_selected("array writes occurring","write_array");
  }

  inline void write_selected(const std::string& header,const std::string& selector)
  {
    app_log()<<"    "<<header<<":";
    int n=0;
    std::map<std::string,TraceQuantity>::iterator it;
    for(it=quantities.begin();it!=quantities.end();++it)
    {
      TraceQuantity& q = it->second;
      bool selected;
      if(selector=="scalar_available")
        selected = q.scalar_available;
      else if(selector=="array_available")
        selected = q.array_available;
      else if(selector=="scalar_stream_requested")
        selected = q.scalar_stream_requested;
      else if(selector=="array_stream_requested")
        selected = q.array_stream_requested;
      else if(selector=="scalar_write_requested")
        selected = q.scalar_write_requested;
      else if(selector=="array_write_requested")
        selected = q.array_write_requested;
      else if(selector=="stream_scalar")
        selected = q.stream_scalar;
      else if(selector=="stream_array")
        selected = q.stream_array;
      else if(selector=="write_scalar")
        selected = q.write_scalar;
      else if(selector=="write_array")
        selected = q.write_array;
      else
        APP_ABORT("TraceRequest::write_selected  unrecognized selector: "+selector);
      if(selected)
      {
        if(n%5==0)
          app_log()<< std::endl<<"      ";
        n++;
        app_log()<<" "<<q.name;
      }
    }
    app_log()<< std::endl;
  }


//private (internal) API
  //query whether a quantity is present
  inline bool quantity_present(const std::string& name)
  {
    return quantities.find(name)!=quantities.end();
  }

  //create a quantity if it is not present
  inline void guarantee_presence(const std::string& name,bool combined=false)
  {
    if(!quantity_present(name))
    {
      TraceQuantity q;
      q.name = name;
      quantities[name] = q;
    }
    if(combined)
      if(combined_dependencies.find(name)==combined_dependencies.end())
      {
        std::set<std::string> stmp;
        combined_dependencies[name] = stmp;
      }
  }

  //abort if a quantity is not present
  inline void check_presence(const std::string& name)
  {
    if(!quantity_present(name))
    {
      APP_ABORT("TraceRequest::check_presence  quantity "+name+" is not present");
    }
  }

  //query whether any quantities are streaming
  inline bool streaming()
  {
    return streaming_default_scalars || streaming_default_arrays;
  }

  //query whether any quantities are writing
  inline bool writing()
  {
    return writing_default_scalars || writing_default_arrays;
  }

  //query whether any scalar quantities are streaming
  inline bool streaming_scalars()
  {
    return streaming_default_scalars;
  }

  //query whether any array quantities are streaming
  inline bool streaming_arrays()
  {
    return streaming_default_arrays;
  }
};



template<typename T>
struct TraceSample
{
  std::string domain;
  std::string name;
  int index;
  bool array_trace;
  int dimension;
  int size;
  int unit_size;
  int data_size;
  TinyVector<int,DMAX> shape;
  std::vector<T>&  sample;
  bool write;
  int buffer_start,buffer_end;
  std::map<std::string,TraceInt>  meta_int;
  std::map<std::string,TraceReal> meta_real;
  std::map<std::string,std::string> meta_string;
  bool verbose;

  inline TraceSample(const std::string& sdomain,const std::string& sname,int sindex,int sdim,std::vector<T>& ssample)
    : sample(ssample), verbose(false)
  {
    initialize(sdomain,sname,sindex,sdim);
  }


  inline TraceSample(const std::string& sdomain,const std::string& sname,int sindex,int sdim,TinyVector<int,DMAX> sshape,std::vector<T>& ssample)
    : sample(ssample), verbose(false)
  {
    initialize(sdomain,sname,sindex,sdim);
    shape = sshape;
    size  = sample.size();
    check_shape();
  }


  inline void initialize(const std::string& sdomain,const std::string& sname,int sindex,int sdim)
  {
    domain       = sdomain,
    name         = sname;
    dimension    = sdim;
    index        = sindex;
    array_trace  = false;
    write        = false;
    buffer_start = -1;
    buffer_end   = -1;
  }


  inline void set_unit_size(int usize)
  {
    unit_size = usize;
  }


  inline void set_data_size()
  {
    data_size = size*unit_size;
  }


  inline void check_shape()
  {
    bool correct_shape     = true;
    bool correct_dimension = dimension<=DMAX;
    if(correct_dimension)
    {
      int tsize = 1;
      for(int d=0; d<dimension; ++d)
      {
        tsize*=shape[d];
        correct_dimension=correct_dimension && shape[d]>0;
      }
      correct_shape = tsize==size;
    }
    if(!correct_dimension)
      APP_ABORT("TraceSample::check_shape dimension of sample array is incorrect");
    if(!correct_shape)
      APP_ABORT("TraceSample::check_shape shape and size of sample array do not match");
  }


  inline bool same_shape(TraceSample<T>* other)
  {
    bool same = dimension==other->dimension && size==other->size;
    if(same)
      for(int d=0; d<dimension; ++d)
        same = same && shape[d]==other->shape[d];
    return same;
  }

  virtual bool is_combined()
  {
    return false;
  }

  inline void set_buffer_range(int& bstart)
  {
    set_data_size();
    if(write)
    {
      buffer_start = bstart;
      buffer_end   = bstart + data_size;
      bstart = buffer_end;
    }
  }


  inline T sum()
  {
    T s(0);
    for(int i=0; i<sample.size(); ++i)
      s+=sample[i];
    return s;
  }

  inline void write_summary(int ind=-1,std::string pad="  ")
  {
    std::string pad2 = pad+"  ";
    if(ind==-1)
      app_log()<<pad<<" TraceSample "<<name<< std::endl;
    else
      app_log()<<pad<<ind<<" TraceSample "<<name<< std::endl;
    app_log()<<pad2<<"domain       = "<< domain        << std::endl;
    app_log()<<pad2<<"name         = "<< name          << std::endl;
    app_log()<<pad2<<"index        = "<< index         << std::endl;
    app_log()<<pad2<<"array_trace  = "<< array_trace<< std::endl;
    app_log()<<pad2<<"dimension    = "<< dimension     << std::endl;
    app_log()<<pad2<<"size         = "<< size          << std::endl;
    app_log()<<pad2<<"unit_size    = "<< unit_size     << std::endl;
    app_log()<<pad2<<"data_size    = "<< data_size     << std::endl;
    app_log()<<pad2<<"shape        = "<< shape         << std::endl;
    app_log()<<pad2<<"write        = "<< write         << std::endl;
    app_log()<<pad2<<"buffer range = ["<<buffer_start<<","<<buffer_end<<")"<< std::endl;
  }

};




template<typename T>
struct CombinedTraceSample : public TraceSample<T>
{
  bool combined;
  std::vector<TraceReal>          weights;
  std::vector<TraceSample<T>*> components;


  inline CombinedTraceSample(const std::string& sdomain,const std::string& sname,int sindex,int sdim,std::vector<T>& ssample)
    : TraceSample<T>(sdomain,sname,sindex,sdim,ssample)
  {
    reset();
  }


  inline CombinedTraceSample(const std::string& sdomain,const std::string& sname,int sindex,int sdim,TinyVector<int,DMAX> sshape,std::vector<T>& ssample)
    : TraceSample<T>(sdomain,sname,sindex,sdim,sshape,ssample)
  {
    reset();
  }

  virtual bool is_combined()
  {
    return true;
  }

  inline void reset()
  {
    combined = false;
  }


  inline void add_component(TraceSample<T>* component,TraceReal weight)
  {
    if(components.size()==0)
    {
      this->dimension   = component->dimension;
      this->size        = component->size;
      this->shape       = component->shape;
      this->data_size   = component->data_size;
      this->array_trace = component->array_trace;
      this->sample.resize(component->size);
    }
    else if(!this->same_shape(component))
    {
      APP_ABORT("CombinedTraceSample::add_component  attempted to add a different shape component\n  my domain: "+this->domain+"\n  my name: "+this->name+"\n  component domain: "+component->domain+"\n  component name: "+component->name);
    }
    else if(this->domain!=component->domain)
      APP_ABORT("CombinedTraceSample::add_component  attempted to add a different domain component\n  my domain: "+this->domain+"\n  my name: "+this->name+"\n  component domain: "+component->domain+"\n  component name: "+component->name);
    weights.push_back(weight);
    components.push_back(component);
  }


  inline void combine()
  {
    fill(this->sample.begin(),this->sample.end(),T(0));
    for(int i=0; i<components.size(); ++i)
    {
      T weight = weights[i];
      std::vector<T>& component = components[i]->sample;
      for(int j=0; j<this->sample.size(); ++j)
        this->sample[j] += weight*component[j];
    }
    combined = true;
  }


  inline void write_summary_combined(int ind,std::string pad="  ")
  {
    std::string pad2 = pad+"  ";
    std::string pad3 = pad2+"  ";
    app_log()<<pad<<ind<<" CombinedTraceSample "<<this->name<< std::endl;
    app_log()<<pad2<<"domain      = "<<this->domain<< std::endl;
    app_log()<<pad2<<"ncomponents = "<<components.size()<< std::endl;
    app_log()<<pad2<<"components"<< std::endl;
    for(int i=0; i<components.size(); ++i)
    {
      TraceSample<T>& c = *components[i];
      app_log()<<pad3<<c.name<<" "<<c.index<<" "<<weights[i] << std::endl;
    }
    app_log()<<pad2<<"end components"<< std::endl;
    app_log()<<pad2<<"vector address = "<<(size_t)&this->sample<< std::endl;
  }
};




template<typename T>
bool TraceSample_comp(TraceSample<T>* left,TraceSample<T>* right)
{
  return left->data_size < right->data_size;
}





template<typename T>
struct TraceSamples
{
  std::vector<TraceSample<T>*>         samples;
  std::map<std::string,std::map<std::string,int> >    sample_indices;
  std::vector<TraceSample<T>*>         ordered_samples;
  std::vector<CombinedTraceSample<T>*> combined_samples;
  std::vector<std::vector<T>*>              combined_sample_vectors;
  bool verbose;

  inline TraceSamples():verbose(false) {}

  inline ~TraceSamples()
  {
    finalize();
  }


  inline void set_verbose(bool v)
  {
    verbose = v;
  }


  inline int size()
  {
    return samples.size();
  }


  inline void assign_sample_index(const std::string& domain,const std::string& name,int index,std::string label="")
  {
    if(sample_indices.count(domain)>0 && sample_indices[domain].count(name)>0)
    {
      APP_ABORT("TraceSamples::checkout "+label+" variable "+name+" already exists in domain "+domain);
    }
    else
    {
      sample_indices[domain][name]=index;
    }
  }

  template<int D>
  inline Array<T,D>* checkout_array(const std::string& domain,const std::string& name,TinyVector<int,DMAX> shape)
  {
    int index = samples.size();
    assign_sample_index(domain,name,index,"array");
    Array<T,D>* a = new Array<T,D>(shape.data());
    TraceSample<T>* s = new TraceSample<T>(domain,name,index,D,shape,a->storage());
    samples.push_back(s);
    if(verbose)
      app_log()<<"TraceSamples::checkout_array  "<<domain<<" "<<name<<" "<<index<< std::endl;
    return a;
  }


  template<int D>
  inline Array<T,D>* checkout_array(const ParticleSet& P,const std::string& name,TinyVector<int,DMAX> shape)
  {
    const std::string& domain = P.parentName();
    int index = samples.size();
    assign_sample_index(domain,name,index,"array");
    Array<T,D>* a = new Array<T,D>(shape.data());
    TraceSample<T>* s = new TraceSample<T>(domain,name,index,D,shape,a->storage());
    samples.push_back(s);
    s->array_trace = true;
    if(verbose)
      app_log()<<"TraceSamples::checkout_array  "<<domain<<" "<<name<<" "<<index<< std::endl;
    return a;
  }


  inline TraceSample<T>* get_trace(const std::string& domain, const std::string& name)
  {
    TraceSample<T>* ts=NULL;
    for(int i=0; i<samples.size(); ++i)
    {
      TraceSample<T>& tsc = *samples[i];
      if(tsc.domain==domain && tsc.name==name)
      {
        ts = &tsc;
        break;
      }
    }
    if(ts==NULL)
      APP_ABORT("TraceSamples::get_trace  failed to get trace for quantity "+name+" in domain "+domain);
    return ts;
  }


  inline CombinedTraceSample<T>* get_combined_trace(const std::string& domain, const std::string& name)
  {
    CombinedTraceSample<T>* ts=NULL;
    for(int i=0; i<combined_samples.size(); ++i)
    {
      CombinedTraceSample<T>& tsc = *combined_samples[i];
      if(tsc.domain==domain && tsc.name==name)
      {
        ts = &tsc;
        break;
      }
    }
    if(ts==NULL)
      APP_ABORT("TraceSamples::get_combined_trace  failed to get trace for quantity "+name+" in domain "+domain);
    return ts;
  }


  inline bool make_combined_trace(const std::string& name,std::vector<std::string>& names,std::vector<TraceReal>& weights)
  {
    bool created = false;
    if(names.size()!=weights.size())
      APP_ABORT("TraceSamples::make_combined_trace  names and weights must be the same size");
    std::map<std::string,std::map<std::string,int> >::iterator it;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      std::string domain = it->first;
      std::map<std::string,int>& indices = it->second;
      bool any_present = false;
      for(int i=0; i<names.size(); ++i)
        any_present = any_present || indices.count(names[i])>0;
      if(any_present)
      {
        int index = samples.size();
        std::vector<T>* sample = new std::vector<T>;
        CombinedTraceSample<T>* combined = new CombinedTraceSample<T>(domain,name,index,0,*sample);
        for(int i=0; i<names.size(); ++i)
        {
          if(indices.count(names[i])>0)
          {
            TraceSample<T>* component = samples[indices[names[i]]];
            combined->add_component(component,weights[i]);
          }
        }
        assign_sample_index(domain,name,index);
        samples.push_back(combined);
        combined_samples.push_back(combined);
        combined_sample_vectors.push_back(sample);
        created = true;
      }
    }
    return created;
  }


  inline void set_unit_size(int usize)
  {
    for(int i=0; i<samples.size(); i++)
      samples[i]->set_unit_size(usize);
  }


  inline void screen_writes(TraceRequest& request)
  {
    for(int i=0; i<samples.size(); i++)
    {
      TraceSample<T>& s = *samples[i];
      bool stream = request.screen_sample(s.domain,s.name,s.write);
      if(verbose)
        app_log()<<"TraceRequest screening "<<s.name<<" in domain "<<s.domain<<". stream: "<<stream<<" write: "<<s.write<< std::endl;
      if(!stream && !s.is_combined())
        app_log()<<"warning: quantity "+s.name+" in domain "+s.domain+" was not requested but is streaming anyway"<< std::endl;
    }
  }


  inline void order_by_size()
  {
    for(int i=0; i<samples.size(); i++)
      samples[i]->set_data_size();
    ordered_samples.resize(samples.size());
    copy(samples.begin(),samples.end(),ordered_samples.begin());
    sort(ordered_samples.begin(),ordered_samples.end(),TraceSample_comp<T>);
  }


  inline void set_buffer_ranges(int& starting_index)
  {
    for(int i=0; i<ordered_samples.size(); i++)
    {
      TraceSample<T>& sample = *ordered_samples[i];
      sample.set_buffer_range(starting_index);
    }
  }


  inline int total_size()
  {
    int s=0;
    for(int i=0; i<samples.size(); i++)
      s+=samples[i]->sample.size()*samples[i]->unit_size;
    return s;
  }


  inline int min_buffer_index()
  {
    int min_index = 2000000000;
    for(int i=0; i<samples.size(); i++)
      min_index = std::min(min_index,samples[i]->buffer_start);
    return min_index;
  }


  inline int max_buffer_index()
  {
    int max_index = -1;
    for(int i=0; i<samples.size(); i++)
      max_index = std::max(max_index,samples[i]->buffer_end);
    return max_index;
  }


  inline void combine_samples()
  {
    for(int i=0; i<combined_samples.size(); ++i)
      combined_samples[i]->combine();
  }


  inline void reset_combined_samples()
  {
    for(int i=0; i<combined_samples.size(); ++i)
      combined_samples[i]->reset();
  }


  inline void finalize()
  {
    //note combined_samples pointers are a subset of those in samples
    //  and so only samples needs to be deleted
    delete_iter(samples.begin(),samples.end());
    delete_iter(combined_sample_vectors.begin(),combined_sample_vectors.end());
    samples.resize(0);
    ordered_samples.resize(0);
    combined_samples.resize(0);
    combined_sample_vectors.resize(0);
    sample_indices.clear();
  }


  inline void register_hdf_data(hdf_archive& f)
  {
    std::map<std::string,std::map<std::string,int> >::iterator it;
    std::map<std::string,int>::iterator it2;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      const std::string& domain = it->first;
      std::map<std::string,int>& indices = it->second;
      f.push(domain);
      for(it2=indices.begin(); it2!=indices.end(); ++it2)
      {
        const std::string& quantity = it2->first;
        TraceSample<T>& sample = *samples[it2->second];
        if(sample.write)
        {
          f.push(quantity);
          f.write(sample.dimension,   "dimension");
          f.write(sample.shape,       "shape"    );
          f.write(sample.size,        "size"     );
          f.write(sample.unit_size,   "unit_size");
          f.write(sample.buffer_start,"row_start");
          f.write(sample.buffer_end,  "row_end"  );
          f.pop();
        }
      }
      f.pop();
    }
  }


#ifdef HAVE_ADIOS
  // determine what adios type to use for output (of one variable)
  inline std::string get_adios_type ( std::string type, std::string domain, std::string name, bool complex=0)
  {
        /* If we keep double/float, do this */
        /*
        if (!type.compare("real") && sizeof(T) == sizeof(double))
            vartype = "double";
        else
            vartype = type;

        if (complex) { 
            if (!type.compare("real") && sizeof(T) == sizeof(std::complex<double>))
                vartype = "complex double";
            else
                vartype = "complex";
        }
        */ 
        std::string vartype;
        if (domain == "scalars") {
            /* /scalars/... variables keep their types */
            if (type == "real" && sizeof(T) == sizeof(double))
                vartype = "double";
            else
                vartype = type;
        } else if (complex) { 
            vartype = "complex"; /* complex float even if TraceCompl is std::complex<double> */
        } else {
            vartype = type; /* int or real (i.e. 4 byte ADIOS types, even if TraceReal is double */
        }
        return vartype;
  }

  inline void register_adios_data(ADIOS::Trace& at, std::string type, bool complex=0)
  {
    std::string varpath;
    std::string vartype;
    int dims[4] = {0,0,0,0};

    for(int i=0; i<ordered_samples.size(); ++i) 
    {
        const struct TraceSample<T>& s = *ordered_samples[i];
        /*
        app_log()<<"domain         = "<< s.domain        << std::endl;
        app_log()<<"name           = "<< s.name          << std::endl;
        app_log()<<"index          = "<< s.index         << std::endl;
        app_log()<<"array_trace    = "<< s.array_trace<< std::endl;
        app_log()<<"dimension      = "<< s.dimension     << std::endl;
        app_log()<<"size           = "<< s.size          << std::endl;
        app_log()<<"unit_size      = "<< s.unit_size     << std::endl;
        app_log()<<"data_size      = "<< s.data_size     << std::endl;
        app_log()<<"shape          = "<< s.shape         << std::endl;
        app_log()<<"buffer range   = ["<<s.buffer_start<<","<<s.buffer_end<<")"<< std::endl;
        */

        //jtk: attempt to restrict adios writes in the correct way, please check!
        if(s.write)
        {
          varpath = "/"+s.domain+"/"+s.name;
          for (int j=0; j < s.dimension; j++)
          {
            dims[j] = s.shape[j]; // * s.unit_size;
          }

          vartype = get_adios_type (type, s.domain, s.name, complex);
          at.define_var (varpath, s.dimension, dims, vartype);
        }
    }
  }
#endif


  inline void write_summary( std::string type,std::string pad="  ")
  {
    std::string pad2 = pad +"  ";
    std::string pad3 = pad2+"  ";
    std::string pad4 = pad3+"  ";
    app_log()<<pad<<"TraceSamples<"<<type<<">"<< std::endl;
    app_log()<<pad2<<"nsamples          = "<<samples.size()<< std::endl;
    app_log()<<pad2<<"ncombined_samples = "<<combined_samples.size()<< std::endl;
    app_log()<<pad2<<"sample_indices"<< std::endl;
    std::map<std::string,std::map<std::string,int> >::iterator it;
    std::map<std::string,int>::iterator it2;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      const std::string& domain = it->first;
      std::map<std::string,int>& indices = it->second;
      app_log()<<pad3<<domain<< std::endl;
      for(it2=indices.begin(); it2!=indices.end(); ++it2)
        app_log()<<pad4<<it2->first<<" = "<<it2->second<< std::endl;
    }
    app_log()<<pad2<<"end sample_indices"<< std::endl;
    app_log()<<pad2<<"combined_sample_vectors = ";
    for(int i=0; i<combined_sample_vectors.size(); ++i)
      app_log()<<(size_t)combined_sample_vectors[i]<<" ";
    app_log()<< std::endl;
    app_log()<<pad2<<"combined_samples"<< std::endl;
    for(int i=0; i<combined_samples.size(); ++i)
      combined_samples[i]->write_summary_combined(i,pad3);
    app_log()<<pad2<<"end combined_samples"<< std::endl;
    app_log()<<pad2<<"samples"<< std::endl;
    for(int i=0; i<ordered_samples.size(); ++i)
      ordered_samples[i]->write_summary(i,pad3);
    //for(int i=0; i<samples.size(); ++i)
    //  samples[i]->write_summary(i,pad3);
    app_log()<<pad2<<"end samples"<< std::endl;
    app_log()<<pad<<"end TraceSamples<"<<type<<">"<< std::endl;
  }


  inline void user_report(const std::string& type,const std::string& pad="  ")
  {
    std::string pad2 = pad+"  ";
    app_log()<<pad<<type<<" traces provided by estimators"<< std::endl;
    std::map<std::string,std::map<std::string,int> >::iterator it;
    std::map<std::string,int>::iterator it2;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      const std::string& domain = it->first;
      std::map<std::string,int>& indices = it->second;
      app_log()<<pad2<<"domain "<<domain<<":  ";
      int n=0;
      for(it2=indices.begin(); it2!=indices.end(); ++it2)
      {
        if(n%5==0)
          app_log()<< std::endl<<pad2<<"  ";
        n++;
        const std::string& quantity = it2->first;
        app_log()<<quantity<<" ";
      }
      app_log()<< std::endl;
    }
  }
};




template<typename T>
struct TraceBuffer
{
  bool has_complex;
  TraceSamples<T>* samples;
  TraceSamples<std::complex<T> >* complex_samples;
  std::string type;
  Array<T,2> buffer;
  bool verbose;

  //hdf variables
  std::string top;
  hsize_t dims[2];
  hsize_t hdf_file_pointer;


  TraceBuffer()
    : samples(0),complex_samples(0),verbose(false)
  {
    type = "?";
    has_complex = false;
    reset();
  }


  inline void set_verbose(bool v)
  {
    verbose = v;
  }


  inline void set_type( std::string stype)
  {
    type = stype;
    top = type + "_data";
  }


  inline void reset()
  {
    buffer.resize(0,buffer.size(1));
  }


  inline void set_samples(TraceSamples<T>& s)
  {
    samples = &s;
  }


  inline void set_samples(TraceSamples<std::complex<T> >& s)
  {
    complex_samples = &s;
    has_complex = true;
  }


  inline void make_combined_trace(const std::string& name,std::vector<std::string>& names,std::vector<TraceReal>& weights)
  {
    bool created_real = samples->make_combined_trace(name,names,weights);
    if(has_complex)
    {
      bool created_complex = complex_samples->make_combined_trace(name,names,weights);
      if(created_real && created_complex)
        APP_ABORT("TraceBuffer<"+type+">::make_combined_trace\n  cannot create real and complex combined traces for the same quantity\n  attempted for quantity "+name);
    }
  }


  inline void order_and_resize()
  {
    //put the sample data in size order
    samples->set_unit_size(1);
    samples->order_by_size();
    if(has_complex)
    {
      complex_samples->set_unit_size(2);
      complex_samples->order_by_size();
    }
    //assign buffer ranges to each sample
    int sample_size = 0;
    samples->set_buffer_ranges(sample_size);
    if(has_complex)
      complex_samples->set_buffer_ranges(sample_size);
#if defined(TRACE_CHECK)
    test_buffer_write(sample_size);
#endif
    //resize the buffer
    int nsamples_init = 1;
    buffer.resize(nsamples_init,sample_size);
  }


  inline bool same_as(TraceBuffer<T>& ref)
  {
    return buffer.size(1)==ref.buffer.size(1);
  }


  inline void collect_sample()
  {
    if(verbose)
      app_log()<<" TraceBuffer<"<<type<<">::collect_sample()"<< std::endl;
    //make more room, if necessary
    int nrows    = buffer.size(0);
    int row_size = buffer.size(1);
    if(row_size>0)
    {
      //make space for the row, if necessary
      int current_row = nrows;
      nrows++;
      buffer.resize(nrows,row_size);
      if(verbose)
        app_log()<<"  increasing # of rows to "<<nrows<< std::endl;
      //combine samples
      samples->combine_samples();
      if(has_complex)
        complex_samples->combine_samples();
      //collect data from all samples into the buffer row
      int offset = current_row*row_size;
      {
        int boffset;
        std::vector<TraceSample<T>*>& ordered_samples = samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<T>& tsample = *ordered_samples[s];
          if(tsample.write)
          {
            std::vector<T>& sample = tsample.sample;
            boffset = offset + tsample.buffer_start;
            for(int i=0; i<sample.size(); ++i)
            {
              buffer(boffset+i) = sample[i];
            }
          }
        }
      }
      if(has_complex)
      {
        int boffset;
        std::vector<TraceSample<std::complex<T> >*>& ordered_samples = complex_samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<std::complex<T> >& tsample = *ordered_samples[s];
          if(tsample.write)
          {
            std::vector<std::complex<T> >& sample = tsample.sample;
            boffset = offset + tsample.buffer_start;
            for(int i=0,ib=0; i<sample.size(); ++i,ib+=2)
            {
              buffer(boffset+ib)   = sample[i].real();
              buffer(boffset+ib+1) = sample[i].imag();
            }
          }
        }
      }
      //reset combined samples so they can be recombined on the next step
      samples->reset_combined_samples();
      if(has_complex)
        complex_samples->reset_combined_samples();
#if defined(TRACE_CHECK)
      test_buffer_collect(current_row);
#endif
    }
  }


  inline void write()
  {
    APP_ABORT("TraceBuffer::write has not yet been implemented");
  }


  inline void write_summary( std::string pad="  ")
  {
    std::string pad2=pad+"  ";
    app_log()<<pad<<"TraceBuffer<"<<type<<">"<< std::endl;
    app_log()<<pad2<<"nrows       = "<< buffer.size(0)<< std::endl;
    app_log()<<pad2<<"row_size    = "<< buffer.size(1)<< std::endl;
    app_log()<<pad2<<"has_complex = "<< has_complex   << std::endl;
    samples->write_summary(type,pad2);
    if(has_complex)
      complex_samples->write_summary("complex "+type,pad2);
    app_log()<<pad<<"end TraceBuffer<"<<type<<">"<< std::endl;
  }


  inline void user_report(const std::string& pad="  ")
  {
    samples->user_report(type,pad);
    if(has_complex)
      complex_samples->user_report("complex "+type,pad);
  }

  inline void register_hdf_data(hdf_archive& f)
  {
    f.push(top);
    f.push("layout");
    samples->register_hdf_data(f);
    if(has_complex)
      complex_samples->register_hdf_data(f);
    f.pop();
    f.pop();
    if(!f.group_id.empty())
      APP_ABORT("TraceBuffer<"+type+">::register_hdf_data() some hdf groups are still open at the end of registration");
    hdf_file_pointer = 0;
  }


#ifdef HAVE_ADIOS
  inline void register_adios_data(ADIOS::Trace& at)
  {
    samples->register_adios_data(at, type, false);
    if(has_complex) 
    {
      complex_samples->register_adios_data(at,type,true);
    }
  }
#endif


  inline void write_hdf(hdf_archive& f)
  {
    write_hdf(f,hdf_file_pointer);
  }


  inline void write_hdf(hdf_archive& f,hsize_t& file_pointer)
  {
    if(verbose)
      app_log()<<"TraceBuffer<"<<type<<">::write_hdf() "<<file_pointer<<" "<<buffer.size(0)<<" "<<buffer.size(1)<< std::endl;
    dims[0] = buffer.size(0);
    dims[1] = buffer.size(1);
    if(dims[0]>0)
    {
      f.push(top);
      h5d_append(f.top(), "traces", file_pointer,
                 buffer.dim(), dims, buffer.data());
      f.pop();
    }
    f.flush();
  }




  inline void test_buffer_write(int sample_size)
  {
    //check that the size is correct
    int ssize = samples->total_size();
    if(has_complex)
      ssize+=complex_samples->total_size();
    if(sample_size!=ssize)
    {
      app_log()<<"sample_size = "<<sample_size<<"\ntotal_size = "<<ssize<< std::endl;
      APP_ABORT("TraceBuffer::test_buffer_write sample_size and total_size do not match");
    }
    //check that buffer indices fall in the expected range
    int nsamples=samples->size();
    int min_index = samples->min_buffer_index();
    int max_index = samples->max_buffer_index();
    if(has_complex)
    {
      nsamples+=complex_samples->size();
      min_index = std::min(min_index,complex_samples->min_buffer_index());
      max_index = std::max(max_index,complex_samples->max_buffer_index());
    }
    if(nsamples>0)
    {
      if(min_index!=0)
        APP_ABORT("TraceBuffer::test_buffer_write min_index!=0\n  min_index="<<min_index);
      if(max_index!=sample_size)
        APP_ABORT("TraceBuffer::test_buffer_write max_index!=sample_size");
      //check that no overlap exists in writes to buffer
      Array<int,2> test_buffer;
      test_buffer.resize(1,sample_size);
      fill(test_buffer.begin(),test_buffer.end(),0);
      int row = 0;
      int row_size = test_buffer.size(1);
      int offset = row*row_size;
      int* loc1 = &test_buffer(offset);
      int* loc2 = &test_buffer(row,0);
      if(loc1!=loc2)
        APP_ABORT("TraceBuffer::test_buffer_write serialized buffer offset improperly computed");
      {
        int boffset;
        std::vector<TraceSample<T>*>& ordered_samples = samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<T>& tsample = *ordered_samples[s];
          std::vector<T>& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0; i<sample.size(); ++i)
            test_buffer(boffset+i) = 1;
        }
      }
      if(has_complex)
      {
        int boffset;
        std::vector<TraceSample<std::complex<T> >*>& ordered_samples = complex_samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<std::complex<T> >& tsample = *ordered_samples[s];
          std::vector<std::complex<T> >& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0,ib=0; i<sample.size(); ++i,ib+=tsample.unit_size)
          {
            test_buffer(boffset+ib)   = 1;
            test_buffer(boffset+ib+1) = 1;
          }
        }
      }
      //app_log()<<"test_buffer:"<< std::endl;
      //for(int i=0;i<row_size;++i)
      //  app_log()<<"  "<<i<<"  "<<test_buffer(row,i)<< std::endl;
      for(int i=0; i<row_size; ++i)
        if(!test_buffer(row,i))
          APP_ABORT("TraceBuffer::test_buffer_write write to row is not contiguous");
    }
  }


  inline void test_buffer_collect(int current_row)
  {
    if(verbose)
      app_log()<<"TraceBuffer::test_buffer_collect"<< std::endl;
    std::string scalars = "scalars";
    std::map<std::string,std::map<std::string,int> >::iterator dom;
    std::map<std::string,int>::iterator var;
    std::map<std::string,std::map<std::string,int> >& domains = samples->sample_indices;
    std::vector<TraceSample<T>*>& tsamples = samples->samples;
    std::map<std::string,int>& scalar_vars = domains[scalars];
    for(var=scalar_vars.begin(); var!=scalar_vars.end(); var++)
    {
      const std::string& name = var->first;
      TraceSample<T>& sample = *tsamples[var->second];
      T value = buffer(current_row,sample.buffer_start);
      T svalue = 0;
      bool any_present = false;
      for(dom=domains.begin(); dom!=domains.end(); dom++)
      {
        const std::string& domain = dom->first;
        if(domain!=scalars)
        {
          std::map<std::string,int>& vars = dom->second;
          if(vars.count(name)>0)
          {
            any_present = true;
            TraceSample<T>& ssample = *tsamples[vars[name]];
            int start = ssample.buffer_start;
            int end   = ssample.buffer_end;
            for(int i=start; i<end; i++)
              svalue += buffer(current_row,i);
          }
        }
      }
      if(any_present)
      {
        if(verbose)
          app_log()<<"  "<<name<<" "<<value<<" "<<svalue<< std::endl;
      }
    }
  }

};




class TraceManager
{
private:
  //collections of samples for a single walker step
  //  the associated arrays will be updated following evaluate
  TraceSamples<TraceInt>   int_samples;
  TraceSamples<TraceReal>  real_samples;
  TraceSamples<TraceComp>  comp_samples;

  //buffers for storing samples
  // single row of buffer is a single sample from one walker
  // number of rows adjusts to accomodate walker samples
  TraceBuffer<TraceInt>  int_buffer;
  TraceBuffer<TraceReal> real_buffer;

public:
  static double trace_tol;

  TraceRequest request;

  bool master_copy;
  std::string default_domain;
  bool method_allows_traces;
  bool streaming_traces;
  bool writing_traces;
  int  throttle;
  bool verbose;
  std::string format;
  bool hdf_format;
  bool adios_format;
  std::string file_root;
  Communicate* communicator;
  hdf_archive* hdf_file;
#ifdef HAVE_ADIOS
  ADIOS::Trace* adios_trace;
#endif
  xmlNodePtr adios_options;

  TraceManager(Communicate* comm=0)
    : hdf_file(0),verbose(false)
  {
    reset_permissions();
    master_copy    = true;
    communicator   = comm;
    throttle       = 1;
    format         = "hdf";
    default_domain = "scalars";
    request.set_scalar_domain(default_domain);
    int_buffer.set_type("int");
    real_buffer.set_type("real");
    int_buffer.set_samples( int_samples);
    real_buffer.set_samples(real_samples);
    real_buffer.set_samples(comp_samples);
#ifdef HAVE_ADIOS
    adios_trace = 0;
#endif
    adios_options = 0;
  }


  inline TraceManager* makeClone()
  {
    if(verbose)
      app_log()<<"TraceManager::makeClone "<<master_copy<< std::endl;
    if(!master_copy)
      APP_ABORT("TraceManager::makeClone  only the master copy should call this function");
    TraceManager* tm = new TraceManager();
    tm->master_copy  = false;
    tm->transfer_state_from(*this);
    tm->distribute();
    return tm;
  }


  inline void transfer_state_from(const TraceManager& tm)
  {
    method_allows_traces = tm.method_allows_traces;
    request              = tm.request;
    streaming_traces     = tm.streaming_traces;
    writing_traces       = tm.writing_traces;
    throttle             = tm.throttle;
    verbose              = tm.verbose;
    format               = tm.format;
    hdf_format           = tm.hdf_format;
    adios_format         = tm.adios_format;
    default_domain       = tm.default_domain;
#ifdef HAVE_ADIOS
    adios_trace = tm.adios_trace;
#endif
    adios_options = tm.adios_options;
  }


  inline void distribute()
  {
    int_samples.set_verbose(verbose);
    real_samples.set_verbose(verbose);
    comp_samples.set_verbose(verbose);
    int_buffer.set_verbose(verbose);
    real_buffer.set_verbose(verbose);
  }


  inline void reset_permissions()
  {
    method_allows_traces = false;
    streaming_traces     = false;
    writing_traces       = false;
    verbose              = false;
    hdf_format           = false;
    adios_format         = false;
    request.reset();
  }


  inline void put(xmlNodePtr cur,bool allow_traces,std::string series_root)
  {
    reset_permissions();
    method_allows_traces = allow_traces;
    file_root            = series_root;
    bool traces_requested = cur!=NULL;
    streaming_traces = traces_requested && method_allows_traces;
    if( streaming_traces)
    {
      if (omp_get_thread_num()==0)
      {
        app_log()<<"\n  TraceManager::put() "<<master_copy<< std::endl;
        app_log()<<"    traces requested          : "<<traces_requested<< std::endl;
        app_log()<<"    method allows traces      : "<<method_allows_traces<< std::endl;
        app_log()<<"    streaming traces          : "<<streaming_traces<< std::endl;
        app_log()<< std::endl;
      }
      //read trace attributes
      std::string writing           = "yes";
      std::string scalar            = "yes";
      std::string array             = "yes";
      std::string scalar_defaults   = "yes";
      std::string array_defaults    = "yes";
      std::string verbose_write     = "no";
      OhmmsAttributeSet attrib;
      attrib.add(writing,        "write"          );
      attrib.add(scalar,         "scalar"         );
      attrib.add(array,          "array"          );
      attrib.add(scalar_defaults,"scalar_defaults");
      attrib.add(array_defaults, "array_defaults" );
      attrib.add(format,         "format"         );
      attrib.add(throttle,       "throttle"       );
      attrib.add(verbose_write,  "verbose"        );
      attrib.add(array,          "particle"          );//legacy
      attrib.add(array_defaults, "particle_defaults" );//legacy
      attrib.put(cur);
      writing_traces           = writing         == "yes";
      bool scalars_on          = scalar          == "yes";
      bool arrays_on           = array           == "yes";
      bool use_scalar_defaults = scalar_defaults == "yes";
      bool use_array_defaults  = array_defaults  == "yes";
      verbose                  = verbose_write   == "yes";
      tolower(format);
      if(format=="hdf")
      {
        hdf_format = true;
      }
#ifdef HAVE_ADIOS
      else if(format=="adios")
      {
        adios_format = true;
      }
      else if(format=="both")
      {
        adios_format = true;
        hdf_format = true;
      }
      else
      {
        APP_ABORT("TraceManager::put "+format+" is not a valid file format for traces\n  valid options are: hdf, adios, both");
      }
#else
      else
      {
        APP_ABORT("TraceManager::put "+format+" is not a valid file format for traces\n  valid options is: hdf");
      }
#endif


      //read scalar and array elements
      //  each requests that certain traces be computed
      std::set<std::string> scalar_requests;
      std::set<std::string> array_requests;
      xmlNodePtr element = cur->children;
      while(element!=NULL)
      {
        std::string name((const char*)element->name);
        if(name=="scalar_traces")
        {
          std::string defaults = "no";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults,"defaults");
          eattrib.put(element);
          use_scalar_defaults = use_scalar_defaults && defaults=="yes";
          if(!use_scalar_defaults)
          {
            std::vector<std::string> scalar_list;
            putContent(scalar_list,element);
            scalar_requests.insert(scalar_list.begin(),scalar_list.end());
          }
        }
        else if(name=="array_traces" || name=="particle_traces")
        {
          std::string defaults = "no";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults,"defaults");
          eattrib.put(element);
          use_array_defaults = use_array_defaults && defaults=="yes";
          if(!use_array_defaults)
          {
            std::vector<std::string> array_list;
            putContent(array_list,element);
            array_requests.insert(array_list.begin(),array_list.end());
          }
        }
        else if(name=="adios_options")
        {
          adios_options = element;
        }
        else if(name!="text")
        {
          APP_ABORT("TraceManager::put "+name+" is not a valid sub-element of <trace/>\n  valid options are: scalar_traces, array_traces");
        }
        element=element->next;
      }

      writing_traces &= method_allows_traces;

      //input user quantity requests into the traces request
      request.allow_streams      = method_allows_traces;
      request.allow_writes       = writing_traces;
      request.scalars_on         = scalars_on;
      request.arrays_on          = arrays_on;
      request.stream_all_scalars = use_scalar_defaults;
      request.stream_all_arrays  = use_array_defaults;        
      request.write_all_scalars  = request.stream_all_scalars && writing_traces;
      request.write_all_arrays   = request.stream_all_arrays && writing_traces;
      request.request_scalar(scalar_requests,writing_traces);
      request.request_array(array_requests,writing_traces);

      //distribute verbosity level to buffer and sample objects
      distribute();
    }

    //streaming_traces = false;
    //writing_traces   = false;
  }


  inline void update_status()
  {
    streaming_traces = request.streaming();
    writing_traces   = request.writing();
  }


  inline void screen_writes()
  {
    int_samples.screen_writes(request);
    real_samples.screen_writes(request);
    comp_samples.screen_writes(request);
  }

  inline void initialize_traces()
  {
    if( streaming_traces)
    {
      if(verbose)
        app_log()<<"TraceManager::initialize_traces "<<master_copy<< std::endl;
      //organize trace samples and initialize buffers
      if(writing_traces)
      {
        int_buffer.order_and_resize();
        real_buffer.order_and_resize();
      }
    }
  }


  inline void finalize_traces()
  {
    if(verbose)
      app_log()<<"TraceManager::finalize_traces "<<master_copy<< std::endl;
    int_samples.finalize();
    real_samples.finalize();
    comp_samples.finalize();
  }


  //checkout functions to be used by any QMCHamiltonianBase or Estimator
  //  the array checked out should be updated during evaluate
  //  object calling checkout is responsible for deleting the new array

  // checkout integer arrays
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_int<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const std::string& domain,const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return int_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const std::string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return int_samples.checkout_array<D>(P,name,shape);
  }


  // checkout real arrays
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_real<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const std::string& domain,const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return real_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const std::string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return real_samples.checkout_array<D>(P,name,shape);
  }


  // checkout complex arrays
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_complex<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const std::string& domain,const std::string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return comp_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const std::string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return comp_samples.checkout_array<D>(P,name,shape);
  }


  //get trace functions
  //  used by estimators that require trace information
  inline TraceSample<TraceInt>* get_int_trace(const std::string& name)
  {
    return get_int_trace(default_domain,name);
  }
  inline TraceSample<TraceInt>* get_int_trace(const std::string& domain, const std::string& name)
  {
    return int_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceInt>* get_int_trace(const ParticleSet& P, const std::string& name)
  {
    return int_samples.get_trace(P.parentName(),name);
  }

  inline TraceSample<TraceReal>* get_real_trace(const std::string& name)
  {
    return get_real_trace(default_domain,name);
  }
  inline TraceSample<TraceReal>* get_real_trace(const std::string& domain, const std::string& name)
  {
    return real_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceReal>* get_real_trace(const ParticleSet& P, const std::string& name)
  {
    return real_samples.get_trace(P.parentName(),name);
  }

  inline TraceSample<TraceComp>* get_complex_trace(const std::string& name)
  {
    return get_complex_trace(default_domain,name);
  }
  inline TraceSample<TraceComp>* get_complex_trace(const std::string& domain, const std::string& name)
  {
    return comp_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceComp>* get_complex_trace(const ParticleSet& P, const std::string& name)
  {
    return comp_samples.get_trace(P.parentName(),name);
  }


  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const std::string& name)
  {
    return get_int_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const std::string& domain, const std::string& name)
  {
    return int_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const ParticleSet& P, const std::string& name)
  {
    return int_samples.get_combined_trace(P.parentName(),name);
  }

  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const std::string& name)
  {
    return get_real_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const std::string& domain, const std::string& name)
  {
    return real_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const ParticleSet& P, const std::string& name)
  {
    return real_samples.get_combined_trace(P.parentName(),name);
  }

  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const std::string& name)
  {
    return get_complex_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const std::string& domain, const std::string& name)
  {
    return comp_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const ParticleSet& P, const std::string& name)
  {
    return comp_samples.get_combined_trace(P.parentName(),name);
  }




  //create combined trace out of existing traces
  inline void make_combined_trace(const std::string& name,std::vector<std::string>& names)
  {
    std::vector<TraceReal> weights;
    weights.resize(names.size());
    fill(weights.begin(),weights.end(),1.0);
    make_combined_trace(name,names,weights);
  }

  inline void make_combined_trace(const std::string& name,std::vector<std::string>& names,std::vector<TraceReal>& weights)
  {
    if( streaming_traces)
    {
      if(verbose)
        app_log()<<"TraceManager::make_combined_trace "<<master_copy<<"  "<<name<< std::endl;
      real_buffer.make_combined_trace(name,names,weights);
    }
  }


  inline void check_clones(std::vector<TraceManager*>& clones)
  {
    if(writing_traces && clones.size()>0)
    {
      if(verbose)
        app_log()<<"TraceManager::check_clones"<< std::endl;
      bool all_same=true;
      bool int_same;
      bool real_same;
      TraceManager& ref = *clones[0];
      for(int i=0; i<clones.size(); ++i)
      {
        TraceManager& tm = *clones[i];
        int_same = tm.int_buffer.same_as(ref.int_buffer);
        real_same = tm.real_buffer.same_as(ref.real_buffer);
        all_same = all_same && int_same && real_same;
      }
      if(!all_same)
      {
        for(int i=0; i<clones.size(); ++i)
          clones[i]->write_summary();
        APP_ABORT("TraceManager::check_clones  trace buffer widths of clones do not match\n  contiguous write is impossible\n  this was first caused by clones contributing array traces from identical, but differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the TraceManager summaries printed above");
      }
    }
  }


  inline void reset_buffers()
  {
    if(writing_traces)
    {
      if(verbose)
        app_log()<<"TraceManager::reset_buffers "<<master_copy<< std::endl;
      int_buffer.reset();
      real_buffer.reset();
    }
  }


  //store the full sample from a single walker step in buffers
  inline void buffer_sample(int current_step)
  {
    if(writing_traces && current_step%throttle==0)
    {
      if(verbose)
        app_log()<<" TraceManager::buffer_sample() "<<master_copy<< std::endl;
      int_buffer.collect_sample();
      real_buffer.collect_sample();
    }
  }


  //write buffered trace data to file
  inline void write_buffers(std::vector<TraceManager*>& clones, int block)
  {
    if(master_copy)
    {
      if(writing_traces)
      {
        //double tstart = MPI_Wtime();
        if(verbose)
          app_log()<<"TraceManager::write_buffers "<<master_copy<< std::endl;
        if(hdf_format)
        {
          write_buffers_hdf(clones);
        }
        if(adios_format)
        {
#ifdef HAVE_ADIOS
          write_buffers_adios(clones, block);
#else
          APP_ABORT("TraceManager::write_buffers (adios) ADIOS is not found");
#endif
          //app_log()<<"TraceManager::write_buffers (adios) has not yet been implemented"<< std::endl;
          //app_log()<<" write_buffers() total time "<<MPI_Wtime()-tstart<< std::endl;
        }
      }
    }
    else
      APP_ABORT("TraceManager::write_buffers should not be called from non-master copy");
  }


  inline void open_file(std::vector<TraceManager*>& clones)
  {
    if(master_copy)
    {
      if(writing_traces)
      {
        if(verbose)
          app_log()<<"TraceManager::open_file "<<master_copy<< std::endl;
        if(verbose)
          clones[0]->write_summary();
        if(hdf_format)
        {
          open_hdf_file(clones);
        }
        if(adios_format)
        {
#ifdef HAVE_ADIOS
          initialize_adios(clones);
#else
          APP_ABORT("TraceManager::open_file (adios) ADIOS is not found");
#endif
        }
      }
    }
    else
      APP_ABORT("TraceManager::open_file should not be called from non-master copy");
  }


  inline void close_file()
  {
    if(master_copy)
    {
      if(writing_traces)
      {
        if(verbose)
          app_log()<<"TraceManager::close_file "<<master_copy<< std::endl;
        if(hdf_format)
        {
          close_hdf_file();
        }
        if(adios_format)
        {
#ifdef HAVE_ADIOS
          finalize_adios();
#else
          APP_ABORT("TraceManager::close_file (adios) ADIOS is not found");
#endif
        }
      }
    }
    else
      APP_ABORT("TraceManager::close_file should not be called from non-master copy");
  }



  inline void startRun(int blocks,std::vector<TraceManager*>& clones)
  {
    if(verbose)
      app_log()<<"TraceManager::startRun "<<master_copy<< std::endl;
    if(master_copy)
    {
      initialize_traces();
      check_clones(clones);
      open_file(clones);
    }
    else
      APP_ABORT("TraceManager::startRun should not be called from non-master copy");
  }


  inline void stopRun()
  {
    if(verbose)
      app_log()<<"TraceManager::stopRun "<<master_copy<< std::endl;
    if(master_copy)
    {
      close_file();
      finalize_traces();
    }
    else
      APP_ABORT("TraceManager::stopRun should not be called from non-master copy");
  }


  inline void startBlock(int nsteps)
  {
    if(verbose)
      app_log()<<"TraceManager::startBlock "<<master_copy<< std::endl;
    reset_buffers();
  }


  inline void stopBlock()
  {
    if(verbose)
      app_log()<<"TraceManager::stopBlock "<<master_copy<< std::endl;
  }


  inline void write_summary( std::string pad="  ")
  {
    std::string pad2 = pad+"  ";
    app_log()<< std::endl;
    app_log()<<pad<<"TraceManager (detailed summary)"<< std::endl;
    app_log()<<pad2<<"master_copy             = "<<master_copy              << std::endl;
    app_log()<<pad2<<"method_allows_traces    = "<<method_allows_traces     << std::endl;
    app_log()<<pad2<<"streaming_traces        = "<<streaming_traces         << std::endl;
    app_log()<<pad2<<"writing_traces          = "<<writing_traces           << std::endl;
    app_log()<<pad2<<"format                  = "<<format                   << std::endl;
    app_log()<<pad2<<"hdf format              = "<<hdf_format               << std::endl;
    app_log()<<pad2<<"adios format            = "<<adios_format             << std::endl;
    app_log()<<pad2<<"default_domain          = "<<default_domain           << std::endl;
    int_buffer.write_summary(pad2);
    real_buffer.write_summary(pad2);
    app_log()<<pad<<"end TraceManager"<< std::endl;
  }


  inline void user_report( std::string pad="  ")
  {
    std::string pad2 = pad+"  ";
    std::string pad3 = pad2+"  ";
    app_log()<< std::endl;
    app_log()<<pad<<"Traces report"<< std::endl;
    request.report();
    app_log()<<pad2<<"Type and domain breakdown of streaming quantities:"<< std::endl;
    std::set<std::string>::iterator req;
    int_buffer.user_report(pad3);
    real_buffer.user_report(pad3);
    app_log()<< std::endl;
    //if(verbose)
    //  write_summary(pad);
  }

  //hdf file operations
  inline void open_hdf_file(std::vector<TraceManager*>& clones)
  {
    if(clones.size()==0)
      APP_ABORT("TraceManager::open_hdf_file  no trace clones exist, cannot open file");
    int nprocs = communicator->size();
    int rank = communicator->rank();
    char ptoken[32];
    std::string file_name = file_root;
    if(nprocs>1)
    {
      if(nprocs>10000)
        sprintf(ptoken,".p%05d",rank);
      else if(nprocs>1000)
        sprintf(ptoken,".p%04d",rank);
      else
        sprintf(ptoken,".p%03d",rank);
      file_name += ptoken;
    }
    file_name += ".traces.h5";
    if(verbose)
      app_log()<<"TraceManager::open_hdf_file  opening traces hdf file "<<file_name<< std::endl;
    hdf_file = new hdf_archive(communicator,false);
    bool successful = hdf_file->create(file_name);
    if(!successful)
      APP_ABORT("TraceManager::open_hdf_file  failed to open hdf file "+file_name);
    // only clones have active buffers and associated data
    TraceManager& tm = *clones[0];
    //tm.write_summary();
    tm.int_buffer.register_hdf_data(*hdf_file);
    tm.real_buffer.register_hdf_data(*hdf_file);
  }


  inline void write_buffers_hdf(std::vector<TraceManager*>& clones)
  {
    if(verbose)
      app_log()<<"TraceManager::write_buffers_hdf "<<master_copy<< std::endl;
    for(int ip=0; ip<clones.size(); ++ip)
    {
      TraceManager& tm = *clones[ip];
      tm.int_buffer.write_hdf(*hdf_file,int_buffer.hdf_file_pointer);
      tm.real_buffer.write_hdf(*hdf_file,real_buffer.hdf_file_pointer);
    }
  }


  inline void close_hdf_file()
  {
    delete hdf_file;
  }


#ifdef HAVE_ADIOS /* Norbert's new code */

    template<typename TraceT, typename ADIOST>
    void
    convert_and_write_adios (ADIOS::Trace& at, std::vector<TraceManager*>& clones, TraceSample<TraceT> s, 
                             std::string type, int max_nrows, int ncols)
    {
        int dims[4] = {0,0,0,0};
        std::string varpath = "/"+s.domain+"/"+s.name;
        int nelems = 1;
        for (int j=0; j < s.dimension; j++)
        {
            dims[j] = s.shape[j]; 
            nelems *= dims[j];
        }

        int bs = s.buffer_start; // start offset of this variable in int samples 
        //app_log()<<"      variable "<<varpath<<" nelems="<<nelems<<" buffer offset="<<bs<< std::endl;

        ADIOST *datap;
        if (type == "complex") {
            /* nelems is the total number of COMPLEX values 
             * i.e. 2x as many real values are in the REAL buffer
             * ADIOST is a real type in this case, need to double the allocation here */
            datap = (ADIOST *) calloc (2 * nelems * max_nrows, sizeof(ADIOST));
        } else {
            datap = (ADIOST *) calloc (nelems * max_nrows, sizeof(ADIOST));
        }

        int row_start = 0; // cumulative row count
        for(int i=0; i<clones.size(); i++)
        {
            TraceManager& tmc = *clones[i];
            int nrows;
            if (type == "int") {
                Array<TraceInt,2>& buffer = tmc.int_buffer.buffer;
                nrows = buffer.size(0);
                for(int j=0; j<buffer.size(0); j++){
                    for(int k=0; k<nelems; k++){
                        // trace[j,k] -> adios[k, j+row_start] 
                        datap[(row_start+j)*nelems+k] = buffer(j*ncols+bs+k);
                    }
                }
            } else if (type == "complex") {
                // The TraceComp type sample is stored in the TraceReal typed array
                Array<TraceReal,2>& buffer = tmc.real_buffer.buffer;
                nrows = buffer.size(0);
                for(int j=0; j<buffer.size(0); j++){
                    for(int k=0; k<nelems; k++){
                        // trace[j,k] -> adios[k, j+row_start] 
                        datap[(row_start+j)*2*nelems+2*k] = buffer(j*ncols+bs+2*k);
                        datap[(row_start+j)*2*nelems+2*k+1] = buffer(j*ncols+bs+2*k+1);
                    }
                }
            } else { 
                Array<TraceReal,2>& buffer = tmc.real_buffer.buffer;
                nrows = buffer.size(0);
                for(int j=0; j<buffer.size(0); j++){
                    for(int k=0; k<nelems; k++){
                        // trace[j,k] -> adios[k, j+row_start] 
                        datap[(row_start+j)*nelems+k] = buffer(j*ncols+bs+k);
                    }
                }
            }

            row_start += nrows;
        }
        adios_trace->write (varpath, datap);
        free (datap);
    }

    /* Called once at each QMC application run (i.e. a few times during the run)
     * Here we define an ADIOS group of variables for output
     * Name of the group is same as file_root, e.g. bench.s000
     */
    inline void initialize_adios(std::vector<TraceManager*>& clones)
    {
        TraceManager& tm = *clones[0];

        //double tstart = MPI_Wtime();
        if(verbose)
            app_log()<<"TraceManager::initialize_adios "<<master_copy<< std::endl;

        adios_trace = new ADIOS::Trace (file_root.c_str(), communicator->getMPI(), adios_options);

        tm.int_buffer.register_adios_data(*adios_trace);
        tm.real_buffer.register_adios_data(*adios_trace);
        //app_log()<<"TraceManager::initialize_adios() total time="<<MPI_Wtime()-tstart<< std::endl;
    }


    inline void write_buffers_adios(std::vector<TraceManager*>& clones, int block)
    {
        if(verbose)
            app_log()<<"TraceManager::write_buffers_adios "<<master_copy<< std::endl;
        int64_t f;
        MPI_Comm comm = communicator->getMPI();

        /* I. open the file, so if that's threaded, we can prepare the 
         * buffering concurrently */
        std::string s = boost::lexical_cast<std::string>(block);
        std::string file_name = file_root+".b"+s+".bp";
        //double io_open_start = MPI_Wtime();
        adios_trace->open (file_name);
        //app_log()<<" io open time on rank 0: "<<MPI_Wtime()-io_open_start<< std::endl;

        int actual_nrows_local = 0; // number of rows in this process' tables
        int max_nrows = 0;          // max of local nrows numbers
        int nrows_total;            // number of rows in the global output
        int offset;                 // this process' offset in global output of samples

        /* II. Get the maximum number of rows present in any process.
         * Each process will pad it's output arrays with 0 to have the same size
         */
        for(int ip=0; ip<clones.size(); ip++){
            TraceManager& tmc = *clones[ip];
            actual_nrows_local += tmc.real_buffer.buffer.size(0);
            if (tmc.int_buffer.buffer.size(0) != tmc.real_buffer.buffer.size(0))
            {
                std::cerr <<"ERROR in ADIOS trace assumptions: number of rows in the integer" 
                    "vs real trace tables do not match on process "
                    <<communicator->rank()
                    <<": clone="<<ip<<" int_nrows="<<tmc.int_buffer.buffer.size(0)
                    <<" real_nrows="<<tmc.real_buffer.buffer.size(0)<< std::endl;
            }
        }
        //double comm_start = MPI_Wtime();
        MPI_Allreduce(&actual_nrows_local, &max_nrows, 1, MPI_INT, MPI_MAX, comm);
        //app_log()<<" Allreduce comm time for nrows "<<MPI_Wtime()-comm_start<< std::endl;

        nrows_total = max_nrows * communicator->size();
        offset      = max_nrows * communicator->rank();

        /* III. Calculate adios group size */
        //double io_group_start = MPI_Wtime();
        adios_trace->set_group_size (max_nrows);
        //app_log()<<" io group time on rank 0: "<<MPI_Wtime()-io_group_start<< std::endl;

        /* IV. Write known adios scalars */
        //double io_write_start = MPI_Wtime();
        adios_trace->write("/aux/max_nrows", &max_nrows);
        adios_trace->write("/aux/actual_nrows_local", &actual_nrows_local);
        adios_trace->write("/aux/offset", &offset);
        adios_trace->write("nrows", &nrows_total);

        /* V. prepare and write output variables 
         * We need to a buffering each entity separately here to have contiguous 
         * pieces for each trace variable. 
         */
        TraceManager& tm = *clones[0];
        int int_cols = tm.int_buffer.buffer.size(1); // number of columns in int table
        int real_cols = tm.real_buffer.buffer.size(1); // number of columns in real table

        /* 
         * process INTEGER data 
         */ 
        //app_log()<<"    Copy integer buffer to ADIOS variables..."<< std::endl;
        struct TraceSamples<TraceInt>& int_samples = *tm.int_buffer.samples;
        for(int i=0; i<int_samples.ordered_samples.size(); ++i) 
        {
            const struct TraceSample<TraceInt>& s = *int_samples.ordered_samples[i];
            //jtk: attempt to restrict adios writes in the correct way, please check!
            if(s.write)
            {
              convert_and_write_adios<TraceInt,int> (
                      *adios_trace, clones, s, "int", max_nrows, int_cols);
            }
        }

        /* 
         * process REAL data 
         */ 
#define ADIOSTraceReal float
//#define ADIOSTraceReal TraceReal
        //app_log()<<"    Copy real buffer to ADIOS variables..."<< std::endl;
        struct TraceSamples<TraceReal>& real_samples = *tm.real_buffer.samples;
        for(int i=0; i<real_samples.ordered_samples.size(); ++i) 
        {
            const struct TraceSample<TraceReal>& s = *real_samples.ordered_samples[i];
            //jtk: attempt to restrict adios writes in the correct way, please check!
            if(s.write)
            {
              std::string vartype = real_samples.get_adios_type ("real", s.domain, s.name, false);
              // we may write some variables as doubles (/scalars/*)
              if (vartype == "double") {
                convert_and_write_adios<TraceReal,double> (
                        *adios_trace, clones, s, "real", max_nrows, real_cols);
              } else {
                convert_and_write_adios<TraceReal,ADIOSTraceReal> (
                        *adios_trace, clones, s, "real", max_nrows, real_cols);
              }
            }
        }


        /* 
         * process COMPLEX data 
         */ 
        //app_log()<<"    Copy complex buffer to ADIOS variables..."<< std::endl;
        struct TraceSamples<TraceComp>& complex_samples = *tm.real_buffer.complex_samples;
        for(int i=0; i<complex_samples.ordered_samples.size(); ++i) 
        {
            const struct TraceSample<TraceComp>& s = *complex_samples.ordered_samples[i];
            //jtk: attempt to restrict adios writes in the correct way, please check!
            if(s.write)
            {
              convert_and_write_adios<TraceComp,ADIOSTraceReal> (
                      *adios_trace, clones, s, "complex", max_nrows, real_cols);
            }
        }


        //app_log()<<" io write time on rank 0: "<<MPI_Wtime()-io_write_start<< std::endl;
        //double io_close_start = MPI_Wtime();
        adios_trace->close();
        //app_log()<<" io close time on rank 0: "<<MPI_Wtime()-io_close_start<< std::endl;

    }

    inline void finalize_adios()
    {
        if(verbose)
            app_log()<<"TraceManager::finalize_adios "<< std::endl;
        delete(adios_trace);
    }

#else
#  ifdef HAVE_ADIOS_OLD


  inline void initialize_adios(std::vector<TraceManager*>& clones)
  {
      //adios_init("qmc_adios.xml", communicator->getMPI());
      TraceManager& tm = *clones[0];
      //tm.int_buffer.register_adios_data("w");
      //tm.real_buffer.register_adios_data("a");
  }


  inline void print_adios(std::vector<TraceManager*>& clones)
  {
      app_log()<<"TraceManager::write_buffers_adios "<<master_copy<< std::endl;
      for(int ip=0; ip<clones.size(); ++ip)
      {
          TraceManager& tm = *clones[ip];
          //tm.int_buffer.write_summary();
          //tm.real_buffer.write_summary();
          app_log() << "Type: " << tm.real_buffer.type << std::endl;
          app_log() << "Buffer Size: " << tm.real_buffer.buffer.size() << " | " << tm.real_buffer.buffer.size(0) << "|" <<  tm.real_buffer.buffer.size(1) << std::endl;
      }
  }

  inline void write_buffers_adios(std::vector<TraceManager*>& clones, int block)
  {
      MPI_Comm comm = communicator->getMPI();

      static bool write_global; 
#if (defined WRITE_GLOBAL) //passed from CMakeList file to switch between local array and global array
      write_global = true;
#else
      write_global = false;
#endif

      if(write_global){
          //preparing data to be written out
          int rows[]={0,0};
          int cols[]={0,0};
          for(int ip=0; ip<clones.size(); ip++){
              TraceManager& tm = *clones[ip];
              rows[0] += tm.int_buffer.buffer.size(0);
              rows[1] += tm.real_buffer.buffer.size(0);
              cols[0] = tm.int_buffer.buffer.size(1);
              cols[1] = tm.real_buffer.buffer.size(1);
          }
          int allreduced_rows[2];
          int total_rows[2];
          //double comm_start = MPI_Wtime();
          MPI_Allreduce(rows, allreduced_rows, 2, MPI_INT, MPI_MAX, comm);
          //app_log()<<" comm time "<<MPI_Wtime()-comm_start<< std::endl;
          //double prep_start = MPI_Wtime();

          int int_max_rows = allreduced_rows[0];
          int int_rows_total = int_max_rows*communicator->size();
          int int_offset_rows = communicator->rank()*int_max_rows;
          int int_cols = cols[0];

          int real_max_rows = allreduced_rows[1];
          int real_rows_total = real_max_rows*communicator->size();
          int real_offset_rows = communicator->rank()*real_max_rows;
          int real_cols = cols[1];

          int* int_buffer = (int *)malloc(int_max_rows*int_cols*sizeof(int));
          double* real_buffer = (double*)malloc(real_max_rows*real_cols*sizeof(double));
          if(!int_buffer || !real_buffer){
              APP_ABORT("not enough memory\n");
          }
          bzero(int_buffer, int_max_rows*int_cols*sizeof(int));
          bzero(real_buffer, real_max_rows*real_cols*sizeof(double));

          int row_count = 0;
          for(int i=0; i<clones.size(); i++)
          {
              TraceManager& tm = *clones[i];
              int rows = tm.int_buffer.buffer.size(0);
              for(int j=0; j<rows; j++){
                  for(int k=0; k<int_cols; k++){
                      int_buffer[row_count*int_cols+k] = tm.int_buffer.buffer(j*int_cols+k);
                  }
                  row_count++;
              }
          }

          row_count = 0;
          for(int i=0; i<clones.size(); ++i)
          {
              TraceManager& tm = *clones[i];
              int rows = tm.real_buffer.buffer.size(0);
              for(int j=0; j<rows; j++){
                  for(int k=0; k<real_cols; k++){
                      real_buffer[row_count*real_cols+k] = tm.real_buffer.buffer(j*real_cols+k);
                  }
                  row_count++;
              }
          }
          //app_log()<<" preparation time "<<MPI_Wtime()-prep_start<< std::endl;

          //open adios file, create adios group, write out variable, close adios file
          int         err;
          uint64_t    group_size, total_size;
          int64_t     handle;
          std::string file_name = file_root;
          char *fileName = (char*)file_name.c_str();

#  if 0
          /* Current trunk */
          //string s = boost::lexical_cast<std::string>(block);
          //file_name = file_name + ".b"+s+".trace.bp";
          //double io_open_start = MPI_Wtime();
          if(ADIOS::getFirstOpen())
          {
              if(!ADIOS::get_adios_init()){
                  APP_ABORT("adios is not properly initialized, check adiosinit in XML");
              }
              adios_open(&handle, "Traces-global", ADIOS::getTraceFileName().c_str(), "w", comm);
              ADIOS::setFirstOpen(false);
          } 
          else 
          {
              adios_open(&handle, "Traces-global", ADIOS::getTraceFileName().c_str(), "a", comm);
          }
#  else
          /* Norbert */
          std::string s = boost::lexical_cast<std::string>(block);
          file_name = "trace."+s+".bp";
          //double io_open_start = MPI_Wtime();
          adios_open(&handle, "Traces-global", file_name.c_str(), "w", comm);
#  endif
          //app_log()<<" io open time "<<MPI_Wtime()-io_open_start<< std::endl;

          //double io_group_start = MPI_Wtime();
          group_size = strlen(fileName) + 8*sizeof(int) + int_max_rows*int_cols*sizeof(int)+real_max_rows*real_cols*sizeof(double);
          adios_group_size (handle, group_size, &total_size);
          //app_log()<<" io group time "<<MPI_Wtime()-io_group_start<< std::endl;

          //double io_write_start = MPI_Wtime();
          adios_write(handle, "filename", fileName);
          adios_write(handle, "int_rows_total", &int_rows_total);
          adios_write(handle, "int_max_rows", &int_max_rows);
          adios_write(handle, "int_offset_rows", &int_offset_rows);
          adios_write(handle, "int_cols", &int_cols);
          adios_write(handle, "real_rows_total", &real_rows_total);
          adios_write(handle, "real_max_rows", &real_max_rows);
          adios_write(handle, "real_offset_rows", &real_offset_rows);
          adios_write(handle, "real_cols", &real_cols);
          adios_write(handle, "int_buffer", int_buffer);
          adios_write(handle, "real_buffer", real_buffer);
          //app_log()<<" io write time "<<MPI_Wtime()-io_write_start<< std::endl;
          //double io_close_start = MPI_Wtime();
          adios_close(handle);
          //app_log()<<" io close time "<<MPI_Wtime()-io_close_start<< std::endl;


          free(int_buffer);
          free(real_buffer);
#ifdef IO_PROFILE
          ADIOS_PROFILE::profile_adios_size(communicator, ADIOS_PROFILE::TRACES, group_size, total_size);
#endif
#ifdef ADIOS_VERIFY
#endif
          //cout<<"write global"<<rows[0]<<" "<<rows[1]<< std::endl;
          //cout<<"allreduced"<<int_max_rows<<" "<<real_max_rows<<" "<<int_rows_total<<" "<<real_rows_total<<" "<<int_offset_rows<<" "<<real_offset_rows<<" "<<int_cols<<" "<<real_cols<< std::endl;
          //cout<<"###############################"<<block<< std::endl;
      } else {
          //write in local arrays. Can use the aggregate method
          int total_size = 0;
          for(int ip=0; ip<clones.size(); ++ip)
          {
              TraceManager& tm = *clones[ip];
              total_size += tm.real_buffer.buffer.size();
          }
          double * adios_buffer = (double *)malloc(total_size*sizeof(double));
          if(!adios_buffer){
              APP_ABORT("not enough memory\n");
          }
          int curr_pos = 0;
          for(int ip=0; ip<clones.size(); ++ip)
          {
              TraceManager& tm = *clones[ip];
              for(std::vector<TraceReal>::iterator  iter= tm.real_buffer.buffer.begin(); iter != tm.real_buffer.buffer.end(); iter++)
              {
                  adios_buffer[curr_pos++] = *iter;
              }
          }
          int         adios_err;
          uint64_t    adios_groupsize, adios_totalsize;
          int64_t     adios_handle;
          std::string file_name = file_root;
          std::string s = boost::lexical_cast<std::string>(block);
          file_name = file_name + ".b"+s+".trace.bp";
          adios_open(&adios_handle, "Traces", file_name.c_str(), "w", comm);
          adios_groupsize = 4 + (total_size * 8);
          adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
          adios_write(adios_handle, "total_size", &total_size);
          adios_write(adios_handle, "buffer_contents", adios_buffer);
          adios_close(adios_handle);
          free(adios_buffer);
#ifdef IO_PROFILE
          ADIOS_PROFILE::profile_adios_size(communicator, ADIOS_PROFILE::TRACES, adios_groupsize, adios_totalsize);
#endif
#ifdef ADIOS_VERIFY
          ADIOS_FILE *fp = adios_read_open_file(file_name.c_str(),
                  ADIOS_READ_METHOD_BP,
                  OHMMS::Controller->getMPI());
          IO_VERIFY::adios_checkpoint_verify_variables(fp, "total_size", &total_size);
          IO_VERIFY::adios_trace_verify_local_variables(fp, "buffer_contents", adios_buffer);
          adios_read_close(fp);
#endif
      }
  }


  inline void finalize_adios()
  {
      //adios_finalize(communicator->rank());
  }

#  endif /* HAVE_ADIOS_OLD */
#endif   /* HAVE_ADIOS */


};


}



#else


// make a vacuous class for TraceManager for lighter compilation
//   disabling TraceManager should not affect other runtime behavior

#include <Particle/ParticleSet.h>

namespace qmcplusplus
{

typedef long   TraceInt;
typedef double TraceReal;
typedef std::complex<TraceReal> TraceComp;

struct TraceRequest
{
  bool streaming_default_scalars;

  inline TraceRequest()
  {
    streaming_default_scalars = false;
  }

  //inline bool screen_sample(const std::string& domain,const std::string& name,bool& write) { return false; }
  inline bool streaming_scalar(const std::string& name)                               { return false; }
  inline bool streaming_array(const std::string& name)                                { return false; }
  inline bool streaming(const std::string& name)                                      { return false; }
  //inline bool quantity_present(const std::string& name)                               { return false; }
  inline bool streaming()                                                        { return false; }
  //inline bool writing()                                                          { return false; }
  //inline bool streaming_scalars()                                                { return false; }
  //inline bool streaming_arrays()                                                 { return false; }      


  //inline void set_scalar_domain(const std::string& domain)                           { }
  inline void request_scalar(const std::string& name,bool write=false)               { }
  inline void request_array(const std::string& name,bool write=false)                { }
  //inline void request_scalar(const std::set<std::string>& names,bool write=false)         { }
  //inline void request_array(const std::set<std::string>& names,bool write=false)          { }
  inline void incorporate(TraceRequest& other)                                  { }
  inline void determine_stream_write()                                          { }
  inline void relay_stream_info(TraceRequest& other)                            { }
  //inline void report()                                                          { }
  //inline void write_selected(const std::string& header,const std::string& selector)       { }
  //inline void guarantee_presence(const std::string& name,bool combined=false)        { }
  //inline void check_presence(const std::string& name)                                { }
  inline void contribute_scalar(const std::string& name,bool default_quantity=false) { }
  inline void contribute_array(const std::string& name,bool default_quantity=false)  { }
  inline void contribute_combined(const std::string& name,std::vector<std::string>&deps,
                                  bool scalar=false,bool array=false,
                                  bool default_quantity=false)                  { }



};



template<typename T>
struct TraceSample
{
  std::vector<T>&  sample;
};


template<typename T>
struct CombinedTraceSample : public TraceSample<T>
{
  inline void combine() { }
};


struct TraceManager
{
  TraceRequest request;
  bool streaming_traces;


  TraceManager(Communicate* comm=0)
  {
    streaming_traces = false;
  }

  inline TraceManager* makeClone() { return new TraceManager(); }

  inline void transfer_state_from(const TraceManager& tm)              { }    
  //inline void distribute()                                             { }
  //inline void reset_permissions()                                      { }
  inline void put(xmlNodePtr cur,bool allow_traces,std::string series_root) { }
  inline void update_status()                                          { }
  inline void screen_writes()                                          { }
  inline void initialize_traces()                                      { }
  inline void finalize_traces()                                        { }  

  template<int D> inline Array<TraceInt,D>*  checkout_int(    const std::string& name,                       int n1=1,int n2=0,int n3=0,int n4=0) { return 0; }
  //template<int D> inline Array<TraceInt,D>*  checkout_int(    const std::string& domain,const std::string& name,  int n1=1,int n2=0,int n3=0,int n4=0) { return 0; }
  template<int D> inline Array<TraceInt,D>*  checkout_int(    const std::string& name,  const ParticleSet& P,         int n2=0,int n3=0,int n4=0) { return 0; }
  template<int D> inline Array<TraceReal,D>* checkout_real(   const std::string& name,                       int n1=1,int n2=0,int n3=0,int n4=0) { return 0; }
  //template<int D> inline Array<TraceReal,D>* checkout_real(   const std::string& domain,const std::string& name,  int n1=1,int n2=0,int n3=0,int n4=0) { return 0; } 
  template<int D> inline Array<TraceReal,D>* checkout_real(   const std::string& name,  const ParticleSet& P,         int n2=0,int n3=0,int n4=0) { return 0; }
  //template<int D> inline Array<TraceComp,D>* checkout_complex(const std::string& name,                       int n1=1,int n2=0,int n3=0,int n4=0) { return 0; }
  //template<int D> inline Array<TraceComp,D>* checkout_complex(const std::string& domain,const std::string& name,  int n1=1,int n2=0,int n3=0,int n4=0) { return 0; } 
  template<int D> inline Array<TraceComp,D>* checkout_complex(const std::string& name,  const ParticleSet& P,         int n2=0,int n3=0,int n4=0) { return 0; }

  //inline TraceSample<TraceInt>* get_int_trace(const std::string& name)                            { return 0; }    
  //inline TraceSample<TraceInt>* get_int_trace(const std::string& domain, const std::string& name)      { return 0; }
  //inline TraceSample<TraceInt>* get_int_trace(const ParticleSet& P, const std::string& name)      { return 0; }
  inline TraceSample<TraceReal>* get_real_trace(const std::string& name)                          { return 0; }
  //inline TraceSample<TraceReal>* get_real_trace(const std::string& domain, const std::string& name)    { return 0; }
  inline TraceSample<TraceReal>* get_real_trace(const ParticleSet& P, const std::string& name)    { return 0; }
  //inline TraceSample<TraceComp>* get_complex_trace(const std::string& name)                       { return 0; }
  //inline TraceSample<TraceComp>* get_complex_trace(const std::string& domain, const std::string& name) { return 0; }
  inline TraceSample<TraceComp>* get_complex_trace(const ParticleSet& P, const std::string& name) { return 0; } 

  //inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const std::string& name)                            { return 0; }  
  //inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const std::string& domain, const std::string& name)      { return 0; }
  //inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const ParticleSet& P, const std::string& name)      { return 0; }
  //inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const std::string& name)                          { return 0; }
  //inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const std::string& domain, const std::string& name)    { return 0; }
  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const ParticleSet& P, const std::string& name)    { return 0; }
  //inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const std::string& name)                       { return 0; }
  //inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const std::string& domain, const std::string& name) { return 0; }
  //inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const ParticleSet& P, const std::string& name) { return 0; }

  inline void make_combined_trace(const std::string& name,std::vector<std::string>& names)                            { }     
  //inline void make_combined_trace(const std::string& name,std::vector<std::string>& names,std::vector<TraceReal>& weights) { }
  //inline void check_clones(std::vector<TraceManager*>& clones)                                              { }
  //inline void reset_buffers()                                                                          { }
  inline void buffer_sample(int current_step)                                                          { }
  inline void write_buffers(std::vector<TraceManager*>& clones, int block)                                  { }
  //inline void open_file(std::vector<TraceManager*>& clones)                                                 { }
  //inline void close_file()                                                                             { }
  inline void startRun(int blocks,std::vector<TraceManager*>& clones)                                       { }
  inline void stopRun()                                                                                { }
  inline void startBlock(int nsteps)                                                                   { }
  inline void stopBlock()                                                                              { }
  //inline void write_summary( std::string pad="  ")                                                           { }
  inline void user_report( std::string pad="  ")                                                             { }     

};


}

#endif






#endif






