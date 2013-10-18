//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_TRACEMANAGER_H
#define QMCPLUSPLUS_TRACEMANAGER_H



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
#include "adios.h"
#include "adios_read.h"
#include "ADIOS/ADIOS_config.h"
#ifdef IO_PROFILE
#include "ADIOS/ADIOS_profile.h"
#endif
#ifdef ADIOS_VERIFY
#include "ADIOS/ADIOS_verify.h"
#endif
#endif


namespace qmcplusplus
{

//#define TRACE_CHECK

const unsigned int DMAX=4;
typedef long   TraceInt;
typedef double TraceReal;
typedef complex<TraceReal> TraceComp;


struct TraceRequest
{
  bool any;
  bool scalars;
  bool particles;

  inline TraceRequest(bool a=false,bool s=false, bool p=false)
  {
    any = a;
    scalars = s;
    particles = p;
  }

  inline void clear()
  {
    any = false;
    scalars = false;
    particles = false;
  }

  inline void update()
  {
    any = scalars || particles;
  }


  inline bool operator!=(TraceRequest& other)
  {
    this->update();
    other.update();
    return any!=other.any || scalars!=other.scalars || particles!=other.particles;
  }

  inline bool operator==(TraceRequest& other)
  {
    return !(*this!=other);
  }
};




template<typename T>
struct TraceSample
{
  string domain;
  string name;
  int index;
  bool particle_trace;
  int dimension;
  int size;
  int unit_size;
  int data_size;
  TinyVector<int,DMAX> shape;
  vector<T>&  sample;
  int buffer_start,buffer_end;
  map<string,TraceInt>  meta_int;
  map<string,TraceReal> meta_real;
  map<string,string> meta_string;
  bool verbose;

  inline TraceSample(const string& sdomain,const string& sname,int sindex,int sdim,vector<T>& ssample)
    : sample(ssample), verbose(false)
  {
    initialize(sdomain,sname,sindex,sdim);
  }


  inline TraceSample(const string& sdomain,const string& sname,int sindex,int sdim,TinyVector<int,DMAX> sshape,vector<T>& ssample)
    : sample(ssample), verbose(false)
  {
    initialize(sdomain,sname,sindex,sdim);
    shape = sshape;
    size  = sample.size();
    check_shape();
  }


  inline void initialize(const string& sdomain,const string& sname,int sindex,int sdim)
  {
    domain    = sdomain,
    name      = sname;
    dimension = sdim;
    index     = sindex;
    particle_trace = false;
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


  inline void set_buffer_range(int& bstart)
  {
    set_data_size();
    buffer_start = bstart;
    buffer_end   = bstart + data_size;
    bstart = buffer_end;
  }


  inline T sum()
  {
    T s(0);
    for(int i=0; i<sample.size(); ++i)
      s+=sample[i];
    return s;
  }

  inline void write_summary(int ind=-1,string pad="  ")
  {
    string pad2 = pad+"  ";
    string pad3 = pad+"  ";
    if(ind==-1)
      app_log()<<pad<<" TraceSample "<<name<<endl;
    else
      app_log()<<pad<<ind<<" TraceSample "<<name<<endl;
    app_log()<<pad2<<"domain         = "<< domain        <<endl;
    app_log()<<pad2<<"name           = "<< name          <<endl;
    app_log()<<pad2<<"index          = "<< index         <<endl;
    app_log()<<pad2<<"particle_trace = "<< particle_trace<<endl;
    app_log()<<pad2<<"dimension      = "<< dimension     <<endl;
    app_log()<<pad2<<"size           = "<< size          <<endl;
    app_log()<<pad2<<"unit_size      = "<< unit_size     <<endl;
    app_log()<<pad2<<"data_size      = "<< data_size     <<endl;
    app_log()<<pad2<<"shape          = "<< shape         <<endl;
    app_log()<<pad2<<"buffer range   = ["<<buffer_start<<","<<buffer_end<<")"<<endl;
  }

};




template<typename T>
struct CombinedTraceSample : public TraceSample<T>
{
  bool combined;
  vector<TraceReal>          weights;
  vector<TraceSample<T>*> components;


  inline CombinedTraceSample(const string& sdomain,const string& sname,int sindex,int sdim,vector<T>& ssample)
    : TraceSample<T>(sdomain,sname,sindex,sdim,ssample)
  {
    reset();
  }


  inline CombinedTraceSample(const string& sdomain,const string& sname,int sindex,int sdim,TinyVector<int,DMAX> sshape,vector<T>& ssample)
    : TraceSample<T>(sdomain,sname,sindex,sdim,sshape,ssample)
  {
    reset();
  }


  inline void reset()
  {
    combined = false;
  }


  inline void add_component(TraceSample<T>* component,TraceReal weight)
  {
    if(components.size()==0)
    {
      this->dimension      = component->dimension;
      this->size           = component->size;
      this->shape          = component->shape;
      this->data_size      = component->data_size;
      this->particle_trace = component->particle_trace;
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
      vector<T>& component = components[i]->sample;
      for(int j=0; j<this->sample.size(); ++j)
        this->sample[j] += weight*component[j];
    }
    combined = true;
  }


  inline void write_summary_combined(int ind,string pad="  ")
  {
    string pad2 = pad+"  ";
    string pad3 = pad2+"  ";
    app_log()<<pad<<ind<<" CombinedTraceSample "<<this->name<<endl;
    app_log()<<pad2<<"domain      = "<<this->domain<<endl;
    app_log()<<pad2<<"ncomponents = "<<components.size()<<endl;
    app_log()<<pad2<<"components"<<endl;
    for(int i=0; i<components.size(); ++i)
    {
      TraceSample<T>& c = *components[i];
      app_log()<<pad3<<c.name<<" "<<c.index<<" "<<weights[i] <<endl;
    }
    app_log()<<pad2<<"end components"<<endl;
    app_log()<<pad2<<"vector address = "<<(size_t)&this->sample<<endl;
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
  vector<TraceSample<T>*>         samples;
  map<string,map<string,int> >    sample_indices;
  vector<TraceSample<T>*>         ordered_samples;
  vector<CombinedTraceSample<T>*> combined_samples;
  vector<vector<T>*>              combined_sample_vectors;
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


  inline void assign_sample_index(const string& domain,const string& name,int index,string label="")
  {
    if(sample_indices.count(domain)>0 && sample_indices[domain].count(name)>0)
    {
      //cynthia: remove for output
      APP_ABORT("TraceSamples::checkout "+label+" variable "+name+" already exists in domain "+domain);
    }
    else
    {
      sample_indices[domain][name]=index;
    }
  }

  template<int D>
  inline Array<T,D>* checkout_array(const string& domain,const string& name,TinyVector<int,DMAX> shape)
  {
    int index = samples.size();
    assign_sample_index(domain,name,index,"array");
    Array<T,D>* a = new Array<T,D>(shape.data());
    TraceSample<T>* s = new TraceSample<T>(domain,name,index,D,shape,a->storage());
    samples.push_back(s);
    if(verbose)
      app_log()<<"TraceSamples::checkout_array  "<<domain<<" "<<name<<" "<<index<<endl;
    return a;
  }


  template<int D>
  inline Array<T,D>* checkout_array(const ParticleSet& P,const string& name,TinyVector<int,DMAX> shape)
  {
    const string& domain = P.parentName();
    int index = samples.size();
    assign_sample_index(domain,name,index,"particle");
    Array<T,D>* a = new Array<T,D>(shape.data());
    TraceSample<T>* s = new TraceSample<T>(domain,name,index,D,shape,a->storage());
    samples.push_back(s);
    s->particle_trace = true;
    if(verbose)
      app_log()<<"TraceSamples::checkout_array  "<<domain<<" "<<name<<" "<<index<<endl;
    return a;
  }


  inline TraceSample<T>* get_trace(const string& domain, const string& name)
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
    return ts;
  }


  inline CombinedTraceSample<T>* get_combined_trace(const string& domain, const string& name)
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
    return ts;
  }


  inline void make_combined_trace(const string& name,vector<string>& names,vector<TraceReal>& weights)
  {
    if(names.size()!=weights.size())
      APP_ABORT("TraceSamples::make_combined_trace  names and weights must be the same size");
    map<string,map<string,int> >::iterator it;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      string domain = it->first;
      map<string,int>& indices = it->second;
      bool any_present = false;
      for(int i=0; i<names.size(); ++i)
        any_present = any_present || indices.count(names[i])>0;
      if(any_present)
      {
        int index = samples.size();
        vector<T>* sample = new vector<T>;
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
      }
    }
  }


  inline void set_unit_size(int usize)
  {
    for(int i=0; i<samples.size(); i++)
      samples[i]->set_unit_size(usize);
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
      min_index = min(min_index,samples[i]->buffer_start);
    return min_index;
  }


  inline int max_buffer_index()
  {
    int max_index = -1;
    for(int i=0; i<samples.size(); i++)
      max_index = max(max_index,samples[i]->buffer_end);
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

#ifdef HAVE_ADIOS
  inline void register_adios_data(const char* write_mode, const char* path)
  {
    MPI_Comm comm = OHMMS::Controller->getMPI();
    int adios_err;
    uint64_t adios_groupsize, adios_totalsize;
    int64_t adios_handle;
    map<string,map<string,int> >::iterator it;
    map<string,int>::iterator it2;
    int num_quantities = 0;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      map<string,int>& indices = it->second;
      num_quantities += indices.size();
    }
    adios_open(&adios_handle, "Structure", "structure.bp", write_mode, comm);
    adios_set_path(adios_handle, path);
    if (num_quantities > 0)
    {
      int dim_buf[num_quantities];
      int shape_buf[num_quantities];
      int size_buf[num_quantities];
      int unit_buf[num_quantities];
      int start_buf[num_quantities];
      int end_buf[num_quantities];
      int curr_quantity = 0;
      for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
      {
        const string& domain = it->first;
        map<string,int>& indices = it->second;
        for(it2=indices.begin(); it2!=indices.end(); ++it2)
        {
          const string& quantity = it2->first;
          const TraceSample<T>& sample = *samples[it2->second];
          dim_buf[curr_quantity] = sample.dimension;
          /* shape_buf[curr_quantity] = new int[sample.dimension]; */
          size_buf[curr_quantity] = sample.size;
          unit_buf[curr_quantity] = sample.unit_size;
          start_buf[curr_quantity] = sample.buffer_start;
          end_buf[curr_quantity] = sample.buffer_end;
          curr_quantity++;
        }
      }
      adios_groupsize = (5 * num_quantities * 4) + 4;
      adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
      //adios_write(adios_handle, "domain", (void*)(domain.c_str()));
      //adios_write(adios_handle, "quantity", (void*)(quantity.c_str()));
      /* adios_write(adios_handle, "mpi_rank", (void*)&mpi_rank); */
      /* adios_write(adios_handle, "mpi_size", (void*)&mpi_size); */
      adios_write(adios_handle, "num_quantities", &num_quantities);
      adios_write(adios_handle, "size", (void*)&size_buf);
      adios_write(adios_handle, "dimension", (void*)&dim_buf);
      //adios_write(adios_handle, "shape", (void*)&shape_buf);
      adios_write(adios_handle, "unit_size", (void*)&unit_buf);
      adios_write(adios_handle, "row_start", (void*)&start_buf);
      adios_write(adios_handle, "row_end", (void*)&end_buf);
    } //if(num_quantities > 0)
    else
    {
      adios_groupsize = 4;
      adios_group_size(adios_handle, adios_groupsize, &adios_totalsize);
      adios_write(adios_handle, "num_quantities", &num_quantities);
    }
    adios_close(adios_handle);
  }
#endif //HAVE_ADIOS


  inline void register_hdf_data(hdf_archive& f)
  {
    map<string,map<string,int> >::iterator it;
    map<string,int>::iterator it2;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      const string& domain = it->first;
      map<string,int>& indices = it->second;
      f.push(domain);
      for(it2=indices.begin(); it2!=indices.end(); ++it2)
      {
        const string& quantity = it2->first;
        TraceSample<T>& sample = *samples[it2->second];
        f.push(quantity);
        f.write(sample.dimension,   "dimension");
        f.write(sample.shape,       "shape"    );
        f.write(sample.size,        "size"     );
        f.write(sample.unit_size,   "unit_size");
        f.write(sample.buffer_start,"row_start");
        f.write(sample.buffer_end,  "row_end"  );
        f.pop();
      }
      f.pop();
    }
  }


  inline void write_summary(string type,string pad="  ")
  {
    string pad2 = pad +"  ";
    string pad3 = pad2+"  ";
    string pad4 = pad3+"  ";
    app_log()<<pad<<"TraceSamples<"<<type<<">"<<endl;
    app_log()<<pad2<<"nsamples          = "<<samples.size()<<endl;
    app_log()<<pad2<<"ncombined_samples = "<<combined_samples.size()<<endl;
    app_log()<<pad2<<"sample_indices"<<endl;
    map<string,map<string,int> >::iterator it;
    map<string,int>::iterator it2;
    for(it=sample_indices.begin(); it!=sample_indices.end(); it++)
    {
      const string& domain = it->first;
      map<string,int>& indices = it->second;
      app_log()<<pad3<<domain<<endl;
      for(it2=indices.begin(); it2!=indices.end(); ++it2)
        app_log()<<pad4<<it2->first<<" = "<<it2->second<<endl;
    }
    app_log()<<pad2<<"end sample_indices"<<endl;
    app_log()<<pad2<<"combined_sample_vectors = ";
    for(int i=0; i<combined_sample_vectors.size(); ++i)
      app_log()<<(size_t)combined_sample_vectors[i]<<" ";
    app_log()<<endl;
    app_log()<<pad2<<"combined_samples"<<endl;
    for(int i=0; i<combined_samples.size(); ++i)
      combined_samples[i]->write_summary_combined(i,pad3);
    app_log()<<pad2<<"end combined_samples"<<endl;
    app_log()<<pad2<<"samples"<<endl;
    for(int i=0; i<ordered_samples.size(); ++i)
      ordered_samples[i]->write_summary(i,pad3);
    app_log()<<pad2<<"end samples"<<endl;
    app_log()<<pad<<"end TraceSamples<"<<type<<">"<<endl;
  }
};




template<typename T>
struct TraceBuffer
{
  bool has_complex;
  TraceSamples<T>* samples;
  TraceSamples<complex<T> >* complex_samples;
  string type;
  Array<T,2> buffer;
  bool verbose;

  //hdf variables
  string top;
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


  inline void set_type(string stype)
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


  inline void set_samples(TraceSamples<complex<T> >& s)
  {
    complex_samples = &s;
    has_complex = true;
  }


  inline void make_combined_trace(const string& name,vector<string>& names,vector<TraceReal>& weights)
  {
    samples->make_combined_trace(name,names,weights);
    if(has_complex)
      complex_samples->make_combined_trace(name,names,weights);
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
      app_log()<<" TraceBuffer<"<<type<<">::collect_sample()"<<endl;
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
        app_log()<<"  increasing # of rows to "<<nrows<<endl;
      //combine samples
      samples->combine_samples();
      if(has_complex)
        complex_samples->combine_samples();
      //collect data from all samples into the buffer row
      int offset = current_row*row_size;
      {
        int boffset;
        vector<TraceSample<T>*>& ordered_samples = samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<T>& tsample = *ordered_samples[s];
          vector<T>& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0; i<sample.size(); ++i)
          {
            buffer(boffset+i) = sample[i];
          }
        }
      }
      if(has_complex)
      {
        int boffset;
        vector<TraceSample<complex<T> >*>& ordered_samples = complex_samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<complex<T> >& tsample = *ordered_samples[s];
          vector<complex<T> >& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0,ib=0; i<sample.size(); ++i,ib+=2)
          {
            buffer(boffset+ib)   = sample[i].real();
            buffer(boffset+ib+1) = sample[i].imag();
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


  inline void write_summary(string pad="  ")
  {
    string pad2=pad+"  ";
    app_log()<<pad<<"TraceBuffer<"<<type<<">"<<endl;
    app_log()<<pad2<<"nrows       = "<< buffer.size(0)<<endl;
    app_log()<<pad2<<"row_size    = "<< buffer.size(1)<<endl;
    app_log()<<pad2<<"has_complex = "<< has_complex   <<endl;
    samples->write_summary(type,pad2);
    if(has_complex)
      complex_samples->write_summary("complex "+type,pad2);
    app_log()<<pad<<"end TraceBuffer<"<<type<<">"<<endl;
  }

#ifdef HAVE_ADIOS
  inline void register_adios_data(const char* write_mode)
  {
    string path = "/" + type;
    samples->register_adios_data(write_mode, path.c_str());
    if(has_complex)
    {
      path = "/" + type + "_complex";
      complex_samples->register_adios_data("a", path.c_str());
    }
  }
#endif

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


  inline void write_hdf(hdf_archive& f)
  {
    write_hdf(f,hdf_file_pointer);
  }


  inline void write_hdf(hdf_archive& f,hsize_t& file_pointer)
  {
    if(verbose)
      app_log()<<"TraceBuffer<"<<type<<">::write_hdf() "<<file_pointer<<" "<<buffer.size(0)<<" "<<buffer.size(1)<<endl;
    dims[0] = buffer.size(0);
    dims[1] = buffer.size(1);
    if(dims[0]>0)
    {
      f.push(top);
      h5d_append(f.top(), "traces", file_pointer,
                 buffer.dim(), dims, buffer.data());
      f.pop();
    }
  }




  inline void test_buffer_write(int sample_size)
  {
    //check that the size is correct
    int ssize = samples->total_size();
    if(has_complex)
      ssize+=complex_samples->total_size();
    if(sample_size!=ssize)
    {
      app_log()<<"sample_size = "<<sample_size<<"\ntotal_size = "<<ssize<<endl;
      APP_ABORT("TraceBuffer::test_buffer_write sample_size and total_size do not match");
    }
    //check that buffer indices fall in the expected range
    int nsamples=samples->size();
    int min_index = samples->min_buffer_index();
    int max_index = samples->max_buffer_index();
    if(has_complex)
    {
      nsamples+=complex_samples->size();
      min_index = min(min_index,complex_samples->min_buffer_index());
      max_index = max(max_index,complex_samples->max_buffer_index());
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
        vector<TraceSample<T>*>& ordered_samples = samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<T>& tsample = *ordered_samples[s];
          vector<T>& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0; i<sample.size(); ++i)
            test_buffer(boffset+i) = 1;
        }
      }
      if(has_complex)
      {
        int boffset;
        vector<TraceSample<complex<T> >*>& ordered_samples = complex_samples->ordered_samples;
        for(int s=0; s<ordered_samples.size(); s++)
        {
          TraceSample<complex<T> >& tsample = *ordered_samples[s];
          vector<complex<T> >& sample = tsample.sample;
          boffset = offset + tsample.buffer_start;
          for(int i=0,ib=0; i<sample.size(); ++i,ib+=tsample.unit_size)
          {
            test_buffer(boffset+ib)   = 1;
            test_buffer(boffset+ib+1) = 1;
          }
        }
      }
      //app_log()<<"test_buffer:"<<endl;
      //for(int i=0;i<row_size;++i)
      //  app_log()<<"  "<<i<<"  "<<test_buffer(row,i)<<endl;
      for(int i=0; i<row_size; ++i)
        if(!test_buffer(row,i))
          APP_ABORT("TraceBuffer::test_buffer_write write to row is not contiguous");
    }
  }


  inline void test_buffer_collect(int current_row)
  {
    if(verbose)
      app_log()<<"TraceBuffer::test_buffer_collect"<<endl;
    string scalars = "scalars";
    map<string,map<string,int> >::iterator dom;
    map<string,int>::iterator var;
    map<string,map<string,int> >& domains = samples->sample_indices;
    vector<TraceSample<T>*>& tsamples = samples->samples;
    map<string,int>& scalar_vars = domains[scalars];
    for(var=scalar_vars.begin(); var!=scalar_vars.end(); var++)
    {
      const string& name = var->first;
      TraceSample<T>& sample = *tsamples[var->second];
      T value = buffer(current_row,sample.buffer_start);
      T svalue = 0;
      bool any_present = false;
      for(dom=domains.begin(); dom!=domains.end(); dom++)
      {
        const string& domain = dom->first;
        if(domain!=scalars)
        {
          map<string,int>& vars = dom->second;
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
        //jtk mark: replace this with a tolerance check and report
        //jtk mark: IonIon is zero
        //jtk mark LocalEnergy and LocalPotential are missing (likely from sample_indices)
        if(verbose)
          app_log()<<"  "<<name<<" "<<value<<" "<<svalue<<endl;
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

  bool master_copy;
  string default_domain;
  bool traces_requested;
  bool method_allows_traces;
  bool traces_available;
  bool scalar_traces_requested;
  bool particle_traces_requested;
  bool writing_traces;
  bool verbose;
  string format;
  bool hdf_format;
  bool adios_format;
  bool scalar_defaults_set;
  bool particle_defaults_set;
  set<string> scalar_requests;
  set<string> particle_requests;
  set<string> requests;
  TraceRequest request;
  string file_root;
  Communicate* communicator;
  hdf_archive* hdf_file;


  TraceManager(Communicate* comm=0)
    : hdf_file(0),verbose(false)
  {
    reset_permissions();
    master_copy    = true;
    communicator   = comm;
    format         = "hdf";
    default_domain = "scalars";
    int_buffer.set_type("int");
    real_buffer.set_type("real");
    int_buffer.set_samples( int_samples);
    real_buffer.set_samples(real_samples);
    real_buffer.set_samples(comp_samples);
  }


  inline TraceManager* makeClone()
  {
    if(verbose)
      app_log()<<"TraceManager::makeClone "<<master_copy<<endl;
    if(!master_copy)
      APP_ABORT("TraceManager::makeClone  only the master copy should call this function");
    TraceManager* tm = new TraceManager();
    tm->master_copy               = false;
    tm->transfer_state_from(*this);
    tm->scalar_requests.insert(scalar_requests.begin(),scalar_requests.end());
    tm->particle_requests.insert(particle_requests.begin(),particle_requests.end());
    tm->requests.insert(requests.begin(),requests.end());
    tm->distribute();
    return tm;
  }


  inline void transfer_state_from(const TraceManager& tm)
  {
    traces_requested          = tm.traces_requested;
    scalar_traces_requested   = tm.scalar_traces_requested;
    particle_traces_requested = tm.particle_traces_requested;
    method_allows_traces      = tm.method_allows_traces;
    traces_available          = tm.traces_available;
    writing_traces            = tm.writing_traces;
    verbose                   = tm.verbose;
    format                    = tm.format;
    hdf_format                = tm.hdf_format;
    adios_format              = tm.adios_format;
    default_domain            = tm.default_domain;
    scalar_defaults_set       = tm.scalar_defaults_set;
    particle_defaults_set     = tm.particle_defaults_set;
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
    method_allows_traces      = false;
    traces_requested          = false;
    scalar_traces_requested   = false;
    particle_traces_requested = false;
    traces_available          = false;
    writing_traces            = false;
    verbose                   = false;
    hdf_format                = false;
    adios_format              = false;
    scalar_defaults_set       = false;
    particle_defaults_set     = false;
  }


  inline void put(xmlNodePtr cur,bool allow_traces,string series_root)
  {
    reset_permissions();
    method_allows_traces = allow_traces;
    file_root            = series_root;
    traces_requested     = cur!=NULL;
    traces_available     = traces_requested && method_allows_traces;
    app_log()<<"\n TraceManager::put() "<<master_copy<<endl;
    app_log()<<"  traces requested          : "<<traces_requested<<endl;
    app_log()<<"  method allows traces      : "<<method_allows_traces<<endl;
    app_log()<<"  traces available          : "<<traces_available<<endl;
    bool use_scalar_defaults = false;
    bool use_particle_defaults = false;
    if(traces_available)
    {
      //read trace attributes
      string writing           = "yes";
      string scalar            = "yes";
      string particle          = "yes";
      string scalar_defaults   = "yes";
      string particle_defaults = "yes";
      string verbose_write     = "no";
      OhmmsAttributeSet attrib;
      attrib.add(writing,          "write"            );
      attrib.add(scalar,           "scalar"           );
      attrib.add(particle,         "particle"         );
      attrib.add(scalar_defaults,  "scalar_defaults"  );
      attrib.add(particle_defaults,"particle_defaults");
      attrib.add(format,           "format"           );
      attrib.add(verbose_write,    "verbose"          );
      attrib.put(cur);
      writing_traces            = writing           == "yes";
      scalar_traces_requested   = scalar            == "yes";
      particle_traces_requested = particle          == "yes";
      use_scalar_defaults       = scalar_defaults   == "yes";
      use_particle_defaults     = particle_defaults == "yes";
      verbose                   = verbose_write     == "yes";
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
      //read scalar and particle elements
      //  each requests that certain traces be computed
      xmlNodePtr element = cur->children;
      while(element!=NULL)
      {
        string name((const char*)element->name);
        if(name=="scalar_traces")
        {
          string defaults = "yes";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults,"defaults");
          eattrib.put(element);
          use_scalar_defaults = use_scalar_defaults && defaults=="yes";
          vector<string> scalar_list;
          putContent(scalar_list,element);
          scalar_requests.insert(scalar_list.begin(),scalar_list.end());
        }
        else if(name=="particle_traces")
        {
          string defaults = "yes";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults,"defaults");
          eattrib.put(element);
          use_particle_defaults = use_particle_defaults && defaults=="yes";
          vector<string> particle_list;
          putContent(particle_list,element);
          particle_requests.insert(particle_list.begin(),particle_list.end());
        }
        else if(name!="text")
        {
          APP_ABORT("TraceManager::put "+name+" is not a valid sub-element of <trace/>\n  valid options are: scalar_traces, particle_traces");
        }
        element=element->next;
      }
      app_log()<<"  writing traces            : "<< writing_traces            <<endl;
      app_log()<<"  trace output file format  : "<< format                    <<endl;
      app_log()<<"  hdf format                : "<< hdf_format                <<endl;
      app_log()<<"  adios format              : "<< adios_format              <<endl;
      app_log()<<"  scalar traces requested   : "<< scalar_traces_requested   <<endl;
      app_log()<<"  particle traces requested : "<< particle_traces_requested <<endl;
      app_log()<<"  using scalar defaults     : "<< use_scalar_defaults       <<endl;
      app_log()<<"  using particle defaults   : "<< use_particle_defaults     <<endl;
      set_default_requests(use_scalar_defaults,use_particle_defaults);
      set<string>::iterator it;
      app_log()<<"  scalar requests:"<<endl;
      if(!scalar_requests.empty())
      {
        app_log()<<"    ";
        for(it=scalar_requests.begin(); it!=scalar_requests.end(); ++it)
          app_log()<<*it<<" ";
        app_log()<<endl;
      }
      app_log()<<"  particle requests:"<<endl;
      if(!particle_requests.empty())
      {
        app_log()<<"    ";
        for(it=particle_requests.begin(); it!=particle_requests.end(); ++it)
          app_log()<<*it<<" ";
        app_log()<<endl;
      }
    }
    app_log()<<endl;
    distribute();
  }


  inline void set_default_requests(bool set_scalar, bool set_particle)
  {
    if(!scalar_defaults_set && set_scalar)
    {
      vector<string> scalar_defaults;
      scalar_defaults.push_back("LocalEnergy");
      scalar_defaults.push_back("Kinetic");
      scalar_defaults.push_back("LocalPotential");
      scalar_requests.insert(scalar_defaults.begin(),scalar_defaults.end());
      scalar_defaults_set = true;
    }
    if(!particle_defaults_set && set_particle)
    {
      vector<string> particle_defaults;
      particle_defaults.push_back("LocalEnergy");
      particle_defaults.push_back("Kinetic");
      particle_defaults.push_back("LocalPotential");
      particle_requests.insert(particle_defaults.begin(),particle_defaults.end());
      particle_defaults_set = true;
    }
  }


  inline void set_requests()
  {
    if(verbose)
      app_log()<<"TraceManager::set_requests "<<master_copy<<endl;
    requests.insert(scalar_requests.begin(),scalar_requests.end());
    requests.insert(particle_requests.begin(),particle_requests.end());
  }


  inline void add_trace_request(TraceRequest& req)
  {
    get_trace_request();
    if(req!=request)
    {
      traces_requested          = traces_requested || req.any;
      traces_available          = traces_requested && method_allows_traces;
      scalar_traces_requested   = scalar_traces_requested   || req.scalars;
      particle_traces_requested = particle_traces_requested || req.particles;
      if(traces_available)
      {
        set_default_requests(req.scalars,req.particles);
      }
    }
  }


  inline TraceRequest& get_trace_request()
  {
    if(verbose)
      app_log()<<"TraceManager::get_trace_request "<<master_copy<<endl;
    request.any       = traces_available;
    request.scalars   = scalar_traces_requested   && traces_available;
    request.particles = particle_traces_requested && traces_available;
    return request;
  }


  inline TraceRequest& get_trace_request(const string& name)
  {
    //may explicitly check if quantity has been requested later
    //  for now, just allow or disallow all based on scalars/particles from xml
    return get_trace_request();
  }


  inline void initialize_traces()
  {
    if(traces_available)
    {
      if(verbose)
        app_log()<<"TraceManager::initialize_traces "<<master_copy<<endl;
      //initialize combined traces
      //  need to add code for user-defined combinations (specified in input file)
      //check that trace dependencies of estimators are met (e.g. Energy Density)
      //  (this will be implemented later)
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
      app_log()<<"TraceManager::finalize_traces "<<master_copy<<endl;
    int_samples.finalize();
    real_samples.finalize();
    comp_samples.finalize();
  }


  //checkout functions to be used by any QMCHamiltonianBase or Estimator
  //  the array checked out should be updated during evaluate
  //  object calling checkout is responsible for deleting the new array

  // checkout integer arrays
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_int<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const string& domain,const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return int_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceInt,D>* checkout_int(const string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return int_samples.checkout_array<D>(P,name,shape);
  }


  // checkout real arrays
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_real<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const string& domain,const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return real_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceReal,D>* checkout_real(const string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return real_samples.checkout_array<D>(P,name,shape);
  }


  // checkout complex arrays
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    return checkout_complex<D>(default_domain,name,n1,n2,n3,n4);
  }
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const string& domain,const string& name,int n1=1,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(n1,n2,n3,n4);
    return comp_samples.checkout_array<D>(domain,name,shape);
  }
  template<int D>
  inline Array<TraceComp,D>* checkout_complex(const string& name,const ParticleSet& P,int n2=0,int n3=0,int n4=0)
  {
    TinyVector<int,DMAX> shape(P.getTotalNum(),n2,n3,n4);
    return comp_samples.checkout_array<D>(P,name,shape);
  }


  //get trace functions
  //  used by estimators that require trace information
  inline TraceSample<TraceInt>* get_int_trace(const string& name)
  {
    return get_int_trace(default_domain,name);
  }
  inline TraceSample<TraceInt>* get_int_trace(const string& domain, const string& name)
  {
    return int_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceInt>* get_int_trace(const ParticleSet& P, const string& name)
  {
    return int_samples.get_trace(P.parentName(),name);
  }

  inline TraceSample<TraceReal>* get_real_trace(const string& name)
  {
    return get_real_trace(default_domain,name);
  }
  inline TraceSample<TraceReal>* get_real_trace(const string& domain, const string& name)
  {
    return real_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceReal>* get_real_trace(const ParticleSet& P, const string& name)
  {
    return real_samples.get_trace(P.parentName(),name);
  }

  inline TraceSample<TraceComp>* get_complex_trace(const string& name)
  {
    return get_complex_trace(default_domain,name);
  }
  inline TraceSample<TraceComp>* get_complex_trace(const string& domain, const string& name)
  {
    return comp_samples.get_trace(domain,name);
  }
  inline TraceSample<TraceComp>* get_complex_trace(const ParticleSet& P, const string& name)
  {
    return comp_samples.get_trace(P.parentName(),name);
  }


  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const string& name)
  {
    return get_int_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const string& domain, const string& name)
  {
    return int_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceInt>* get_int_combined_trace(const ParticleSet& P, const string& name)
  {
    return int_samples.get_combined_trace(P.parentName(),name);
  }

  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const string& name)
  {
    return get_real_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const string& domain, const string& name)
  {
    return real_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceReal>* get_real_combined_trace(const ParticleSet& P, const string& name)
  {
    return real_samples.get_combined_trace(P.parentName(),name);
  }

  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const string& name)
  {
    return get_complex_combined_trace(default_domain,name);
  }
  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const string& domain, const string& name)
  {
    return comp_samples.get_combined_trace(domain,name);
  }
  inline CombinedTraceSample<TraceComp>* get_complex_combined_trace(const ParticleSet& P, const string& name)
  {
    return comp_samples.get_combined_trace(P.parentName(),name);
  }




  //create combined trace out of existing traces
  inline void make_combined_trace(const string& name,vector<string>& names,vector<TraceReal>& weights)
  {
    if(traces_available && requests.find(name)!=requests.end())
    {
      if(verbose)
        app_log()<<"TraceManager::make_combined_trace "<<master_copy<<"  "<<name<<endl;
      real_buffer.make_combined_trace(name,names,weights);
    }
  }


  inline void check_clones(vector<TraceManager*>& clones)
  {
    if(writing_traces && clones.size()>0)
    {
      if(verbose)
        app_log()<<"TraceManager::check_clones"<<endl;
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
        APP_ABORT("TraceManager::check_clones  trace buffer widths of clones do not match\n  contiguous write is impossible\n  this was first caused by clones contributing particle traces from identical, but differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the TraceManager summaries printed above");
      }
    }
  }


  inline void reset_buffers()
  {
    if(writing_traces)
    {
      if(verbose)
        app_log()<<"TraceManager::reset_buffers "<<master_copy<<endl;
      int_buffer.reset();
      real_buffer.reset();
    }
  }


  //store the full sample from a single walker step in buffers
  inline void buffer_sample()
  {
    if(writing_traces)
    {
      if(verbose)
        app_log()<<" TraceManager::buffer_sample() "<<master_copy<<endl;
      int_buffer.collect_sample();
      real_buffer.collect_sample();
    }
  }


  //write buffered trace data to file
  inline void write_buffers(vector<TraceManager*>& clones)
  {
    if(master_copy)
    {
      if(writing_traces)
      {
        if(verbose)
          app_log()<<"TraceManager::write_buffers "<<master_copy<<endl;
        if(hdf_format)
        {
          write_buffers_hdf(clones);
        }
        if(adios_format)
        {
#ifdef HAVE_ADIOS
          write_buffers_adios(clones);
#else
          APP_ABORT("TraceManager::write_buffers (adios) ADIOS is not found");
#endif
          //app_log()<<"TraceManager::write_buffers (adios) has not yet been implemented"<<endl;
        }
      }
    }
    else
      APP_ABORT("TraceManager::write_buffers should not be called from non-master copy");
  }


  inline void open_file(vector<TraceManager*>& clones)
  {
    if(master_copy)
    {
      if(writing_traces)
      {
        if(verbose)
          app_log()<<"TraceManager::open_file "<<master_copy<<endl;
        //if(verbose)
        //  clones[0]->write_summary();
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
          app_log()<<"TraceManager::close_file "<<master_copy<<endl;
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



  inline void startRun(int blocks,vector<TraceManager*>& clones)
  {
    if(verbose)
      app_log()<<"TraceManager::startRun "<<master_copy<<endl;
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
      app_log()<<"TraceManager::stopRun "<<master_copy<<endl;
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
      app_log()<<"TraceManager::startBlock "<<master_copy<<endl;
    reset_buffers();
  }


  inline void stopBlock()
  {
    if(verbose)
      app_log()<<"TraceManager::stopBlock "<<master_copy<<endl;
  }


  inline void write_summary(string pad="  ")
  {
    string pad2 = pad+"  ";
    app_log()<<pad<<"TraceManager"<<endl;
    app_log()<<pad2<<"master_copy               = "<<master_copy              <<endl;
    app_log()<<pad2<<"traces_requested          = "<<traces_requested         <<endl;
    app_log()<<pad2<<"method_allows_traces      = "<<method_allows_traces     <<endl;
    app_log()<<pad2<<"traces_available          = "<<traces_available         <<endl;
    app_log()<<pad2<<"scalar_traces_requested   = "<<scalar_traces_requested  <<endl;
    app_log()<<pad2<<"particle_traces_requested = "<<particle_traces_requested<<endl;
    app_log()<<pad2<<"writing_traces            = "<<writing_traces           <<endl;
    app_log()<<pad2<<"format                    = "<<format                   <<endl;
    app_log()<<pad2<<"hdf format                = "<<hdf_format               <<endl;
    app_log()<<pad2<<"adios format              = "<<adios_format             <<endl;
    app_log()<<pad2<<"default_domain            = "<<default_domain           <<endl;
    app_log()<<pad2<<"use_scalar_defaults       = "<<scalar_defaults_set      <<endl;
    app_log()<<pad2<<"use_particle_defaults     = "<<particle_defaults_set    <<endl;
    set<string>::iterator it;
    app_log()<<pad2<<"scalar_requests           = ";
    if(!scalar_requests.empty())
    {
      for(it=scalar_requests.begin(); it!=scalar_requests.end(); ++it)
        app_log()<<*it<<" ";
    }
    app_log()<<endl;
    app_log()<<pad2<<"particle_requests         = ";
    if(!particle_requests.empty())
    {
      for(it=particle_requests.begin(); it!=particle_requests.end(); ++it)
        app_log()<<*it<<" ";
    }
    app_log()<<endl;
    app_log()<<pad2<<"requests                  = ";
    if(!requests.empty())
    {
      for(it=requests.begin(); it!=requests.end(); ++it)
        app_log()<<*it<<" ";
    }
    app_log()<<endl;
    int_buffer.write_summary(pad2);
    real_buffer.write_summary(pad2);
    app_log()<<pad<<"end TraceManager"<<endl;
  }


  //hdf file operations
  inline void open_hdf_file(vector<TraceManager*>& clones)
  {
    if(clones.size()==0)
      APP_ABORT("TraceManager::open_hdf_file  no trace clones exist, cannot open file");
    int nprocs = communicator->size();
    int rank = communicator->rank();
    char ptoken[32];
    string file_name = file_root;
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
      app_log()<<"TraceManager::open_hdf_file  opening traces hdf file "<<file_name<<endl;
    hdf_file = new hdf_archive(communicator,false);
    bool successful = hdf_file->create(file_name);
    if(!successful)
      APP_ABORT("TraceManager::open_hdf_file  failed to open hdf file "+file_name);
    // only clones have active buffers and associated data
    TraceManager& tm = *clones[0];
    tm.write_summary();
    tm.int_buffer.register_hdf_data(*hdf_file);
    tm.real_buffer.register_hdf_data(*hdf_file);
  }


  inline void write_buffers_hdf(vector<TraceManager*>& clones)
  {
    if(verbose)
      app_log()<<"TraceManager::write_buffers_hdf "<<master_copy<<endl;
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


#ifdef HAVE_ADIOS

  inline void initialize_adios(vector<TraceManager*>& clones)
  {
    //adios_init("qmc_adios.xml", communicator->getMPI());
    TraceManager& tm = *clones[0];
    tm.int_buffer.register_adios_data("w");
    tm.real_buffer.register_adios_data("a");
  }


  inline void print_adios(vector<TraceManager*>& clones)
  {
    app_log()<<"TraceManager::write_buffers_adios "<<master_copy<<endl;
    for(int ip=0; ip<clones.size(); ++ip)
    {
      TraceManager& tm = *clones[ip];
      //tm.int_buffer.write_summary();
      //tm.real_buffer.write_summary();
      app_log() << "Type: " << tm.real_buffer.type << endl;
      app_log() << "Buffer Size: " << tm.real_buffer.buffer.size() << " | " << tm.real_buffer.buffer.size(0) << "|" <<  tm.real_buffer.buffer.size(1) << endl;
    }
  }

  inline void write_buffers_adios(vector<TraceManager*>& clones)
  {
    //print_adios(clones);
    MPI_Comm comm = communicator->getMPI();
    int total_size = 0;
    for(int ip=0; ip<clones.size(); ++ip)
    {
      TraceManager& tm = *clones[ip];
      total_size += tm.real_buffer.buffer.size();
    }
    double adios_buffer[total_size];
    int curr_pos = 0;
    for(int ip=0; ip<clones.size(); ++ip)
    {
      TraceManager& tm = *clones[ip];
      for(vector<TraceReal>::iterator  iter= tm.real_buffer.buffer.begin(); iter != tm.real_buffer.buffer.end(); iter++)
      {
        adios_buffer[curr_pos++] = *iter;
      }
    }
    static bool write_flag = true;
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    if(write_flag)
    {
      adios_open(&adios_handle, "Traces", "traces.bp", "w", comm);
      //      write_flag = false;
    }
    else
    {
      adios_open(&adios_handle, "Traces", "traces.bp", "a", comm);
    }
    adios_groupsize = 4 + (total_size * 8);
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
    adios_write(adios_handle, "total_size", &total_size);
    adios_write(adios_handle, "buffer_contents", adios_buffer);
    adios_close(adios_handle);
#ifdef IO_PROFILE
    ADIOS_PROFILE::profile_adios_size(communicator, ADIOS_PROFILE::TRACES, adios_groupsize, adios_totalsize);
#endif
#ifdef ADIOS_VERIFY
    ADIOS_FILE *fp = adios_read_open_file("traces.bp",
                                          ADIOS_READ_METHOD_BP,
                                          OHMMS::Controller->getMPI());
    IO_VERIFY::adios_checkpoint_verify_variables(fp, "total_size", &total_size);
    IO_VERIFY::adios_trace_verify_local_variables(fp, "buffer_contents", adios_buffer);
    adios_read_close(fp);
#endif
  }


  inline void finalize_adios()
  {
    //adios_finalize(communicator->rank());
  }

#endif

};


}
#endif






