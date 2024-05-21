//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TRACEMANAGERNEW_H
#define QMCPLUSPLUS_TRACEMANAGERNEW_H


#if !defined(DISABLE_TRACEMANAGER)


#include <Configuration.h>
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "Utilities/IteratorUtility.h"
#include "ModernStringUtils.hpp"
#include "Message/Communicate.h"
#include "hdf/hdf_archive.h"
#include "Concurrency/OpenMP.h"

#include <algorithm>
#include <array>
#include <map>
#include <memory>
#include <set>

namespace qmcplusplus
{
//#define TRACE_CHECK

const unsigned int DMAXNEW = 4;
using TraceIntNew          = long;
using TraceRealNew         = OHMMS_PRECISION;
using TraceCompNew         = std::complex<TraceRealNew>;


struct TraceQuantityNew
{
  std::string name;
  bool default_quantity;
  bool scalar_available;
  bool array_available;
  bool scalar_stream_requested;
  bool array_stream_requested;
  bool scalar_write_requested;
  bool array_write_requested;
  bool stream_scalar;
  bool stream_array;
  bool write_scalar;
  bool write_array;

  inline TraceQuantityNew()
  {
    default_quantity        = false;
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

  inline void incorporate(const TraceQuantityNew& other)
  {
    if (name != other.name)
    {
      APP_ABORT("TraceQuantityNew::incorporate\n  cannot merge quantities with differing names\n  names: " + name + " " +
                other.name);
    }
    default_quantity |= other.default_quantity;
    scalar_available |= other.scalar_available;
    array_available |= other.array_available;
    scalar_stream_requested |= other.scalar_stream_requested;
    array_stream_requested |= other.array_stream_requested;
  }
};


//means of control of traces
//  which quantities are available, streaming, writing, etc
//  handles global (user + all observable) trace request
//medium of exchange of trace information between trace manager and observables
//  handles individual trace requests from observables
struct TraceRequestNew
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
  bool streaming_default_arrays;
  bool writing_default_scalars;
  bool writing_default_arrays;


  //quantities made available or requested for streaming or writing
  std::map<std::string, TraceQuantityNew> quantities;

  //used to screen checked out quantities for writing
  std::string scalar_domain;

  inline TraceRequestNew() { reset(); }

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
  }

  inline void set_scalar_domain(const std::string& domain) { scalar_domain = domain; }

  inline bool screen_sample(const std::string& domain, const std::string& name, bool& write)
  {
    bool scalar  = domain == scalar_domain;
    bool present = quantities.find(name) != quantities.end();
    bool stream  = false;
    write        = false;
    if (present)
    {
      TraceQuantityNew& q = quantities[name];
      if (scalar)
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
    app_log()<<"JTK: TraceRequestNew::screen_sample "<<name<<" "<<present<<" "<<scalar<<" "<<stream<<" "<<write<<std::endl;
    return stream;
  }

  //Contributor API (OperatorBase and others)
  //declare that scalars are available for a quantity
  inline void contribute_scalar(const std::string& name, bool default_quantity = false)
  {
    app_log()<<"JTK: TraceRequestNew::contribute_scalar (op) "<<name<<std::endl;
    guarantee_presence(name);
    quantities[name].scalar_available = true;
    if (default_quantity)
      quantities[name].default_quantity = true;
  }

  //declare that arrays are available for a quantity
  inline void contribute_array(const std::string& name, bool default_quantity = false)
  {
    app_log()<<"JTK: TraceRequestNew::contribute_array (op)"<<name<<std::endl;
    guarantee_presence(name);
    quantities[name].array_available = true;
    if (default_quantity)
      quantities[name].default_quantity = true;
  }

  //declare that a scalar quantity is desired for streaming
  inline void request_scalar(const std::string& name, bool write = false)
  {
    app_log()<<"JTK: TraceRequestNew::request_scalar (op)"<<name<<std::endl;
    guarantee_presence(name);
    quantities[name].scalar_stream_requested = true;
    if (write)
      quantities[name].scalar_write_requested = true;
  }

  //declare that a array quantity is desired for streaming
  inline void request_array(const std::string& name, bool write = false)
  {
    app_log()<<"JTK: TraceRequestNew::request_array (op)"<<name<<std::endl;
    guarantee_presence(name);
    quantities[name].array_stream_requested = true;
    if (write)
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
    TraceQuantityNew& q = quantities[name];
    return q.stream_scalar || q.stream_array;
  }


  //TraceManager API
  //declare that scalar quantities are desired for streaming
  inline void request_scalar(const std::set<std::string>& names, bool write = false)
  {
    app_log()<<"JTK: TraceRequestNew::request_scalar (tm)"<<std::endl;
    std::set<std::string>::iterator name;
    for (name = names.begin(); name != names.end(); ++name)
      request_scalar(*name, write);
  }

  //declare that array quantities are desired for streaming
  inline void request_array(const std::set<std::string>& names, bool write = false)
  {
    app_log()<<"JTK: TraceRequestNew::request_array (tm)"<<std::endl;
    std::set<std::string>::iterator name;
    for (name = names.begin(); name != names.end(); ++name)
      request_array(*name, write);
  }

  //merge in all quantities from a contributor request
  inline void incorporate(TraceRequestNew& other)
  {
    std::map<std::string, TraceQuantityNew>::iterator it;
    for (it = other.quantities.begin(); it != other.quantities.end(); ++it)
    {
      const TraceQuantityNew& q = it->second;
      if (quantity_present(q.name))
        quantities[q.name].incorporate(q);
      else
        quantities[q.name] = q;
    }
  }

  //balance requests with availability and turn streaming/writing on/off
  inline void determine_stream_write()
  {
    std::map<std::string, TraceQuantityNew>::iterator it;
    for (it = quantities.begin(); it != quantities.end(); ++it)
    {
      TraceQuantityNew& q = it->second;

      q.stream_scalar =
          allow_streams && scalars_on && q.scalar_available && (q.scalar_stream_requested || stream_all_scalars);

      q.write_scalar = q.stream_scalar && (q.scalar_write_requested || write_all_scalars);

      q.stream_array =
          allow_streams && arrays_on && q.array_available && (q.array_stream_requested || stream_all_arrays);

      q.write_array = q.stream_array && (q.array_write_requested || write_all_arrays);
    }
    // default quantities stream and write if any others do
    streaming_default_scalars = false;
    writing_default_scalars   = false;
    streaming_default_arrays  = false;
    writing_default_arrays    = false;
    for (it = quantities.begin(); it != quantities.end(); ++it)
    {
      TraceQuantityNew& q = it->second;
      streaming_default_scalars |= q.stream_scalar;
      writing_default_scalars |= q.write_scalar;
      streaming_default_arrays |= q.stream_array;
      writing_default_arrays |= q.write_array;
    }
    for (it = quantities.begin(); it != quantities.end(); ++it)
    {
      TraceQuantityNew& q = it->second;
      if (q.default_quantity)
      {
        q.stream_scalar = streaming_default_scalars;
        q.write_scalar  = writing_default_scalars;
        q.stream_array  = streaming_default_arrays;
        q.write_array   = writing_default_arrays;
      }
    }
  }

  //relay updated streaming information to contributor
  inline void relay_stream_info(TraceRequestNew& other)
  {
    other.allow_streams             = allow_streams;
    other.allow_writes              = allow_writes;
    other.scalars_on                = scalars_on;
    other.arrays_on                 = arrays_on;
    other.stream_all_scalars        = stream_all_scalars;
    other.stream_all_arrays         = stream_all_arrays;
    other.write_all_scalars         = write_all_scalars;
    other.write_all_arrays          = write_all_arrays;
    other.streaming_default_scalars = streaming_default_scalars;
    other.streaming_default_arrays  = streaming_default_arrays;
    other.writing_default_scalars   = writing_default_scalars;
    other.writing_default_arrays    = writing_default_arrays;
    std::map<std::string, TraceQuantityNew>::iterator it;
    for (it = other.quantities.begin(); it != other.quantities.end(); ++it)
    {
      TraceQuantityNew& q = it->second;
      check_presence(q.name);
      q = quantities[q.name];
    }
  }

  inline void report()
  {
    //app_log()<<"\n  TraceRequestNew"<< std::endl;
    app_log() << "    allow_streams             = " << allow_streams << std::endl;
    app_log() << "    allow_writes              = " << allow_writes << std::endl;
    app_log() << "    scalars_on                = " << scalars_on << std::endl;
    app_log() << "    arrays_on                 = " << arrays_on << std::endl;
    app_log() << "    stream_all_scalars        = " << stream_all_scalars << std::endl;
    app_log() << "    stream_all_arrays         = " << stream_all_arrays << std::endl;
    app_log() << "    write_all_scalars         = " << write_all_scalars << std::endl;
    app_log() << "    write_all_arrays          = " << write_all_arrays << std::endl;
    app_log() << "    streaming_default_scalars = " << streaming_default_scalars << std::endl;
    app_log() << "    streaming_default_arrays  = " << streaming_default_arrays << std::endl;
    app_log() << "    writing_default_scalars   = " << writing_default_scalars << std::endl;
    app_log() << "    writing_default_arrays    = " << writing_default_arrays << std::endl;

    write_selected("scalars available", "scalar_available");
    write_selected("arrays available", "array_available");

    write_selected("scalar streams requested", "scalar_stream_requested");
    write_selected("array streams requested", "array_stream_requested");

    write_selected("scalar writes requested", "scalar_write_requested");
    write_selected("array writes requested", "array_write_requested");

    write_selected("scalar streams occurring", "stream_scalar");
    write_selected("array streams occurring", "stream_array");

    write_selected("scalar writes occurring", "write_scalar");
    write_selected("array writes occurring", "write_array");
  }

  inline void write_selected(const std::string& header, const std::string& selector)
  {
    app_log() << "    " << header << ":";
    int n = 0;
    std::map<std::string, TraceQuantityNew>::iterator it;
    for (it = quantities.begin(); it != quantities.end(); ++it)
    {
      TraceQuantityNew& q = it->second;
      bool selected    = false;
      if (selector == "scalar_available")
        selected = q.scalar_available;
      else if (selector == "array_available")
        selected = q.array_available;
      else if (selector == "scalar_stream_requested")
        selected = q.scalar_stream_requested;
      else if (selector == "array_stream_requested")
        selected = q.array_stream_requested;
      else if (selector == "scalar_write_requested")
        selected = q.scalar_write_requested;
      else if (selector == "array_write_requested")
        selected = q.array_write_requested;
      else if (selector == "stream_scalar")
        selected = q.stream_scalar;
      else if (selector == "stream_array")
        selected = q.stream_array;
      else if (selector == "write_scalar")
        selected = q.write_scalar;
      else if (selector == "write_array")
        selected = q.write_array;
      else
        APP_ABORT("TraceRequestNew::write_selected  unrecognized selector: " + selector);
      if (selected)
      {
        if (n % 5 == 0)
          app_log() << std::endl << "      ";
        n++;
        app_log() << " " << q.name;
      }
    }
    app_log() << std::endl;
  }


  //private (internal) API
  //query whether a quantity is present
  inline bool quantity_present(const std::string& name) { return quantities.find(name) != quantities.end(); }

  //create a quantity if it is not present
  inline void guarantee_presence(const std::string& name, bool combined = false)
  {
    if (!quantity_present(name))
    {
      TraceQuantityNew q;
      q.name           = name;
      quantities[name] = q;
    }
  }

  //abort if a quantity is not present
  inline void check_presence(const std::string& name)
  {
    if (!quantity_present(name))
    {
      APP_ABORT("TraceRequestNew::check_presence  quantity " + name + " is not present");
    }
  }

  //query whether any quantities are streaming
  inline bool streaming() { return streaming_default_scalars || streaming_default_arrays; }

  //query whether any quantities are writing
  inline bool writing() { return writing_default_scalars || writing_default_arrays; }

  //query whether any scalar quantities are streaming
  inline bool streaming_scalars() { return streaming_default_scalars; }

  //query whether any array quantities are streaming
  inline bool streaming_arrays() { return streaming_default_arrays; }
};


template<typename T>
struct TraceSampleNew
{
  std::string domain;
  std::string name;
  int index;
  bool array_trace;
  int dimension;
  int size;
  int unit_size;
  int data_size;
  TinyVector<int, DMAXNEW> shape;
  Vector<T>& sample;
  bool write;
  int buffer_start, buffer_end;
  std::map<std::string, TraceIntNew> meta_int;
  std::map<std::string, TraceRealNew> meta_real;
  std::map<std::string, std::string> meta_string;
  bool verbose;

  inline TraceSampleNew(const std::string& sdomain, const std::string& sname, int sindex, int sdim, Vector<T>& ssample)
      : sample(ssample), verbose(false)
  {
    initialize(sdomain, sname, sindex, sdim);
  }


  inline TraceSampleNew(const std::string& sdomain,
                     const std::string& sname,
                     int sindex,
                     int sdim,
                     TinyVector<int, DMAXNEW> sshape,
                     Vector<T>& ssample)
      : sample(ssample), verbose(false)
  {
    initialize(sdomain, sname, sindex, sdim);
    shape = sshape;
    size  = sample.size();
    check_shape();
  }

  inline virtual ~TraceSampleNew() = default;

  inline void initialize(const std::string& sdomain, const std::string& sname, int sindex, int sdim)
  {
    domain = sdomain, name = sname;
    dimension    = sdim;
    index        = sindex;
    array_trace  = false;
    write        = false;
    buffer_start = -1;
    buffer_end   = -1;
  }


  inline void set_unit_size(int usize) { unit_size = usize; }


  inline void set_data_size() { data_size = size * unit_size; }


  inline void check_shape()
  {
    bool correct_shape     = true;
    bool correct_dimension = dimension <= DMAXNEW;
    if (correct_dimension)
    {
      int tsize = 1;
      for (int d = 0; d < dimension; ++d)
      {
        tsize *= shape[d];
        correct_dimension = correct_dimension && shape[d] > 0;
      }
      correct_shape = tsize == size;
    }
    if (!correct_dimension)
      APP_ABORT("TraceSampleNew::check_shape dimension of sample array is incorrect");
    if (!correct_shape)
      APP_ABORT("TraceSampleNew::check_shape shape and size of sample array do not match");
  }


  inline bool same_shape(TraceSampleNew<T>* other)
  {
    bool same = dimension == other->dimension && size == other->size;
    if (same)
      for (int d = 0; d < dimension; ++d)
        same = same && shape[d] == other->shape[d];
    return same;
  }

  inline void set_buffer_range(int& bstart)
  {
    set_data_size();
    if (write)
    {
      buffer_start = bstart;
      buffer_end   = bstart + data_size;
      bstart       = buffer_end;
    }
  }


  inline T sum()
  {
    T s(0);
    for (int i = 0; i < sample.size(); ++i)
      s += sample[i];
    return s;
  }

  inline void write_summary(int ind = -1, std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    if (ind == -1)
      app_log() << pad << " TraceSampleNew " << name << std::endl;
    else
      app_log() << pad << ind << " TraceSampleNew " << name << std::endl;
    app_log() << pad2 << "domain       = " << domain << std::endl;
    app_log() << pad2 << "name         = " << name << std::endl;
    app_log() << pad2 << "index        = " << index << std::endl;
    app_log() << pad2 << "array_trace  = " << array_trace << std::endl;
    app_log() << pad2 << "dimension    = " << dimension << std::endl;
    app_log() << pad2 << "size         = " << size << std::endl;
    app_log() << pad2 << "unit_size    = " << unit_size << std::endl;
    app_log() << pad2 << "data_size    = " << data_size << std::endl;
    app_log() << pad2 << "shape        = " << shape << std::endl;
    app_log() << pad2 << "write        = " << write << std::endl;
    app_log() << pad2 << "buffer range = [" << buffer_start << "," << buffer_end << ")" << std::endl;
  }
};



template<typename T>
bool TraceSampleNew_comp(TraceSampleNew<T>* left, TraceSampleNew<T>* right)
{
  return left->data_size < right->data_size;
}


template<typename T>
struct TraceSampleNews
{
  std::vector<TraceSampleNew<T>*> samples;
  std::map<std::string, std::map<std::string, int>> sample_indices;
  std::vector<TraceSampleNew<T>*> ordered_samples;
  bool verbose;

  inline TraceSampleNews() : verbose(false) {}

  inline ~TraceSampleNews() { finalize(); }


  inline void set_verbose(bool v) { verbose = v; }


  inline int size() { return samples.size(); }


  inline void assign_sample_index(const std::string& domain, const std::string& name, int index, std::string label = "")
  {
    if (sample_indices.count(domain) > 0 && sample_indices[domain].count(name) > 0)
    {
      APP_ABORT("TraceSampleNews::checkout " + label + " variable " + name + " already exists in domain " + domain);
    }
    else
    {
      sample_indices[domain][name] = index;
    }
  }

  template<int D>
  inline Array<T, D>* checkout_array(const std::string& domain, const std::string& name, TinyVector<int, DMAXNEW> shape)
  {
    int index = samples.size();
    assign_sample_index(domain, name, index, "array");
    std::array<size_t, D> subshape;
    for (int idim = 0; idim < D; idim++)
      subshape[idim] = shape[idim];
    Array<T, D>* a    = new Array<T, D>(subshape);
    TraceSampleNew<T>* s = new TraceSampleNew<T>(domain, name, index, D, shape, a->storage());
    samples.push_back(s);
    if (verbose)
      app_log() << "TraceSampleNews::checkout_array  " << domain << " " << name << " " << index << std::endl;
    return a;
  }


  template<int D>
  inline Array<T, D>* checkout_array(const ParticleSet& P, const std::string& name, TinyVector<int, DMAXNEW> shape)
  {
    const std::string& domain = P.parentName();
    int index                 = samples.size();
    assign_sample_index(domain, name, index, "array");
    std::array<size_t, D> subshape;
    for (int idim = 0; idim < D; idim++)
      subshape[idim] = shape[idim];
    Array<T, D>* a    = new Array<T, D>(subshape);
    TraceSampleNew<T>* s = new TraceSampleNew<T>(domain, name, index, D, shape, a->storage());
    samples.push_back(s);
    s->array_trace = true;
    if (verbose)
      app_log() << "TraceSampleNews::checkout_array  " << domain << " " << name << " " << index << std::endl;
    return a;
  }


  inline TraceSampleNew<T>* get_trace(const std::string& domain, const std::string& name)
  {
    TraceSampleNew<T>* ts = NULL;
    for (int i = 0; i < samples.size(); ++i)
    {
      TraceSampleNew<T>& tsc = *samples[i];
      if (tsc.domain == domain && tsc.name == name)
      {
        ts = &tsc;
        break;
      }
    }
    if (ts == NULL)
      APP_ABORT("TraceSampleNews::get_trace  failed to get trace for quantity " + name + " in domain " + domain);
    return ts;
  }



  inline void set_unit_size(int usize)
  {
    for (int i = 0; i < samples.size(); i++)
      samples[i]->set_unit_size(usize);
  }


  inline void screen_writes(TraceRequestNew& request)
  {
    for (int i = 0; i < samples.size(); i++)
    {
      TraceSampleNew<T>& s = *samples[i];
      bool stream       = request.screen_sample(s.domain, s.name, s.write);
      if (verbose)
        app_log() << "TraceRequestNew screening " << s.name << " in domain " << s.domain << ". stream: " << stream
                  << " write: " << s.write << std::endl;
      if (!stream)
        app_log() << "warning: quantity " + s.name + " in domain " + s.domain +
                " was not requested but is streaming anyway"
                  << std::endl;
    }
  }


  inline void order_by_size()
  {
    for (int i = 0; i < samples.size(); i++)
      samples[i]->set_data_size();
    ordered_samples.resize(samples.size());
    copy(samples.begin(), samples.end(), ordered_samples.begin());
    sort(ordered_samples.begin(), ordered_samples.end(), TraceSampleNew_comp<T>);
  }


  inline void set_buffer_ranges(int& starting_index)
  {
    for (int i = 0; i < ordered_samples.size(); i++)
    {
      TraceSampleNew<T>& sample = *ordered_samples[i];
      sample.set_buffer_range(starting_index);
    }
  }


  inline int total_size()
  {
    int s = 0;
    for (int i = 0; i < samples.size(); i++)
      s += samples[i]->sample.size() * samples[i]->unit_size;
    return s;
  }


  inline int min_buffer_index()
  {
    int min_index = 2000000000;
    for (int i = 0; i < samples.size(); i++)
      min_index = std::min(min_index, samples[i]->buffer_start);
    return min_index;
  }


  inline int max_buffer_index()
  {
    int max_index = -1;
    for (int i = 0; i < samples.size(); i++)
      max_index = std::max(max_index, samples[i]->buffer_end);
    return max_index;
  }


  inline void finalize()
  {
    delete_iter(samples.begin(), samples.end());
    samples.resize(0);
    ordered_samples.resize(0);
    sample_indices.clear();
  }


  inline void register_hdf_data(hdf_archive& f)
  {
    std::map<std::string, std::map<std::string, int>>::iterator it;
    std::map<std::string, int>::iterator it2;
    for (it = sample_indices.begin(); it != sample_indices.end(); it++)
    {
      const std::string& domain           = it->first;
      std::map<std::string, int>& indices = it->second;
      f.push(domain);
      for (it2 = indices.begin(); it2 != indices.end(); ++it2)
      {
        const std::string& quantity = it2->first;
        TraceSampleNew<T>& sample      = *samples[it2->second];
        if (sample.write)
        {
          f.push(quantity);
          f.write(sample.dimension, "dimension");
          f.write(sample.shape, "shape");
          f.write(sample.size, "size");
          f.write(sample.unit_size, "unit_size");
          f.write(sample.buffer_start, "row_start");
          f.write(sample.buffer_end, "row_end");
          f.pop();
        }
      }
      f.pop();
    }
  }


  inline void write_summary(std::string type, std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    std::string pad3 = pad2 + "  ";
    std::string pad4 = pad3 + "  ";
    app_log() << pad << "TraceSampleNews<" << type << ">" << std::endl;
    app_log() << pad2 << "nsamples          = " << samples.size() << std::endl;
    app_log() << pad2 << "sample_indices" << std::endl;
    std::map<std::string, std::map<std::string, int>>::iterator it;
    std::map<std::string, int>::iterator it2;
    for (it = sample_indices.begin(); it != sample_indices.end(); it++)
    {
      const std::string& domain           = it->first;
      std::map<std::string, int>& indices = it->second;
      app_log() << pad3 << domain << std::endl;
      for (it2 = indices.begin(); it2 != indices.end(); ++it2)
        app_log() << pad4 << it2->first << " = " << it2->second << std::endl;
    }
    app_log() << pad2 << "end sample_indices" << std::endl;
    app_log() << pad2 << "samples" << std::endl;
    for (int i = 0; i < ordered_samples.size(); ++i)
      ordered_samples[i]->write_summary(i, pad3);
    app_log() << pad2 << "end samples" << std::endl;
    app_log() << pad << "end TraceSampleNews<" << type << ">" << std::endl;
  }


  inline void user_report(const std::string& type, const std::string& pad = "  ")
  {
    std::string pad2 = pad + "  ";
    app_log() << pad << type << " traces provided by estimators" << std::endl;
    std::map<std::string, std::map<std::string, int>>::iterator it;
    std::map<std::string, int>::iterator it2;
    for (it = sample_indices.begin(); it != sample_indices.end(); it++)
    {
      const std::string& domain           = it->first;
      std::map<std::string, int>& indices = it->second;
      app_log() << pad2 << "domain " << domain << ":  ";
      int n = 0;
      for (it2 = indices.begin(); it2 != indices.end(); ++it2)
      {
        if (n % 5 == 0)
          app_log() << std::endl << pad2 << "  ";
        n++;
        const std::string& quantity = it2->first;
        app_log() << quantity << " ";
      }
      app_log() << std::endl;
    }
  }
};


template<typename T>
struct TraceBufferNew
{
  bool has_complex;
  TraceSampleNews<T>* samples;
  TraceSampleNews<std::complex<T>>* complex_samples;
  std::string type;
  Array<T, 2> buffer;
  bool verbose;

  //hdf variables
  std::string top;
  hsize_t dims[2];
  hsize_t hdf_file_pointer;


  TraceBufferNew() : samples(0), complex_samples(0), verbose(false)
  {
    type        = "?";
    has_complex = false;
    reset();
  }


  inline void set_verbose(bool v) { verbose = v; }


  inline void set_type(std::string stype)
  {
    type = stype;
    top  = type + "_data";
  }


  inline void reset() { buffer.resize(0, buffer.size(1)); }


  inline void set_samples(TraceSampleNews<T>& s) { samples = &s; }


  inline void set_samples(TraceSampleNews<std::complex<T>>& s)
  {
    complex_samples = &s;
    has_complex     = true;
  }


  inline void order_and_resize()
  {
    //put the sample data in size order
    samples->set_unit_size(1);
    samples->order_by_size();
    if (has_complex)
    {
      complex_samples->set_unit_size(2);
      complex_samples->order_by_size();
    }
    //assign buffer ranges to each sample
    int sample_size = 0;
    samples->set_buffer_ranges(sample_size);
    if (has_complex)
      complex_samples->set_buffer_ranges(sample_size);
#if defined(TRACE_CHECK)
    test_buffer_write(sample_size);
#endif
    //resize the buffer
    int nsamples_init = 1;
    buffer.resize(nsamples_init, sample_size);
  }


  inline bool same_as(TraceBufferNew<T>& ref) { return buffer.size(1) == ref.buffer.size(1); }


  inline void collect_sample()
  {
    if (verbose)
      app_log() << " TraceBufferNew<" << type << ">::collect_sample()" << std::endl;
    //make more room, if necessary
    int nrows    = buffer.size(0);
    int row_size = buffer.size(1);
    if (row_size > 0)
    {
      //make space for the row, if necessary
      int current_row = nrows;
      nrows++;
      // resizing buffer(type Array) doesn't preserve data. Thus keep old data and copy over
      auto buffer_old(buffer);
      buffer.resize(nrows, row_size);
      std::copy_n(buffer_old.data(), buffer_old.size(), buffer.data());
      if (verbose)
        app_log() << "  increasing # of rows to " << nrows << std::endl;
      //collect data from all samples into the buffer row
      {
        std::vector<TraceSampleNew<T>*>& ordered_samples = samples->ordered_samples;
        for (int s = 0; s < ordered_samples.size(); s++)
        {
          TraceSampleNew<T>& tsample = *ordered_samples[s];
          if (tsample.write)
          {
            auto& sample = tsample.sample;
            for (int i = 0; i < sample.size(); ++i)
              buffer(current_row, tsample.buffer_start + i) = sample[i];
          }
        }
      }
      if (has_complex)
      {
        std::vector<TraceSampleNew<std::complex<T>>*>& ordered_samples = complex_samples->ordered_samples;
        for (int s = 0; s < ordered_samples.size(); s++)
        {
          TraceSampleNew<std::complex<T>>& tsample = *ordered_samples[s];
          if (tsample.write)
          {
            auto& sample = tsample.sample;
            for (int i = 0, ib = 0; i < sample.size(); ++i, ib += 2)
            {
              buffer(current_row, tsample.buffer_start + ib)     = sample[i].real();
              buffer(current_row, tsample.buffer_start + ib + 1) = sample[i].imag();
            }
          }
        }
      }
#if defined(TRACE_CHECK)
      test_buffer_collect(current_row);
#endif
    }
  }


  inline void write() { APP_ABORT("TraceBufferNew::write has not yet been implemented"); }


  inline void write_summary(std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    app_log() << pad << "TraceBufferNew<" << type << ">" << std::endl;
    app_log() << pad2 << "nrows       = " << buffer.size(0) << std::endl;
    app_log() << pad2 << "row_size    = " << buffer.size(1) << std::endl;
    app_log() << pad2 << "has_complex = " << has_complex << std::endl;
    samples->write_summary(type, pad2);
    if (has_complex)
      complex_samples->write_summary("complex " + type, pad2);
    app_log() << pad << "end TraceBufferNew<" << type << ">" << std::endl;
  }


  inline void user_report(const std::string& pad = "  ")
  {
    samples->user_report(type, pad);
    if (has_complex)
      complex_samples->user_report("complex " + type, pad);
  }

  inline void register_hdf_data(hdf_archive& f)
  {
    f.push(top);
    f.push("layout");
    samples->register_hdf_data(f);
    if (has_complex)
      complex_samples->register_hdf_data(f);
    f.pop();
    f.pop();
    if (!f.open_groups())
      APP_ABORT("TraceBufferNew<" + type +
                ">::register_hdf_data() some hdf groups are still open at the end of registration");
    hdf_file_pointer = 0;
  }


  inline void write_hdf(hdf_archive& f) { write_hdf(f, hdf_file_pointer); }


  inline void write_hdf(hdf_archive& f, hsize_t& file_pointer)
  {
    if (verbose)
      app_log() << "TraceBufferNew<" << type << ">::write_hdf() " << file_pointer << " " << buffer.size(0) << " "
                << buffer.size(1) << std::endl;
    dims[0] = buffer.size(0);
    dims[1] = buffer.size(1);
    if (dims[0] > 0)
    {
      f.push(top);
      h5d_append(f.top(), "traces", file_pointer, buffer.dim(), dims, buffer.data());
      f.pop();
    }
    f.flush();
  }


  inline void test_buffer_write(int sample_size)
  {
    //check that the size is correct
    int ssize = samples->total_size();
    if (has_complex)
      ssize += complex_samples->total_size();
    if (sample_size != ssize)
    {
      app_log() << "sample_size = " << sample_size << "\ntotal_size = " << ssize << std::endl;
      APP_ABORT("TraceBufferNew::test_buffer_write sample_size and total_size do not match");
    }
    //check that buffer indices fall in the expected range
    int nsamples  = samples->size();
    int min_index = samples->min_buffer_index();
    int max_index = samples->max_buffer_index();
    if (has_complex)
    {
      nsamples += complex_samples->size();
      min_index = std::min(min_index, complex_samples->min_buffer_index());
      max_index = std::max(max_index, complex_samples->max_buffer_index());
    }
    if (nsamples > 0)
    {
      if (min_index != 0)
        APP_ABORT("TraceBufferNew::test_buffer_write min_index!=0\n  min_index=" << min_index);
      if (max_index != sample_size)
        APP_ABORT("TraceBufferNew::test_buffer_write max_index!=sample_size");
      //check that no overlap exists in writes to buffer
      Array<int, 2> test_buffer;
      test_buffer.resize(1, sample_size);
      std::fill(test_buffer.begin(), test_buffer.end(), 0);
      int row      = 0;
      int row_size = test_buffer.size(1);
      int offset   = row * row_size;
      int* loc1    = &test_buffer(offset);
      int* loc2    = &test_buffer(row, 0);
      if (loc1 != loc2)
        APP_ABORT("TraceBufferNew::test_buffer_write serialized buffer offset improperly computed");
      {
        int boffset;
        std::vector<TraceSampleNew<T>*>& ordered_samples = samples->ordered_samples;
        for (int s = 0; s < ordered_samples.size(); s++)
        {
          TraceSampleNew<T>& tsample = *ordered_samples[s];
          std::vector<T>& sample  = tsample.sample;
          boffset                 = offset + tsample.buffer_start;
          for (int i = 0; i < sample.size(); ++i)
            test_buffer(boffset + i) = 1;
        }
      }
      if (has_complex)
      {
        int boffset;
        std::vector<TraceSampleNew<std::complex<T>>*>& ordered_samples = complex_samples->ordered_samples;
        for (int s = 0; s < ordered_samples.size(); s++)
        {
          TraceSampleNew<std::complex<T>>& tsample = *ordered_samples[s];
          std::vector<std::complex<T>>& sample  = tsample.sample;
          boffset                               = offset + tsample.buffer_start;
          for (int i = 0, ib = 0; i < sample.size(); ++i, ib += tsample.unit_size)
          {
            test_buffer(boffset + ib)     = 1;
            test_buffer(boffset + ib + 1) = 1;
          }
        }
      }
      //app_log()<<"test_buffer:"<< std::endl;
      //for(int i=0;i<row_size;++i)
      //  app_log()<<"  "<<i<<"  "<<test_buffer(row,i)<< std::endl;
      for (int i = 0; i < row_size; ++i)
        if (!test_buffer(row, i))
          APP_ABORT("TraceBufferNew::test_buffer_write write to row is not contiguous");
    }
  }


  inline void test_buffer_collect(int current_row)
  {
    if (verbose)
      app_log() << "TraceBufferNew::test_buffer_collect" << std::endl;
    std::string scalars = "scalars";
    std::map<std::string, std::map<std::string, int>>::iterator dom;
    std::map<std::string, int>::iterator var;
    std::map<std::string, std::map<std::string, int>>& domains = samples->sample_indices;
    std::vector<TraceSampleNew<T>*>& tsamples                     = samples->samples;
    std::map<std::string, int>& scalar_vars                    = domains[scalars];
    for (var = scalar_vars.begin(); var != scalar_vars.end(); var++)
    {
      const std::string& name = var->first;
      TraceSampleNew<T>& sample  = *tsamples[var->second];
      T value                 = buffer(current_row, sample.buffer_start);
      T svalue                = 0;
      bool any_present        = false;
      for (dom = domains.begin(); dom != domains.end(); dom++)
      {
        const std::string& domain = dom->first;
        if (domain != scalars)
        {
          std::map<std::string, int>& vars = dom->second;
          if (vars.count(name) > 0)
          {
            any_present             = true;
            TraceSampleNew<T>& ssample = *tsamples[vars[name]];
            int start               = ssample.buffer_start;
            int end                 = ssample.buffer_end;
            for (int i = start; i < end; i++)
              svalue += buffer(current_row, i);
          }
        }
      }
      if (any_present)
      {
        if (verbose)
          app_log() << "  " << name << " " << value << " " << svalue << std::endl;
      }
    }
  }
};




struct TraceManagerState 
{
  bool method_allows_traces;
  TraceRequestNew request;
  bool streaming_traces;
  bool writing_traces;
  int throttle;
  bool verbose;
  std::string default_domain;
};



class TraceCollector
{
public:
  //collections of samples for a single walker step
  //  the associated arrays will be updated following evaluate
  TraceSampleNews<TraceIntNew> int_samples;
  TraceSampleNews<TraceRealNew> real_samples;
  TraceSampleNews<TraceCompNew> comp_samples;

  //buffers for storing samples
  // single row of buffer is a single sample from one walker
  // number of rows adjusts to accommodate walker samples
  TraceBufferNew<TraceIntNew> int_buffer;
  TraceBufferNew<TraceRealNew> real_buffer;

  TraceRequestNew request;

  std::string default_domain;
  bool method_allows_traces;
  bool streaming_traces;
  bool writing_traces;
  int throttle;
  bool verbose;

  TraceCollector() : verbose(false)
  {
    reset_permissions();
    throttle       = 1;
    default_domain = "scalars";
    request.set_scalar_domain(default_domain);
    int_buffer.set_type("int");
    real_buffer.set_type("real");
    int_buffer.set_samples(int_samples);
    real_buffer.set_samples(real_samples);
    real_buffer.set_samples(comp_samples);
  }


  inline void transfer_state_from(const TraceManagerState& tms)
  {
    method_allows_traces = tms.method_allows_traces;
    request              = tms.request;
    streaming_traces     = tms.streaming_traces;
    writing_traces       = tms.writing_traces;
    throttle             = tms.throttle;
    verbose              = tms.verbose;
    default_domain       = tms.default_domain;
  }


  inline void distribute()
  {
    app_log()<<"JTK: TraceCollector::distribute"<<std::endl;
    int_samples.set_verbose(verbose);
    real_samples.set_verbose(verbose);
    comp_samples.set_verbose(verbose);
    int_buffer.set_verbose(verbose);
    real_buffer.set_verbose(verbose);
  }


  inline void reset_permissions()
  {
    app_log()<<"JTK: TraceCollector::reset_permissions"<<std::endl;
    method_allows_traces = false;
    streaming_traces     = false;
    writing_traces       = false;
    verbose              = false;
    request.reset();
  }


  inline void update_status()
  {
    app_log()<<"JTK: TraceCollector::update_status"<<std::endl;
    streaming_traces = request.streaming();
    writing_traces   = request.writing();
  }


  inline void screen_writes()
  {
    app_log()<<"JTK: TraceCollector::screen_writes"<<std::endl;
    int_samples.screen_writes(request);
    real_samples.screen_writes(request);
    comp_samples.screen_writes(request);
  }

  inline void initialize_traces()
  {
    app_log()<<"JTK: TraceCollector::initialize_traces "<<streaming_traces<<std::endl;
    if (streaming_traces)
    {
      if (verbose)
        app_log() << "TraceCollector::initialize_traces " << std::endl;
      //organize trace samples and initialize buffers
      if (writing_traces)
      {
        int_buffer.order_and_resize();
        real_buffer.order_and_resize();
      }
    }
  }


  inline void finalize_traces()
  {
    app_log()<<"JTK: TraceCollector::finalize_traces "<<std::endl;
    if (verbose)
      app_log() << "TraceCollector::finalize_traces " << std::endl;
    int_samples.finalize();
    real_samples.finalize();
    comp_samples.finalize();
  }


  //checkout functions to be used by any OperatorBase or Estimator
  //  the array checked out should be updated during evaluate
  //  object calling checkout is responsible for deleting the new array

  // checkout integer arrays
  template<int D>
  inline Array<TraceIntNew, D>* checkout_int(const std::string& name, int n1 = 1, int n2 = 0, int n3 = 0, int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_int"<<std::endl;
    return checkout_int<D>(default_domain, name, n1, n2, n3, n4);
  }
  template<int D>
  inline Array<TraceIntNew, D>* checkout_int(const std::string& domain,
                                          const std::string& name,
                                          int n1 = 1,
                                          int n2 = 0,
                                          int n3 = 0,
                                          int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_int"<<std::endl;
    TinyVector<int, DMAXNEW> shape(n1, n2, n3, n4);
    return int_samples.checkout_array<D>(domain, name, shape);
  }
  template<int D>
  inline Array<TraceIntNew, D>* checkout_int(const std::string& name,
                                          const ParticleSet& P,
                                          int n2 = 0,
                                          int n3 = 0,
                                          int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_int"<<std::endl;
    TinyVector<int, DMAXNEW> shape(P.getTotalNum(), n2, n3, n4);
    return int_samples.checkout_array<D>(P, name, shape);
  }


  // checkout real arrays
  template<int D>
  inline Array<TraceRealNew, D>* checkout_real(const std::string& name, int n1 = 1, int n2 = 0, int n3 = 0, int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_real"<<std::endl;
    return checkout_real<D>(default_domain, name, n1, n2, n3, n4);
  }
  template<int D>
  inline Array<TraceRealNew, D>* checkout_real(const std::string& domain,
                                            const std::string& name,
                                            int n1 = 1,
                                            int n2 = 0,
                                            int n3 = 0,
                                            int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_real"<<std::endl;
    TinyVector<int, DMAXNEW> shape(n1, n2, n3, n4);
    return real_samples.checkout_array<D>(domain, name, shape);
  }
  template<int D>
  inline Array<TraceRealNew, D>* checkout_real(const std::string& name,
                                            const ParticleSet& P,
                                            int n2 = 0,
                                            int n3 = 0,
                                            int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_real"<<std::endl;
    TinyVector<int, DMAXNEW> shape(P.getTotalNum(), n2, n3, n4);
    return real_samples.checkout_array<D>(P, name, shape);
  }


  // checkout complex arrays
  template<int D>
  inline Array<TraceCompNew, D>* checkout_complex(const std::string& name, int n1 = 1, int n2 = 0, int n3 = 0, int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_complex"<<std::endl;
    return checkout_complex<D>(default_domain, name, n1, n2, n3, n4);
  }
  template<int D>
  inline Array<TraceCompNew, D>* checkout_complex(const std::string& domain,
                                               const std::string& name,
                                               int n1 = 1,
                                               int n2 = 0,
                                               int n3 = 0,
                                               int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_complex"<<std::endl;
    TinyVector<int, DMAXNEW> shape(n1, n2, n3, n4);
    return comp_samples.checkout_array<D>(domain, name, shape);
  }
  template<int D>
  inline Array<TraceCompNew, D>* checkout_complex(const std::string& name,
                                               const ParticleSet& P,
                                               int n2 = 0,
                                               int n3 = 0,
                                               int n4 = 0)
  {
    app_log()<<"JTK: TraceCollector::checkout_complex"<<std::endl;
    TinyVector<int, DMAXNEW> shape(P.getTotalNum(), n2, n3, n4);
    return comp_samples.checkout_array<D>(P, name, shape);
  }


  inline void reset_buffers()
  {
    app_log() << "JTK: TraceCollector::reset_buffers"<<std::endl;
    if (writing_traces)
    {
      if (verbose)
        app_log() << "TraceCollector::reset_buffers " << std::endl;
      int_buffer.reset();
      real_buffer.reset();
    }
  }


  //store the full sample from a single walker step in buffers
  inline void buffer_sample(int current_step)
  {
    app_log() << "JTK: TraceCollector::buffer_sample "<<writing_traces<<" "<<current_step<<std::endl;
    if (writing_traces && current_step % throttle == 0)
    {
      if (verbose)
        app_log() << " TraceCollector::buffer_sample() " << std::endl;
      int_buffer.collect_sample();
      real_buffer.collect_sample();
    }
  }


  inline void startBlock(int nsteps)
  {
    app_log() << "JTK: TraceCollector::startBlock " << std::endl;
    if (verbose)
      app_log() << "TraceCollector::startBlock " << std::endl;
    reset_buffers();
  }


  inline void write_summary(std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    app_log() << std::endl;
    app_log() << pad << "TraceCollector (detailed summary)" << std::endl;
    app_log() << pad2 << "method_allows_traces    = " << method_allows_traces << std::endl;
    app_log() << pad2 << "streaming_traces        = " << streaming_traces << std::endl;
    app_log() << pad2 << "writing_traces          = " << writing_traces << std::endl;
    app_log() << pad2 << "default_domain          = " << default_domain << std::endl;
    int_buffer.write_summary(pad2);
    real_buffer.write_summary(pad2);
    app_log() << pad << "end TraceCollector" << std::endl;
  }


  inline void user_report(std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    std::string pad3 = pad2 + "  ";
    app_log() << std::endl;
    app_log() << pad << "TraceCollector report" << std::endl;
    request.report();
    app_log() << pad2 << "Type and domain breakdown of streaming quantities:" << std::endl;
    std::set<std::string>::iterator req;
    int_buffer.user_report(pad3);
    real_buffer.user_report(pad3);
    app_log() << std::endl;
  }

};






class TraceManagerNew
{
public:
  static double trace_tol;  // remove this

  TraceRequestNew request;

  std::string default_domain;
  bool method_allows_traces;
  bool streaming_traces;
  bool writing_traces;
  int throttle;
  bool verbose;
  std::string format;
  bool hdf_format;
  std::string file_root;
  Communicate* communicator;
  std::unique_ptr<hdf_archive> hdf_file;

  TraceManagerNew(Communicate* comm = 0) : verbose(false)
  {
    reset_permissions();
    communicator   = comm;
    throttle       = 1;
    format         = "hdf";
    default_domain = "scalars";
    request.set_scalar_domain(default_domain);
  }


  inline TraceManagerState getState()
  {
    TraceManagerState tms;
    tms.method_allows_traces = method_allows_traces;
    tms.request              = request;
    tms.streaming_traces     = streaming_traces;
    tms.writing_traces       = writing_traces;
    tms.throttle             = throttle;
    tms.verbose              = verbose;
    tms.default_domain       = default_domain;
    return tms;
  }

  inline TraceCollector* makeCollector()
  {
    app_log()<<"JTK: TraceManagerNew::makeCollector"<<std::endl;
    if (verbose)
      app_log() << "TraceManagerNew::makeCollector " << std::endl;
    TraceCollector* tc = new TraceCollector();
    tc->transfer_state_from(getState());
    tc->distribute();
    return tc;
  }


  inline void reset_permissions()
  {
    app_log()<<"JTK: TraceManagerNew::reset_permissions"<<std::endl;
    method_allows_traces = false;
    streaming_traces     = false;
    writing_traces       = false;
    verbose              = false;
    hdf_format           = false;
    request.reset();
  }


  inline void put(xmlNodePtr cur, bool allow_traces, std::string series_root)
  {
    app_log()<<"JTK: TraceManagerNew::put"<<std::endl;
    reset_permissions();
    method_allows_traces  = allow_traces;
    file_root             = series_root;
    bool traces_requested = cur != NULL;
    streaming_traces      = traces_requested && method_allows_traces;
    app_log()<<"JTK:    method_allows_traces "<<method_allows_traces<<std::endl;
    app_log()<<"JTK:    traces_requested "<<traces_requested<<std::endl;
    app_log()<<"JTK:    streaming_traces "<<streaming_traces<<std::endl;
    if (streaming_traces)
    {
      if (omp_get_thread_num() == 0)
      {
        app_log() << "\n  TraceManagerNew::put() " << std::endl;
        app_log() << "    traces requested          : " << traces_requested << std::endl;
        app_log() << "    method allows traces      : " << method_allows_traces << std::endl;
        app_log() << "    streaming traces          : " << streaming_traces << std::endl;
        app_log() << std::endl;
      }
      app_log()<<"JTK:    reading xml attributes "<<std::endl;
      //read trace attributes
      std::string writing         = "yes";
      std::string scalar          = "yes";
      std::string array           = "yes";
      std::string scalar_defaults = "yes";
      std::string array_defaults  = "yes";
      std::string verbose_write   = "no";
      OhmmsAttributeSet attrib;
      attrib.add(writing, "write");
      attrib.add(scalar, "scalar");
      attrib.add(array, "array");
      attrib.add(scalar_defaults, "scalar_defaults");
      attrib.add(array_defaults, "array_defaults");
      attrib.add(format, "format");
      attrib.add(throttle, "throttle");
      attrib.add(verbose_write, "verbose");
      attrib.add(array, "particle");                   //legacy
      attrib.add(array_defaults, "particle_defaults"); //legacy
      attrib.put(cur);
      writing_traces           = writing == "yes";
      bool scalars_on          = scalar == "yes";
      bool arrays_on           = array == "yes";
      bool use_scalar_defaults = scalar_defaults == "yes";
      bool use_array_defaults  = array_defaults == "yes";
      verbose                  = verbose_write == "yes";
      format                   = lowerCase(format);
      if (format == "hdf")
      {
        hdf_format = true;
      }
      else
      {
        APP_ABORT("TraceManagerNew::put " + format + " is not a valid file format for traces\n  valid options is: hdf");
      }

      //read scalar and array elements
      //  each requests that certain traces be computed
      std::set<std::string> scalar_requests;
      std::set<std::string> array_requests;
      xmlNodePtr element = cur->children;
      while (element != NULL)
      {
        std::string name((const char*)element->name);
        if (name == "scalar_traces")
        {
          std::string defaults = "no";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults, "defaults");
          eattrib.put(element);
          use_scalar_defaults = use_scalar_defaults && defaults == "yes";
          if (!use_scalar_defaults)
          {
            std::vector<std::string> scalar_list;
            putContent(scalar_list, element);
            scalar_requests.insert(scalar_list.begin(), scalar_list.end());
          }
        }
        else if (name == "array_traces" || name == "particle_traces")
        {
          std::string defaults = "no";
          OhmmsAttributeSet eattrib;
          eattrib.add(defaults, "defaults");
          eattrib.put(element);
          use_array_defaults = use_array_defaults && defaults == "yes";
          if (!use_array_defaults)
          {
            std::vector<std::string> array_list;
            putContent(array_list, element);
            array_requests.insert(array_list.begin(), array_list.end());
          }
        }
        else if (name != "text")
        {
          APP_ABORT("TraceManagerNew::put " + name +
                    " is not a valid sub-element of <trace/>\n  valid options are: scalar_traces, array_traces");
        }
        element = element->next;
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
      request.request_scalar(scalar_requests, writing_traces);
      request.request_array(array_requests, writing_traces);

    }
  }


  inline void check_clones(std::vector<TraceCollector*>& clones)
  {
    app_log() << "JTK: TraceManagerNew::check_clones "<<writing_traces<<" "<<clones.size() << std::endl;
    if (writing_traces && clones.size() > 0)
    {
      if (verbose)
        app_log() << "TraceManagerNew::check_clones" << std::endl;
      bool all_same = true;
      bool int_same;
      bool real_same;
      TraceCollector& ref = *clones[0];
      for (int i = 0; i < clones.size(); ++i)
      {
        TraceCollector& tm = *clones[i];
        int_same         = tm.int_buffer.same_as(ref.int_buffer);
        real_same        = tm.real_buffer.same_as(ref.real_buffer);
        all_same         = all_same && int_same && real_same;
      }
      if (!all_same)
      {
        for (int i = 0; i < clones.size(); ++i)
          clones[i]->write_summary();
        APP_ABORT("TraceManagerNew::check_clones  trace buffer widths of clones do not match\n  contiguous write is "
                  "impossible\n  this was first caused by clones contributing array traces from identical, but "
                  "differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the TraceManagerNew "
                  "summaries printed above");
      }
    }
  }


  //write buffered trace data to file
  inline void write_buffers(std::vector<TraceCollector*>& clones, int block)
  {
    app_log() << "JTK: TraceManagerNew::write_buffers "<< writing_traces<<std::endl;
    if (writing_traces)
    {
      //double tstart = MPI_Wtime();
      if (verbose)
        app_log() << "TraceManagerNew::write_buffers " << std::endl;
      if (hdf_format)
        write_buffers_hdf(clones);
    }
  }


  inline void open_file(std::vector<TraceCollector*>& clones)
  {
    app_log() << "JTK: TraceManagerNew::open_file "<< writing_traces<<std::endl;
    if (writing_traces)
    {
      if (verbose)
        app_log() << "TraceManagerNew::open_file " << std::endl;
      if (verbose)
        clones[0]->write_summary();
      if (hdf_format)
        open_hdf_file(clones);
    }
  }


  inline void close_file()
  {
    app_log() << "JTK: TraceManagerNew::close_file "<< writing_traces<<std::endl;
    if (writing_traces)
    {
      if (verbose)
        app_log() << "TraceManagerNew::close_file " << std::endl;
      if (hdf_format)
        close_hdf_file();
    }
  }


  inline void startRun(int blocks, std::vector<TraceCollector*>& clones)
  {
    app_log() << "JTK: TraceManagerNew::startRun " << std::endl;
    if (verbose)
      app_log() << "TraceManagerNew::startRun " << std::endl;
    check_clones(clones);
    open_file(clones);
  }


  inline void stopRun()
  {
    app_log() << "JTK: TraceManagerNew::stopRun " << std::endl;
    if (verbose)
      app_log() << "TraceManagerNew::stopRun " << std::endl;
    close_file();
  }

  //hdf file operations
  inline void open_hdf_file(std::vector<TraceCollector*>& clones)
  {
    app_log() << "JTK: TraceManagerNew::open_hdf_file " << std::endl;
    if (clones.size() == 0)
      APP_ABORT("TraceManagerNew::open_hdf_file  no trace clones exist, cannot open file");
    int nprocs = communicator->size();
    int rank   = communicator->rank();
    std::array<char, 32> ptoken;
    std::string file_name = file_root;
    if (nprocs > 1)
    {
      int length{0};
      if (nprocs > 10000)
        length = std::snprintf(ptoken.data(), ptoken.size(), ".p%05d", rank);
      else if (nprocs > 1000)
        length = std::snprintf(ptoken.data(), ptoken.size(), ".p%04d", rank);
      else
        length = std::snprintf(ptoken.data(), ptoken.size(), ".p%03d", rank);
      if (length < 0)
        throw std::runtime_error("Error generating filename");
      file_name.append(ptoken.data(), length);
    }
    file_name += ".traces.h5";
    if (verbose)
      app_log() << "TraceManagerNew::open_hdf_file  opening traces hdf file " << file_name << std::endl;
    hdf_file        = std::make_unique<hdf_archive>();
    bool successful = hdf_file->create(file_name);
    if (!successful)
      APP_ABORT("TraceManagerNew::open_hdf_file  failed to open hdf file " + file_name);
    // only clones have active buffers and associated data
    TraceCollector& tm = *clones[0];
    //tm.write_summary();
    tm.int_buffer.register_hdf_data(*hdf_file);
    tm.real_buffer.register_hdf_data(*hdf_file);
  }


  inline void write_buffers_hdf(std::vector<TraceCollector*>& clones)
  {
    app_log() << "JTK: TraceManagerNew::write_buffers_hdf " << std::endl;
    if (verbose)
      app_log() << "TraceManagerNew::write_buffers_hdf " << std::endl;
    TraceCollector& tm_lead = *clones[0];
    for (int ip = 0; ip < clones.size(); ++ip)
    {
      TraceCollector& tm = *clones[ip];
      tm.int_buffer.write_hdf(*hdf_file, tm_lead.int_buffer.hdf_file_pointer);
      tm.real_buffer.write_hdf(*hdf_file, tm_lead.real_buffer.hdf_file_pointer);
    }
  }

  inline void close_hdf_file() { 
    app_log() << "JTK: TraceManagerNew::close_hdf_file " << std::endl;
    hdf_file.reset(); 
  }
};




} // namespace qmcplusplus


#else

// no vacuous class: new tracemanager always supported

#endif


#endif
