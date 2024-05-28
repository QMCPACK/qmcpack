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
#include <unordered_set>



namespace qmcplusplus
{

const unsigned int WDMAX = 4;
using WTraceInt          = long;
using WTraceReal         = OHMMS_PRECISION_FULL;
using WTraceComp         = std::complex<WTraceReal>;
#ifndef QMC_COMPLEX
using WTracePsiVal = WTraceReal;
#else
using WTracePsiVal = WTraceComp;
#endif




struct WalkerQuantityInfo
{
  std::string name;
  size_t dimension;
  size_t size;
  size_t unit_size;
  TinyVector<size_t, WDMAX> shape;
  size_t buffer_start;
  size_t buffer_end;

  WalkerQuantityInfo(const std::string& name_,size_t unit_size_,size_t buffer_start_,size_t n1=1,size_t n2=0,size_t n3=0,size_t n4=0)
  {
    name         = name_;
    unit_size    = unit_size_;
    buffer_start = buffer_start_;
    shape[0]     = n1;
    shape[1]     = n2;
    shape[2]     = n3;
    shape[3]     = n4;

    dimension = 0;
    size      = 1;
    for (size_t d=0;d<WDMAX;++d)
      if (shape[d]>0)
      {
        dimension++;
        size *= shape[d];
      }
    buffer_end = buffer_start + size*unit_size;
  }
};



template<typename T>
struct WalkerTraceBuffer
{
  Array<T, 2> buffer;
  bool verbose;

  size_t walker_data_size;
  std::vector<WalkerQuantityInfo> quantity_info;

  bool first_collect;
  size_t quantity_pointer;
  size_t row_start_pointer;
  size_t row_end_pointer;

  //hdf variables
  std::string label;
  hsize_t dims[2];
  hsize_t hdf_file_pointer;

  WalkerTraceBuffer()
  {
    label             = "?";
    verbose           = false;
    first_collect     = true;
    walker_data_size  = 0;
    quantity_pointer  = 0;
    reset_buffer();
  }


  inline void reset_buffer()
  {
    app_log()<<"WalkerTraceBuffer("<<label<<")::reset_buffer"<<std::endl;
    buffer.resize(0, buffer.size(1));
    row_start_pointer = 0;
    row_end_pointer   = 0;    
  }


  inline void reset_collect()
  {
    if(quantity_pointer!=quantity_info.size())
      throw std::runtime_error("WalkerTraceBuffer quantity_pointer has not been moved through all quantities prior during collect() call.");
    first_collect = false;
    quantity_pointer = 0;
  }

  inline void reset_step()
  {
    if (verbose)
      app_log()<<"WalkerTraceBuffer("<<label<<")::reset_step"<<std::endl;
    row_start_pointer = row_end_pointer;
  }

  inline void reset_block() 
  { 
    if (verbose)
      app_log()<<"WalkerTraceBuffer("<<label<<")::reset_block"<<std::endl;
    reset_buffer();
  }


  inline bool same_as(WalkerTraceBuffer<T>& ref) 
  { 
    return buffer.size(1) == ref.buffer.size(1); 
  }


  inline void reset_rowsize(size_t row_size)
  {
    auto nrows = buffer.size(0);
    if(nrows==0)
      nrows++;
    if(nrows!=1)
      throw std::runtime_error("WalkerTraceBuffer::reset_rowsize  row_size (number of columns) should only be changed during growth of the first row.");
    auto buffer_old(buffer);
    buffer.resize(nrows,row_size);
    std::copy_n(buffer_old.data(), buffer_old.size(), buffer.data());
    if(buffer.size(0)!=1)
      throw std::runtime_error("WalkerTraceBuffer::reset_rowsize  buffer should contain only a single row upon completion.");
    if(buffer.size(1)!=row_size)
      throw std::runtime_error("WalkerTraceBuffer::reset_rowsize  length of buffer row should match the requested row_size following the reset/udpate.");
  }


  inline void make_new_row()
  {
    size_t nrows    = buffer.size(0);
    size_t row_size = buffer.size(1);
    if (row_size==0)
      throw std::runtime_error("WalkerTraceBuffer::make_new_row  Cannot make a new row of size zero.");
    nrows++;
    // resizing buffer(type Array) doesn't preserve data. Thus keep old data and copy over
    auto buffer_old(buffer);
    buffer.resize(nrows, row_size);
    std::copy_n(buffer_old.data(), buffer_old.size(), buffer.data());
  }


  inline void collect(const std::string& name, const T& value)
  {
    //if (verbose)
    //  app_log()<<"WalkerTraceBuffer("<<label<<")::collect"<<std::endl;
    size_t irow=0;
    if( first_collect )
    {
      WalkerQuantityInfo wqi_(name,1,walker_data_size);
      quantity_info.push_back(wqi_);
      walker_data_size = wqi_.buffer_end;
      reset_rowsize(walker_data_size);
    }
    else
    {
      if(quantity_pointer==0)
        make_new_row();
      irow = buffer.size(0)-1;
    }
    auto& wqi = quantity_info[quantity_pointer];
    buffer(irow,wqi.buffer_start) = value;
    quantity_pointer++;
  }


  template<unsigned D>
  inline void collect(const std::string& name, Array<T,D> arr)
  {
    //if (verbose)
    //  app_log()<<"WalkerTraceBuffer("<<label<<")::collect"<<std::endl;
    size_t n1 = arr.size(0);
    size_t n2,n3,n4;
    n2=n3=n4=0;
    if (D>4)
      throw std::runtime_error("WalkerTraceBuffer::collect  Only arrays up to dimension 4 are currently supported.");
    if (D>1) n2 = arr.size(1);
    if (D>2) n3 = arr.size(2);
    if (D>3) n4 = arr.size(3);
    size_t irow=0;
    if( first_collect )
    {
      WalkerQuantityInfo wqi_(name,1,walker_data_size,n1,n2,n3,n4);
      quantity_info.push_back(wqi_);
      walker_data_size = wqi_.buffer_end;
      reset_rowsize(walker_data_size);
    }
    else
    {
      if(quantity_pointer==0)
        make_new_row();
      irow = buffer.size(0)-1;
    }
    auto& wqi = quantity_info[quantity_pointer];
    auto& arr1d = arr.storage();
    for (size_t n=0;n<arr1d.size();++n)
      buffer(irow,wqi.buffer_start+n) = arr1d[n];
    quantity_pointer++;
  }


  template<unsigned D>
  inline void collect(const std::string& name, Array<std::complex<T>,D> arr)
  {
    size_t n1 = arr.size(0);
    size_t n2,n3,n4;
    n2=n3=n4=0;
    if (D>4)
      throw std::runtime_error("WalkerTraceBuffer::collect  Only arrays up to dimension 4 are currently supported.");
    if (D>1) n2 = arr.size(1);
    if (D>2) n3 = arr.size(2);
    if (D>3) n4 = arr.size(3);
    size_t irow=0;
    if( first_collect )
    {
      WalkerQuantityInfo wqi_(name,2,walker_data_size,n1,n2,n3,n4);
      quantity_info.push_back(wqi_);
      walker_data_size = wqi_.buffer_end;
      reset_rowsize(walker_data_size);
    }
    else
    {
      if(quantity_pointer==0)
        make_new_row();
      irow = buffer.size(0)-1;
    }
    auto& wqi = quantity_info[quantity_pointer];
    auto& arr1d = arr.storage();
    size_t n = 0;
    for (size_t i=0;i<arr1d.size();++i)
    {
      buffer(irow,wqi.buffer_start+n) = std::real(arr1d[n]); ++n;
      buffer(irow,wqi.buffer_start+n) = std::imag(arr1d[n]); ++n;
    }
    quantity_pointer++;
  }


  inline void add_row(Array<T, 2> other_buffer, size_t i)
  {
    app_log()<<"WalkerTraceBuffer("<<label<<")::add_row"<<std::endl;
    if(buffer.size(1)!=other_buffer.size(1))
      throw std::runtime_error("WalkerTraceBuffer::add_row  Row sizes must match.");
    make_new_row();
    size_t ib = buffer.size(0)-1;
    for(size_t j=0;j<buffer.size(1);++j)
      buffer(ib,j) = other_buffer(i,j);
  }


  inline void write_summary(std::string pad = "  ")
  {
    std::string pad2 = pad  + "  ";
    std::string pad3 = pad2 + "  ";
    app_log() << std::endl;
    app_log() << pad << "WalkerTraceBuffer(" << label << ")" << std::endl;
    app_log() << pad2 << "nrows       = " << buffer.size(0) << std::endl;
    app_log() << pad2 << "row_size    = " << buffer.size(1) << std::endl;
    for(size_t n=0;n<quantity_info.size();++n)
    {
      auto& wqi = quantity_info[n];
      app_log()<<pad2<<"quantity "<< n <<":  "<<wqi.dimension<<"  "<<wqi.size<<"  "<<wqi.unit_size<<"  "<<wqi.buffer_start<<"  "<<wqi.buffer_end<<" ("<<wqi.name<<")"<<std::endl;
    }
    app_log() << pad << "end WalkerTraceBuffer(" << label << ")" << std::endl;
  }

  inline void register_hdf_data(hdf_archive& f)
  {
    if (verbose)
      app_log()<<"WalkerTraceBuffer("<<label<<")::register_hdf_data"<<std::endl;
    auto& top = label;
    f.push(top);
    f.push("data_layout");
    for(auto& wqi: quantity_info)
    {
      f.push(wqi.name);
      f.write(wqi.dimension,    "dimension"  );
      f.write(wqi.shape,        "shape"      );
      f.write(wqi.size,         "size"       );
      f.write(wqi.unit_size,    "unit_size"  );
      f.write(wqi.buffer_start, "index_start");
      f.write(wqi.buffer_end,   "index_end"  );
      f.pop();
    }
    f.pop();
    f.pop();
    if (!f.open_groups())
      throw std::runtime_error("WalkerTraceBuffer(" + label +
                ")::register_hdf_data() some hdf groups are still open at the end of registration");
    hdf_file_pointer = 0;
  }


  inline void write_hdf(hdf_archive& f, hsize_t& file_pointer)
  {
    if (verbose)
      app_log()<<"WalkerTraceBuffer("<<label<<")::write_hdf "<<file_pointer<<"  "<<buffer.size(0)<<" "<<buffer.size(1)<<std::endl;
    auto& top = label;
    dims[0] = buffer.size(0);
    dims[1] = buffer.size(1);
    app_log()<<"    "<<buffer.dim()<<"  "<<dims[0]<<"  "<<dims[1]<<"  "<<buffer.size()<<"  "<<file_pointer<<std::endl;
    if (dims[0] > 0)
    {
      f.push(top);
      h5d_append(f.top(), "data", file_pointer, buffer.dim(), dims, buffer.data());
      f.pop();
    }
    f.flush();
    app_log()<<"    "<<buffer.dim()<<"  "<<dims[0]<<"  "<<dims[1]<<"  "<<buffer.size()<<"  "<<file_pointer<<std::endl;
  }
};





struct TraceManagerState 
{
  bool method_allows_traces;
  bool streaming_traces;
  bool writing_traces;
  int throttle;
  bool verbose;

  TraceManagerState()
  {
    reset_permissions();
    throttle       = 1;
    verbose        = false;
  }

  inline void reset_permissions()
  {
    method_allows_traces = false;
    streaming_traces     = false;
    writing_traces       = false;
    verbose              = false;
  }
};



template<class QT, class PT> class Walker;
class TrialWaveFunction;
class QMCHamiltonian;
using MCPWalker = Walker<QMCTraits, PtclOnLatticeTraits>;



class WalkerTraceCollector
{
public:

  TraceManagerState state;

  // new walker buffers
  std::unordered_set<std::string> properties_include;
  std::vector<size_t>             property_indices;
  int energy_index;

  std::vector<WTraceReal> energies;
  WalkerTraceBuffer<WTraceInt>  walker_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> walker_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> walker_particle_real_buffer;
  
  // temporary (contiguous) storage for awful ParticleAttrib<> quantities
  Array<WTraceReal  , 2> Rtmp;
  Array<WTraceReal  , 1> Stmp;
  Array<WTracePsiVal, 2> Gtmp;
  Array<WTracePsiVal, 1> Ltmp;


  WalkerTraceCollector()
    : properties_include{"R2Accepted","R2Proposed","LocalEnergy","LocalPotential","Kinetic","ElecElec","ElecIon","LocalECP","NonLocalECP"}
  {
    state.reset_permissions();

    // new walker buffers, etc
    energy_index = -1;
    energies.resize(0);
    walker_property_int_buffer.label  = "walker_property_int";
    walker_property_real_buffer.label = "walker_property_real";
    walker_particle_real_buffer.label = "walker_particle_real";
  }


  inline void set_state(const TraceManagerState& tms)
  {
    // JTK: rename this function as "setState"
    state = tms;
    walker_property_int_buffer.verbose  = state.verbose;
    walker_property_real_buffer.verbose = state.verbose;
    walker_particle_real_buffer.verbose = state.verbose;
  }


  inline void reset_step()
  {
    energies.resize(0);
  }


  inline void reset_buffers()
  {
    if (state.verbose)
      app_log() << "WalkerTraceCollector::reset_buffers"<<std::endl;
    walker_property_int_buffer.reset_buffer();
    walker_property_real_buffer.reset_buffer();
    walker_particle_real_buffer.reset_buffer();
  }


  inline void startBlock(int nsteps)
  {
    if (state.verbose)
      app_log() << "WalkerTraceCollector::startBlock " << std::endl;
    reset_buffers();
  }


  inline void write_summary(std::string pad = "  ")
  {
    std::string pad2 = pad + "  ";
    app_log() << std::endl;
    app_log() << pad << "TraceCollector (detailed summary)" << std::endl;
    app_log() << pad2 << "method_allows_traces    = " << state.method_allows_traces << std::endl;
    app_log() << pad2 << "streaming_traces        = " << state.streaming_traces << std::endl;
    app_log() << pad2 << "writing_traces          = " << state.writing_traces << std::endl;
    app_log() << pad << "end TraceCollector" << std::endl;
  }



  void collect(MCPWalker& walker, ParticleSet& pset, TrialWaveFunction& wfn, QMCHamiltonian& ham);

  void stopStep()
  {
    // find min/max/median walker and collect in buffers
  }

};





class WalkerTraceManager
{
public:

  std::string file_root;
  Communicate* communicator;
  std::unique_ptr<hdf_archive> hdf_file;
  bool registered_hdf;

  TraceManagerState state;

  // new walker buffers
  bool write_particle_data;
  bool write_min_data;
  bool write_max_data;
  bool write_med_data;

  WalkerTraceBuffer<WTraceInt>  wmin_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmin_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmin_particle_real_buffer;

  WalkerTraceBuffer<WTraceInt>  wmax_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmax_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmax_particle_real_buffer;

  WalkerTraceBuffer<WTraceInt>  wmed_property_int_buffer;
  WalkerTraceBuffer<WTraceReal> wmed_property_real_buffer;
  WalkerTraceBuffer<WTraceReal> wmed_particle_real_buffer;


  WalkerTraceManager(Communicate* comm = 0)
  {
    state.reset_permissions();
    communicator   = comm;

    registered_hdf = false;

    // new walker buffers, etc
    write_particle_data = false;
    write_min_data      = true;
    write_max_data      = true;
    write_med_data      = true;

    wmin_property_int_buffer.label  = "wmin_property_int";
    wmin_property_real_buffer.label = "wmin_property_real";
    wmin_particle_real_buffer.label = "wmin_particle_real";

    wmax_property_int_buffer.label  = "wmax_property_int";
    wmax_property_real_buffer.label = "wmax_property_real";
    wmax_particle_real_buffer.label = "wmax_particle_real";

    wmed_property_int_buffer.label  = "wmed_property_int";
    wmed_property_real_buffer.label = "wmed_property_real";
    wmed_particle_real_buffer.label = "wmed_particle_real";
  }


  inline TraceManagerState get_state()
  {
    return state;
  }

  inline WalkerTraceCollector* makeCollector()
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::makeCollector " << std::endl;
    WalkerTraceCollector* tc = new WalkerTraceCollector();
    tc->set_state(get_state());
    return tc;
  }


  inline void put(xmlNodePtr cur, bool allow_traces, std::string series_root)
  {
    if (state.verbose)
      app_log()<<"WalkerTraceManager::put"<<std::endl;
    state.reset_permissions();
    state.method_allows_traces  = allow_traces;
    file_root             = series_root;
    bool traces_requested = cur != NULL;
    state.streaming_traces      = traces_requested && state.method_allows_traces;
    if (state.streaming_traces)
    {
      if (omp_get_thread_num() == 0)
      {
        app_log() << "\n  WalkerTraceManager::put() " << std::endl;
        app_log() << "    traces requested          : " << traces_requested << std::endl;
        app_log() << "    method allows traces      : " << state.method_allows_traces << std::endl;
        app_log() << "    streaming traces          : " << state.streaming_traces << std::endl;
        app_log() << std::endl;
      }
      //read trace attributes
      std::string verbose_write = "no";
      OhmmsAttributeSet attrib;
      attrib.add(state.throttle, "throttle");
      attrib.add(verbose_write, "verbose");
      attrib.put(cur);
      state.verbose = verbose_write == "yes";
    }
  }


  inline void check_clones(std::vector<WalkerTraceCollector*>& clones)
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::check_clones" << std::endl;
    if (state.writing_traces && clones.size() > 0)
    {
      bool all_same = true;
      WalkerTraceCollector& ref = *clones[0];
      for (int i = 0; i < clones.size(); ++i)
      {
        WalkerTraceCollector& tm = *clones[i];
        all_same &= tm.walker_property_int_buffer.same_as(ref.walker_property_int_buffer);
        all_same &= tm.walker_property_real_buffer.same_as(ref.walker_property_real_buffer);
        all_same &= tm.walker_particle_real_buffer.same_as(ref.walker_particle_real_buffer);
      }
      if (!all_same)
      {
        for (int i = 0; i < clones.size(); ++i)
          clones[i]->write_summary();
        throw std::runtime_error("WalkerTraceManager::check_clones  trace buffer widths of clones do not match\n  contiguous write is "
                  "impossible\n  this was first caused by clones contributing array traces from identical, but "
                  "differently named, particlesets such as e, e2, e3 ... (fixed)\n  please check the WalkerTraceManager "
                  "summaries printed above");
      }
    }
  }


  //write buffered trace data to file
  inline void write_buffers(std::vector<WalkerTraceCollector*>& clones, int block)
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::write_buffers "<<std::endl;
    if (state.writing_traces)
    {
      write_buffers_hdf(clones);
    }
  }


  inline void open_file(std::vector<WalkerTraceCollector*>& clones)
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::open_file "<<std::endl;
    if (state.writing_traces)
    {
      if (state.verbose)
        clones[0]->write_summary();
      open_hdf_file(clones);
    }
  }


  inline void close_file()
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::close_file " << std::endl;
    if (state.writing_traces)
      close_hdf_file();
  }


  inline void startRun(int blocks, std::vector<WalkerTraceCollector*>& clones)
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::startRun " << std::endl;
    check_clones(clones);
    open_file(clones);
  }


  inline void stopRun()
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::stopRun " << std::endl;
    close_file();
  }

  //hdf file operations
  inline void open_hdf_file(std::vector<WalkerTraceCollector*>& clones)
  {
    if (state.verbose) 
      app_log() << "WalkerTraceManager::open_hdf_file " << std::endl;
    if (clones.size() == 0)
      throw std::runtime_error("WalkerTraceManager::open_hdf_file  no trace clones exist, cannot open file");
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
    file_name += ".wtraces.h5";
    if (state.verbose)
      app_log() << "WalkerTraceManager::open_hdf_file  opening traces hdf file " << file_name << std::endl;
    hdf_file        = std::make_unique<hdf_archive>();
    bool successful = hdf_file->create(file_name);
    if (!successful)
      throw std::runtime_error("WalkerTraceManager::open_hdf_file  failed to open hdf file " + file_name);
    // only clones have active buffers and associated data
    WalkerTraceCollector& tm = *clones[0];
  }


  inline void write_buffers_hdf(std::vector<WalkerTraceCollector*>& clones)
  {
    if (state.verbose)
      app_log() << "WalkerTraceManager::write_buffers_hdf " << std::endl;
    WalkerTraceCollector& tm_lead = *clones[0];
    if(!registered_hdf)
    {
      tm_lead.walker_property_int_buffer.register_hdf_data(*hdf_file);
      tm_lead.walker_property_real_buffer.register_hdf_data(*hdf_file);
      tm_lead.walker_particle_real_buffer.register_hdf_data(*hdf_file);

      registered_hdf = true;
    }
    for (int ip = 0; ip < clones.size(); ++ip)
    {
      WalkerTraceCollector& tm = *clones[ip];
      tm.walker_property_int_buffer.write_hdf(*hdf_file, tm_lead.walker_property_int_buffer.hdf_file_pointer);
      tm.walker_property_real_buffer.write_hdf(*hdf_file, tm_lead.walker_property_real_buffer.hdf_file_pointer);
      tm.walker_particle_real_buffer.write_hdf(*hdf_file, tm_lead.walker_particle_real_buffer.hdf_file_pointer);
    }
  }

  inline void close_hdf_file() { 
    if (state.verbose)
      app_log() << "WalkerTraceManager::close_hdf_file " << std::endl;
    hdf_file.reset(); 
  }

};




} // namespace qmcplusplus



#endif
