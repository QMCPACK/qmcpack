//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <Configuration.h>
#include <Particle/ParticleSet.h>
#include <Particle/FastParticleOperators.h>
#include <ParticleIO/ParticleLayoutIO.h>
//#include <OhmmsPETE/OhmmsVector.h>
//#include <OhmmsPETE/OhmmsMatrix.h>
//#include <algorithm>
#include <io/hdf_archive.h>
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/OpenMP.h"
#include <Utilities/IteratorUtility.h>
#include <Numerics/VectorViewer.h>
#include <Utilities/unit_conversion.h>
#include <QMCApp/ParticleSetPool.h>

// design 
//   read input file
//   determine process partitioning
//   initialize quantities
//   for s in series
//     for g in groups
//       if master
//         open output hdf file
//       #end if
//       for p in processes
//         read shared quantities
//         determine block count
//         for q in quantities
//           q.reset(nrows,nblocks) reset block count
//           q.load_trace()         read quantity q
//           q.analyze()            analyze/accumulate quantity q
//           q.clear_trace()        free trace data
//         end for
//       end for
//       reduce(min) min block count
//       resize buffer
//       for b in blocks
//         for q in quantities
//           q.pack(b,buffer);
//         end for
//         gather buffers
//         if master
//           for q in quantities
//             q.write(master_buffer)
//           end for
//         end if
//       end for
//       if master
//         close output hdf file
//       end if
//     end for
//   end for  
//  
//  
//  Input
//  QuantityTrace id,weight;
//  QuantityAnalyzer
//    QuantityTrace trace_data;
//    int buffer_index;
//  
//  
//  
//  


/** xml input file for density
<?xml version="1.0"?>
<qmcsystem>
  <simulationcell>
    <parameter name="scale">30 </parameter>
    <parameter name="lattice">
      1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0
    </parameter>
    <parameter name="bconds">p p p</parameter>
  </simulationcell>
  <density file_root="XYZ" first="0" last="1">
  </density>
</qmcsystem>
*/

using namespace qmcplusplus;

typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
typedef ParticleSet::ParticlePos_t ParticlePos_t;
typedef ParticleSet::Tensor_t Tensor_t;

enum{DIM=OHMMS_DIM};
typedef double real_t;
typedef unsigned long mem_t;
typedef TinyVector<real_t,DIM> Pos_t;


const std::string & none_type_string = "none";
const std::string & int_type_string  = "int";
const std::string & real_type_string = "real";

template<typename T>
const std::string& get_type_string()
{
  return none_type_string;
}

template<>
const std::string& get_type_string<int>()
{
  return int_type_string;
}

template<>
const std::string& get_type_string<real_t>()
{
  return real_type_string;
}


void zfill(int number,std::string& str,int digits)
{
  char token[32];
  if(digits==6)
    sprintf(token,"%06d",number);
  else if(digits==5)
    sprintf(token,"%05d",number);
  else if(digits==4)
    sprintf(token,"%04d",number);
  else if(digits==3)
    sprintf(token,"%03d",number);
  else
    APP_ABORT("zfill: only digits 3-6 are supported");
  str = token;
}


template<typename T>
T min(std::vector<T>&v)
{
  if(v.size()>0)
  {
    T m = v[0];
    for(int i=1;i<v.size();++i)
      m= std::min(m,v[i]);
    return m;
  }
  else
  {
    APP_ABORT("min(vector): vector is size 0");
    return T(0);
  }
}

template<typename T>
T max(std::vector<T>&v)
{
  if(v.size()>0)
  {
    T m = v[0];
    for(int i=1;i<v.size();++i)
      m= std::max(m,v[i]);
    return m;
  }
  else
  {
    APP_ABORT("max(vector): vector is size 0");
    return T(0);
  }
}


template<typename T>
T sum(std::vector<T>&v)
{
  T s = 0.0;
  for(int i=0;i<v.size();++i)
    s+=v[i];
  return s;
}

template<typename T>
T sum2(std::vector<T>&v)
{
  T s = 0.0;
  for(int i=0;i<v.size();++i)
    s+=v[i]*v[i];
  return s;
}


template<typename T>
T mean(std::vector<T>&v) 
{
  return sum(v)/v.size();
}

template<typename T>
T mean2(std::vector<T>&v) 
{
  return sum2(v)/v.size();
}


template<typename T>
void autocorr_stats(std::vector<T>& block_data,T& avg,T& var,T& error,T& blocks_eff)
{
  avg = mean(block_data);
  var = mean2(block_data)-avg*avg;
  int blocks = block_data.size();
  real_t kappa=0.0;
  if(var>1e-10)
  {
    real_t corr = .5;
    for(int i=0;i<blocks-1;++i)
    {
      if(corr<0.0)
        break;
      kappa += 2*corr;
      corr = 0.0;
      for(int j=0;j<blocks-i;++j)
        corr += (block_data[j]-avg)*(block_data[i+j]-avg);
      corr /= (var*(blocks-i));
    }
  }
  blocks_eff = blocks/std::max(kappa,1.0);
  error = sqrt(var/blocks_eff);
}


template<typename T>
void autocorr_stats(std::vector<T>& block_data,T& blocks_eff)
{
  T avg,var,error;
  autocorr_stats(block_data,avg,var,error,blocks_eff);
}


struct block_props
{
  int i1;
  int i2;
  real_t w;
  block_props(int n1,int n2) : i1(n1),i2(n2) { w=0.0;}
  block_props(int n1,int n2,real_t wt) : i1(n1),i2(n2),w(wt) { }
};



struct Input
{
  // qmcsystem
  ParticleLayout_t cell;
  ParticleSet* Peln;
  ParticleSet* Pion;
  // run
  std::string path;
  std::string prefix;
  int tot_procs;
  int tot_groups;
  int tot_series;
  std::vector<int> steps;
  std::vector<int> blocks;
  int proc_digits;
  std::string elns;
  std::string ions;
  // controls
  int blocks_exclude;
  bool average_groups;
  bool verbose;
  std::vector<int> series;
  std::vector<int> groups;
  std::vector<int> procs;
  // quantities
  xmlNodePtr quantities;

  //derived parameters
  int nspecies;
  std::vector<int> species_size;
  std::vector<std::string> species_name;
  int nparticles;

  Libxml2Document doc;
  Communicate* comm;
  ParticleSetPool ptclPool;

  Input(Communicate* c)
    : ptclPool(c)
  {
  }

  void read(char* infile)
  {
    // qmcsystem
    // run
    path        = ".";
    prefix      = "";
    tot_procs   = -1;
    tot_groups  = 0;
    tot_series  = -1;
    proc_digits = 3;
    elns        = "e";
    ions        = "ion0";
    // controls
    blocks_exclude = 0;
    verbose = false;
    average_groups = true;
    std::vector<int> series_range;
    std::vector<int> group_range;
    std::vector<int> proc_range;
    std::string vstr  = "no";
    std::string agstr = "yes";
    // quantities
    quantities = NULL;

    bool success    = doc.parse(infile);
    xmlNodePtr elem = doc.getRoot()->children;
    while(elem != NULL)
    {
      std::string ename((const char*)elem->name);
      xmlNodePtr cur=elem->children;
      if(ename=="qmcsystem")
      {
        processP(elem);
        while(cur!=NULL)
        {
          std::string cname((const char*)cur->name);
          if(cname=="simulationcell")
          {
            LatticeParser a(cell);
            a.put(cur);
          }
          //else if(cname=="electrons")
          //{
          //  OhmmsAttributeSet a;
          //  a.add(nup,"nup");
          //  a.add(ndown,"ndown");
          //  a.put(cur);
          //}
          //else if(cname!="text" && cname!="comment")
          //{
          //  APP_ABORT(cname+" is not an element of <qmcsystem>");
          //}
          cur=cur->next;
        }
      }
      else if(ename=="run")
      {
        while(cur!=NULL)
        {
          std::string cname((const char*)cur->name);
          if(cname=="parameter")
          {
            std::string pname((const char*)(xmlGetProp(cur,(const xmlChar*)"name")));
            if(pname=="path")
              putContent(path,cur);
            else if(pname=="prefix")
              putContent(prefix,cur);
            else if(pname=="tot_series")
              putContent(tot_series,cur);
            else if(pname=="tot_groups")
              putContent(tot_groups,cur);
            else if(pname=="tot_procs")
              putContent(tot_procs,cur);
            else if(pname=="steps")
              putContent(steps,cur);
            else if(pname=="blocks")
              putContent(blocks,cur);
            else if(pname=="proc_digits")
              putContent(proc_digits,cur);
            else if(pname=="elns")
              putContent(elns,cur);
            else if(pname=="ions")
              putContent(ions,cur);
            else
              APP_ABORT(pname+" is not a parameter of <run>");
          }
          else if(cname!="text" && cname!="comment")
          {
            APP_ABORT(cname+" is not an element of <run>");
          }
          cur=cur->next;
        }
      }
      else if(ename=="controls")
      {
        while(cur!=NULL)
        {
          std::string cname((const char*)cur->name);
          if(cname=="parameter")
          {
            std::string pname((const char*)(xmlGetProp(cur,(const xmlChar*)"name")));
            if(pname=="series")
              putContent(series,cur);
            else if(pname=="groups")
              putContent(groups,cur);
            else if(pname=="procs")
              putContent(procs,cur);
            else if(pname=="series_range")
              putContent(series_range,cur);
            else if(pname=="group_range")
              putContent(group_range,cur);
            else if(pname=="proc_range")
              putContent(proc_range,cur);
            else if(pname=="verbose")
            {
              putContent(vstr,cur);
              verbose = vstr=="yes";
            }
            else if(pname=="blocks_exclude")
              putContent(blocks_exclude,cur);
            else if(pname=="average_groups")
            {
              putContent(agstr,cur);
              average_groups = agstr=="yes";
            }
            else
              APP_ABORT(pname+" is not a parameter of <controls>");
          }
          else if(cname!="text" && cname!="comment")
          {
            APP_ABORT(cname+" is not an element of <controls>");
          }
          cur=cur->next;
        }
      }
      else if(ename=="quantities")
      {
        quantities = elem;
      }
      else if(ename!="text" && ename!="comment")
      {
        APP_ABORT(ename+" is not an element of <traceanalysis>");
      }
      elem=elem->next;
    }
    set_range("series",series,series_range);
    set_range("groups",groups,group_range );
    set_range("procs" ,procs ,proc_range  );
    if(series.size()==0)
      for(int s=0;s<tot_series;++s)
        series.push_back(s);
    if(procs.size()==0)
      for(int p=0;p<tot_procs;++p)
        procs.push_back(p);
    if(groups.size()==0)
    {
      if(tot_groups>0)
      {
        for(int p=0;p<tot_groups;++p)
          groups.push_back(p);
      }
      else
      {
        groups.push_back(-1);
      }
    }
    Peln = get_particleset(elns);
    Pion = get_particleset(ions);
    if(tot_series<1)
      APP_ABORT("tot_series must be provided");
    if(tot_procs<1)
      APP_ABORT("tot_procs must be provided");
    if(steps.size()!=tot_series)
      APP_ABORT("steps must be given for each series");
    if(blocks.size()!=tot_series)
      APP_ABORT("steps must be given for each series");
    if(quantities==NULL)
      APP_ABORT("quantities must be provided");
    if(Peln==NULL)
      APP_ABORT("electron particleset "+elns+" does not exist");
    if(Pion==NULL)
      APP_ABORT("electron particleset "+ions+" does not exist");

    SpeciesSet& species = Peln->getSpeciesSet();
    nspecies = species.size();
    int isize  = species.addAttribute("membersize");
    if(isize==species.numAttributes())
      APP_ABORT("species set does not have the required attribute 'membersize'");
    for(int s=0;s<nspecies;++s)
      species_size.push_back(species(isize,s));
    for(int s=0;s<nspecies;++s)
      species_name.push_back(species.speciesName[s]);
    nparticles = 0;
    for(int s=0;s<nspecies;++s)
      nparticles += species_size[s];
  }

  void processP(xmlNodePtr qmcsystem)
  {
    xmlNodePtr cur = qmcsystem->children;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == "simulationcell")
        ptclPool.putLattice(cur);
      else if(cname == "particleset")
        ptclPool.put(cur);
      cur=cur->next;
    }
  }

  ParticleSet* get_particleset(const std::string& psname)
  {
    std::map<std::string,ParticleSet*>& pool = ptclPool.getPool();
    if(pool.find(psname)==pool.end())
      return NULL;
    else
      return pool[psname];
  }

  void set_range(const std::string& name,std::vector<int>& vec,std::vector<int>& range)
  {
    if(range.size()==0)
    {
      return;
    }
    else if(range.size()!=2)
    {
      APP_ABORT("range for "+name+" must have only 2 entries");
    }
    for(int i=range[0];i<range[1]+1;++i)
      vec.push_back(i);
  }

  void get_range(std::vector<int>& vec,int& vmin,int& vmax)
  {
    vmin = min(vec);
    vmax = max(vec);
  }

  void report()
  {
    app_log()<< std::endl;
    app_log()<<"Trace analyzer input:"<< std::endl;
    app_log()<<"  path           = "<< path << std::endl;
    app_log()<<"  prefix         = "<< prefix << std::endl;
    app_log()<<"  tot_procs      = "<< tot_procs << std::endl;
    app_log()<<"  tot_groups     = "<< tot_groups << std::endl;
    app_log()<<"  tot_series     = "<< tot_series << std::endl;
    app_log()<<"  steps          = ";
    for(int i=0;i<steps.size();++i)
      app_log()<<steps[i]<<" ";
    app_log()<< std::endl;
    app_log()<<"  blocks         = ";
    for(int i=0;i<blocks.size();++i)
      app_log()<<blocks[i]<<" ";
    app_log()<< std::endl;
    app_log()<<"  proc_digits    = "<< proc_digits << std::endl;
    int smin,smax;
    get_range(series,smin,smax);
    app_log()<<"  series range   = "<<smin<<" "<<smax<< std::endl;
    if(tot_groups>0)
    {
      get_range(groups,smin,smax);
      app_log()<<"  group range    = "<<smin<<" "<<smax<< std::endl;
    }
    get_range(procs,smin,smax);
    app_log()<<"  proc range     = "<<smin<<" "<<smax<< std::endl;
    app_log()<<"  average_groups = "<<average_groups<< std::endl;
    app_log()<<"  verbose        = "<<verbose<< std::endl;
  }
};


struct QuantityTraceBase
{
  static bool verbose;
  static real_t read_time;
  std::string path;
  int nrows;
  int ncols;
  Timer time;

  virtual void clear()=0;
  virtual bool read(hdf_archive& hin,const std::string& inpath)=0;
  virtual mem_t memory()=0;
};
bool QuantityTraceBase::verbose = false;
real_t QuantityTraceBase::read_time = 0.0;


template<typename T>
struct QuantityTrace : public QuantityTraceBase
{
  std::vector<T> trace;

  QuantityTrace()
  {
    clear();
  }

  void clear()
  {
    path  = "";
    nrows = 0;
    ncols = 0;
    trace.clear();
  }


  bool read(hdf_archive& hin,const std::string& inpath)
  {
    clear();
    path = inpath;
    if(verbose)
      app_log()<<"    reading "<<path<< std::endl;
    time.restart();
    hsize_t dims[2];
    int blockoffset=0;
    const std::string& pushpath = "/"+get_type_string<T>()+"_data";
    bool success=hin.push(pushpath,false);
    if(success)
    {
      hid_t gid=hin.top();
      bool gotit=get_space(gid,"traces",2,dims);
      //if(gotit)
      {
        int col_start;
        int col_end;
        success&=hin.read(col_start,"layout/"+path+"/row_start");
        success&=hin.read(col_end,"layout/"+path+"/row_end");
        nrows = dims[0];
        ncols = col_end-col_start;
        trace.resize(nrows*ncols);
        TinyVector<int,2> gcounts(dims[0],dims[1]);
        TinyVector<int,2> counts(nrows,ncols);
        TinyVector<int,2> offsets(blockoffset,col_start);
        hyperslab_proxy<std::vector<T>,2> trace_slab(trace,gcounts,counts,offsets);
        success&=hin.read(trace_slab,"traces");
      }
    }
    hin.pop();
    if(verbose)
    {
      app_log()<<"      successful = "<<success<< std::endl;
      app_log()<<"      wall time  = "<<time.elapsed()<< std::endl;
      app_log()<<"      mbytes     = "<<(memory()>>20)<< std::endl;
      app_log()<<"      nrows      = "<<nrows<< std::endl;
      app_log()<<"      ncols      = "<<ncols<< std::endl;
    }
    return success;
  }

  mem_t memory()
  {
    mem_t mem=trace.size();
    mem *= sizeof(T);
    return mem;
  }
};


struct QuantityAnalyzer
{
  std::string name;
  Input& input;
  int series_index;
  int group_index;
  std::string label;
  int pack_start;
  size_t blocks;
  size_t block_size;
  std::vector<real_t> block_weights;
  std::vector<real_t> block_buffer;
  std::vector<real_t> data_mean;
  std::vector<real_t> data_error;
  std::vector<real_t> stat_blocks;

  QuantityAnalyzer( std::string& n,Input& i)
    : name(n),input(i)
  {
    blocks     = 0;
    block_size = 0;
  }

  virtual ~QuantityAnalyzer() {}

  void init(xmlNodePtr xml,int s)
  {
    series_index = s;
    blocks     = input.blocks[series_index];
    block_size = 0;
    label      = "";

    OhmmsAttributeSet attrib;
    attrib.add(label,"label");
    attrib.put(xml);

    put(xml);

    if(block_size==0)
      APP_ABORT("put must set block_size, this is a developer error");
    block_weights.resize(blocks);
    block_buffer.resize(blocks*block_size);
    data_mean.resize(block_size);
    data_error.resize(block_size);
    stat_blocks.resize(blocks);
    if(input.average_groups)
    {
      fill(block_weights.begin(),block_weights.end(),0.0);
      fill(block_buffer.begin(),block_buffer.end(),0.0);
    }
  }

  mem_t memory()
  {
    mem_t mem = block_buffer.size();
    mem*=sizeof(real_t);
    return mem;
  }

  void set_group(int g)
  {
    group_index = g;
    if(!input.average_groups)
    {
      fill(block_weights.begin(),block_weights.end(),0.0);
      fill(block_buffer.begin(),block_buffer.end(),0.0);
    }
  }

  void get_pack_start(int& bsize)
  {
    pack_start = bsize;
    bsize += block_size;
  }

  void pack(int b,std::vector<real_t>& buffer)
  {
    int boffset = b*block_size;
    for(int n=0;n<block_size;++n)
      buffer[pack_start+n] = block_buffer[boffset+n];
  }

  void unpack(int b,std::vector<real_t>& buffer)
  {
    block_weights[b] = buffer[0];
    int boffset = b*block_size;
    for(int n=0;n<block_size;++n)
      block_buffer[boffset+n] = buffer[pack_start+n]; 
  }

  void analyze()
  {
    int blocks_exclude = input.blocks_exclude;
    int blocks_include = blocks-blocks_exclude;
    //normalize block data
    for(int b=0,boffset=0;b<blocks;++b,boffset+=block_size)
    {
      real_t norm = 1.0/block_weights[b];
      for(int n=0;n<block_size;++n)
        block_buffer[boffset+n] *= norm;
    }
    //find the mean and variance
    fill(data_mean.begin() ,data_mean.end() ,0.0);
    fill(data_error.begin(),data_error.end(),0.0);
    for(int b=blocks_exclude;b<blocks;++b)
    {
      int boffset = b*block_size;
      for(int n=0;n<block_size;++n)
      {
        real_t val = block_buffer[boffset+n];
        data_mean[n]  += val;
        data_error[n] += val*val;
      }
    }
    for(int n=0;n<block_size;++n)
      data_mean[n] /= blocks_include;
    for(int n=0;n<block_size;++n)
      data_error[n] = data_error[n]/blocks_include-data_mean[n]*data_mean[n]; //variance
    //compute blocks_eff where the data is largest
    int    nmax = -1;
    real_t dmax = -1.0;
    for(int n=0;n<block_size;++n)
    {
      real_t d = std::abs(data_mean[n]);
      if(d>dmax)
      {
        dmax=d;
        nmax=n;
      }
    }
    stat_blocks.resize(blocks_include);
    for(int b=blocks_exclude;b<blocks;++b)
    {
      int boffset = b*block_size;
      stat_blocks[b] = block_buffer[boffset+nmax];
    }
    real_t blocks_eff;
    autocorr_stats(stat_blocks,blocks_eff);
    //compute the error
    for(int n=0;n<block_size;++n)
      data_error[n] = sqrt(data_error[n]/blocks_eff); //error
  }

  void get_output_file_prefix( std::string& ofprefix,int s,int g)
  {
    int series = input.series[s];
    int group  = input.groups[g];
    std::string stoken;
    std::string gtoken;
    zfill(series,stoken,3);
    ofprefix = input.path+"/"+input.prefix;
    if(group!=-1 && !input.average_groups)
    {
      zfill(group,gtoken,3);
      ofprefix += ".g"+gtoken;
    }
    ofprefix += ".s"+stoken+"."+name;
    if(label!="")
      ofprefix += "_"+label;
  }
  

  virtual void put(xmlNodePtr qxml)=0;  
  virtual void load_traces(hdf_archive& hin)=0;
  virtual void accumulate(QuantityTrace<real_t>& weights,std::vector<block_props>& block_properties)=0;
  virtual void clear_traces()=0;
  virtual void write()=0;
};


struct SpinDensityAnalyzer : public QuantityAnalyzer
{
  std::string format;
  int npoints;
  TinyVector<int,DIM> grid;
  TinyVector<int,DIM> gdims;
  QuantityTrace<real_t> positions;
  Pos_t corner;

  int nspecies;
  int nparticles;

  SpinDensityAnalyzer( std::string& n,Input& i)
    : QuantityAnalyzer(n,i)
  {
    nspecies   = input.nspecies;
    nparticles = input.nparticles;
  }

  void put(xmlNodePtr qxml)
  {
    format = "xsf";
    bool has_grid = false;
    xmlNodePtr cur=qxml->children;
    while(cur!=NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname=="parameter")
      {
        std::string pname((const char*)(xmlGetProp(cur,(const xmlChar*)"name")));
        if(pname=="grid")
        {
          has_grid = true;
          putContent(grid,cur);
        }
        else if(pname=="format")
          putContent(format,cur);
      }
      cur=cur->next;
    }
    if(!has_grid)
      APP_ABORT("spindensity must have parameter grid");

    if(format!="xsf")
      APP_ABORT("SpinDensityAnalyzer::put\n  file format "+format+" is unknown\n  valid options are: xsf");


    npoints = 1;
    for(int d=0;d<DIM;++d)
      npoints *= grid[d];

    gdims[0] = npoints/grid[0];
    for(int d=1;d<DIM;++d)
      gdims[d] = gdims[d-1]/grid[d];

    for(int d=0;d<DIM;++d)
      corner[d] = .5/grid[d];
    corner = input.cell.toCart(corner);

    block_size = nspecies*npoints;
  }

  void load_traces(hdf_archive& hin)
  {
    positions.read(hin,"e/position");
    if(positions.ncols/DIM!=input.nparticles)
      APP_ABORT("SpinDensityAnalyzer::load_traces number of particles in traces file and input file differ");
  }

  void accumulate(QuantityTrace<real_t>& weights,std::vector<block_props>& block_properties)
  {
    ParticleLayout_t& cell = input.cell;
    ParticlePos_t R(nparticles); //cartesian coordiates
    ParticlePos_t Ru(nparticles);//unit coordinates [0,1)
    const int psize = DIM*nparticles;

    for(int b=0;b<blocks;++b)
    {
      const block_props& bp = block_properties[b];
      block_weights[b] += bp.w;
      int boffset = b*block_size;
      int poffset = bp.i1*psize;
      for(int i=bp.i1;i<bp.i2;++i,poffset+=psize)
      {
        std::copy(positions.trace.data()+poffset,positions.trace.data()+poffset+psize, &(R[0][0]));
        ApplyBConds<ParticlePos_t,Tensor_t,3,false>::Cart2Unit(R,cell.G,Ru,0,nparticles);
        real_t w = weights.trace[i];
        int p=0;
        int offset = boffset;
        for(int s=0;s<input.nspecies;++s,offset+=npoints)
          for(int ps=0;ps<input.species_size[s];++ps,++p)
          {
            const Pos_t& u = Ru[p];
            int point=offset;
            for(int d=0;d<DIM;++d)
              point += gdims[d]*((int)(grid[d]*(u[d]-floor(u[d]))));
            block_buffer[point] += w;
          }
      }
      //real_t norm = 1.0/bp.w;
      //for(int n=0;n<block_size;++n)
      //  block_buffer[boffset+n] *= norm;
    }

    //APP_ABORT("sda::analyze");
  }

  void clear_traces()
  {
    positions.clear();
  }

  void write()
  {
    std::string ofprefix;
    get_output_file_prefix(ofprefix,series_index,group_index);
    
    std::vector<real_t> dwrite(npoints);
    std::vector<real_t> dmean(npoints,0.0);
    std::vector<real_t> derror(npoints,0.0);
    for(int s=0,offset=0;s<input.nspecies;++s,offset+=npoints)
    {
      std::string sfprefix = ofprefix+"_"+input.species_name[s];

      VectorViewer<real_t> dm(&data_mean[offset],npoints);
      VectorViewer<real_t> de(&data_error[offset],npoints);
      VectorViewer<real_t> dwr(&dwrite[0],npoints);

      write_densities(sfprefix,dm,de,dwr);

      for(int p=0;p<npoints;++p)
        dmean[p] += dm[p];
      for(int p=0;p<npoints;++p)
        derror[p] += de[p]*de[p];
    }
    for(int p=0;p<npoints;++p)
      derror[p] += sqrt(derror[p]);

    write_densities(ofprefix,dmean,derror,dwrite);
    

    //vector<real_t> dmn;
    //dmn.resize(block_size);
    //fill(dmn.begin(),dmn.end(),0.0);
    //int boffset = 0;
    //for(int b=0;b<blocks;++b,boffset+=block_size)
    //  for(int n=0;n<block_size;++n)
    //    dmn[n] += block_buffer[boffset+n];
    //for(int n=0;n<block_size;++n)
    //  dmn[n] /= blocks;
    //
    //app_log()<<"\nsums of SpinDensity"<< std::endl;
    ////boffset = 0;
    ////for(int b=0;b<blocks;++b,boffset+=block_size)
    ////{
    ////  real_t bsum = 0.0;
    ////  for(int n=0;n<block_size;++n)
    ////    bsum += block_buffer[boffset+n];
    ////  app_log()<<"  block sum "<<b<<" "<<bsum<< std::endl;
    ////}
    //for(int s=0,n=0;s<input.nspecies;++s)
    //{
    //  app_log()<<"  spin "<<s<< std::endl;
    //  for(int p=0;p<npoints;++p,++n)
    //    fprintf(stdout,"    %d  %16.12f  %16.12f  %16.12f\n",n,dmn[n],data_mean[n],data_error[n]);
    //}
  }

  template<typename VEC>
  void write_densities(const std::string& ofprefix,VEC& dmean,VEC& derror,VEC& dwrite)
  {
    write_density(ofprefix+"_mean",dmean);
    write_density(ofprefix+"_err" ,derror);

    for(int p=0;p<npoints;++p)
      dwrite[p] = dmean[p]+derror[p];
    write_density(ofprefix+"_mean+err",dwrite);

    for(int p=0;p<npoints;++p)
      dwrite[p] = dmean[p]-derror[p];
    write_density(ofprefix+"_mean-err",dwrite);
  }

  template<typename VEC>
  void write_density(const std::string& outfile,VEC& density)
  {
    if(format=="xsf")
      write_density_xsf(outfile,density);
  }

  template<typename VEC>
  void write_density_xsf(const std::string& outfile,VEC& density)
  {
    using Units::convert;
    using Units::B;
    using Units::A;

    ParticleLayout_t& cell = input.cell;
    ParticleSet& Pq = *input.Peln;
    ParticleSet& Pc = *input.Pion;

    std::ofstream file;
    std::string filename = outfile+".xsf";
    file.open(filename.c_str(),std::ios::out | std::ios::trunc);
    if(!file.is_open())
      APP_ABORT("SpinDensityAnalzer::write_density_xsf\n  failed to open file for output: "+outfile+".xsf");

    file.precision(6);
    file<<std::scientific;
    int columns = 5;

    int natoms = Pc.getTotalNum();

    file<<" CRYSTAL"<< std::endl;
    file<<" PRIMVEC"<< std::endl;
    for(int i=0;i<DIM;++i)
    {
      file<<" ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Pq.Lattice.Rv[i][d],B,A);
      file<< std::endl;
    }
    file<<" PRIMCOORD"<< std::endl;
    file<<"   "<<natoms<<"   1"<< std::endl;
    for(int i=0;i<natoms;++i)
    {
      file<<"   "<<Pc.mySpecies.speciesName[Pc.GroupID[i]];
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(Pc.R[i][d],B,A);
      file<< std::endl;
    }
    file<<" BEGIN_BLOCK_DATAGRID_3D"<< std::endl;
    file<<"   "<<outfile<< std::endl;
    file<<"   DATAGRID_3D_SPIN_DENSITY"<< std::endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<grid[d];
    file<< std::endl;
    file<<"   ";
    for(int d=0;d<DIM;++d)
      file<<"  "<<convert(corner[d],B,A);
    file<< std::endl;
    for(int i=0;i<DIM;++i)
    {
      file<<"   ";
      for(int d=0;d<DIM;++d)
        file<<"  "<<convert(cell.Rv[i][d],B,A);
      file<< std::endl;
    }
    file<<"   ";
    for(int p=0;p<npoints;++p)
    {
      file<<"  "<<density[p];
      if((p+1)%columns==0)
        file<< std::endl<<"   ";
    }
    file<< std::endl;
    file<<"   END_DATAGRID_3D"<< std::endl;
    file<<" END_BLOCK_DATAGRID_3D"<< std::endl;
  }

};



struct TraceAnalyzer
{
  Communicate* comm;
  bool master;
  int master_rank;
  std::vector<int> local_procs;
  Input input;
  std::vector<QuantityAnalyzer*> quantities;
  size_t nblocks_min;
  Timer time;

  TraceAnalyzer(Communicate* c)
    : comm(c),input(c),master_rank(0)
  {
  }

  void analyze(char* infile)
  {
    if (comm->rank() != 0) {
      outputManager.shutOff();
    }

    Random.init(0,1,11);
    master = comm->rank()==master_rank;

    //read input file
    input.read(infile);
    input.report();
    QuantityTraceBase::verbose = input.verbose;

    //determine partitioning of files across mpi tasks
    partition_trace_procs();

    //process trace files
    app_log()<<"\nProcessing trace files by series and group"<< std::endl;
    for(int s=0;s<input.series.size();++s)
    {
      // setup quantity analyzers
      initialize_quantities(s);

      for(int g=0;g<input.groups.size();++g)
      {
        // accumulate trace data into each quantity (block form)
        accumulate_traces(s,g);

        // average over blocks and write data per group
        if(!input.average_groups)
          collect_data_and_write(s);
      }
      // average over blocks and write group average data
      if(input.average_groups)
        collect_data_and_write(s);

      // cleanup quantity analyzers
      finalize_quantities();
    }
  }

  void partition_trace_procs()
  {
    app_log()<<"\nPartitioning of trace procs across mpi ranks"<< std::endl;
    int nmpi = comm->size();
    for(int rank=0;rank<nmpi;++rank)
    {
      int nprocs = input.procs.size()/nmpi;
      int nrem   = input.procs.size()%nmpi;
      int pstart;
      if(rank<nrem)
      {
        nprocs++;
        pstart = rank*nprocs;
      }
      else
      {
        pstart = nrem*(nprocs+1) + (rank-nrem)*nprocs;
      }
      if(rank==comm->rank())
      {
        local_procs.resize(nprocs);
        for(int p=0;p<nprocs;++p)
          local_procs[p] = input.procs[pstart+p];
      }

      app_log()<<"  rank "<<rank<<"  procs ";
      for(int p=0;p<nprocs;++p)
        app_log()<<input.procs[pstart+p]<<" ";
      app_log()<< std::endl;
    }
  }

  void initialize_quantities(int s)
  {
    if(input.average_groups)
      nblocks_min = 2000000000;
    app_log()<<"  initializing quantities for series "<<input.series[s]<< std::endl;
    xmlNodePtr cur = input.quantities->children;
    while(cur!=NULL)
    {
      std::string cname((const char*)cur->name);
      QuantityAnalyzer* qa = 0;
      if(cname=="spindensity")
        qa = new SpinDensityAnalyzer(cname,input);
      else if(cname!="text" && cname!="comment")
        APP_ABORT(cname+" is not a valid quantity to analyze");
      if(qa)
      {
        qa->init(cur,s);
        app_log()<<"    found "<<qa->name<<"  mbytes = "<<(qa->memory()>>20)<< std::endl;
        quantities.push_back(qa);
      }
      cur=cur->next;
    }
  }

  void accumulate_traces(int s,int g)
  {
    if(!input.average_groups)
      nblocks_min = 2000000000;
    //accumulate trace data from each file
    app_log()<<"  processing trace files for ";
    if(input.groups[g]==-1)
      app_log()<<"series "<<input.series[s]<< std::endl;
    else
      app_log()<<"series "<<input.series[s]<<" group "<<input.groups[g]<< std::endl;
    for(int q=0;q<quantities.size();++q)
    {
      QuantityAnalyzer& qa = *quantities[q];
      qa.set_group(g);
    }
    real_t read_time = 0.0;
    real_t accum_time = 0.0;
    for(int ip=0;ip<local_procs.size();++ip)
    {
      int p = local_procs[ip];

      // get the trace file name
      std::string h5fname;
      get_trace_file_name(h5fname,s,g,p);

      // open the hdf archive
      if(input.verbose)
        app_log()<<"\n  reading "+h5fname<< std::endl;
      hdf_archive hin(comm);
      hin.open(h5fname,H5F_ACC_RDONLY);
      
      // read shared trace quantities
      if(input.verbose)
        app_log()<<"  reading shared quantities"<< std::endl;
      QuantityTrace<int>    steps;
      QuantityTrace<real_t> weights;
      steps.read(hin,"scalars/step");
      weights.read(hin,"scalars/weight");

      // determine block intervals and block weights
      std::vector<block_props> block_properties;
      determine_block_properties(s,steps,weights,block_properties);

      // accumulate trace data for each quantity
      if(input.verbose)
        app_log()<<"  accumulating quantities"<< std::endl;
      for(int q=0;q<quantities.size();++q)
      {
        QuantityAnalyzer& qa = *quantities[q];
        if(input.verbose)
          app_log()<<"   analyzing "<<qa.name<< std::endl;
        time.restart();
        qa.load_traces(hin);
        read_time += time.elapsed();

        time.restart();
        qa.accumulate(weights,block_properties);
        accum_time += time.elapsed();

        qa.clear_traces();
      }
    }
    app_log()<<"    processed "<<input.procs.size()<<" trace files"<< std::endl;
    app_log()<<"    read    time: "<<read_time<< std::endl;
    app_log()<<"    process time: "<<accum_time<< std::endl;
  }

  void get_trace_file_name( std::string& h5fname,int s,int g,int p)
  {
    int series = input.series[s];
    int group  = input.groups[g];
    int proc   = input.procs[p];
    std::string stoken;
    std::string gtoken;
    std::string ptoken;
    zfill(series,stoken,3);
    zfill(proc,ptoken,input.proc_digits);
    h5fname = input.path+"/"+input.prefix;
    if(group!=-1)
    {
      zfill(group,gtoken,3);
      h5fname+=".g"+gtoken;
    }
    if(input.tot_procs>1)
      h5fname+=".s"+stoken+".p"+ptoken+".traces.h5";
    else
      h5fname+=".s"+stoken+".traces.h5";
  }


  void determine_block_properties(int sindex,QuantityTrace<int>& steps,QuantityTrace<real_t>& weights,std::vector<block_props>& block_properties)
  {
    int steps_per_block = input.steps[sindex];
    int stepmax = steps_per_block;
    if(steps.trace.size()>0)
    {
      int i1 = 0;
      int i2 = i1;
      int step;
      for(int i=0;i<steps.trace.size();++i)
      {
        step = steps.trace[i];
        if(step==stepmax)
        {
          i2 = i;
          block_properties.push_back(block_props(i1,i2));
          i1 = i2;
          stepmax+=steps_per_block;
        }
      }
      if(step==stepmax-1)
      {
        i2 = steps.trace.size();
        block_properties.push_back(block_props(i1,i2));
      }
    }
    nblocks_min = std::min(nblocks_min,block_properties.size());
    if(input.verbose)
    {
      app_log()<<"    steps per block = "<<steps_per_block<< std::endl;
      app_log()<<"    nblocks         = "<<block_properties.size()<< std::endl;
    }
    for(int b=0;b<block_properties.size();++b)
    {
      block_props& bp = block_properties[b];
      bp.w = 0.0;
      for(int i=bp.i1;i<bp.i2;++i)
        bp.w += weights.trace[i];
    }
    //for(int b=0;b<block_properties.size();++b)
    //{
    //  block_props& bp = block_properties[b];
    //  app_log()<<"  block "<<b<<" "<<bp.i1<<" "<<bp.i2<<" "<<bp.w<< std::endl;
    //}
    return;
  }


  void collect_data_and_write(int s)
  {
#if USE_MPI
    //collect quantity data across mpi
    if(comm->size()>1)
    {
      app_log()<<"  collecting blocks across ranks "<< std::endl;
      time.restart();
      // compute mpi buffer offsets for each quantity
      int bsize = 1;
      for(int q=0;q<quantities.size();++q)
        quantities[q]->get_pack_start(bsize);
      std::vector<real_t> buffer(bsize);
      // pack, reduce, unpack each data block
      for(int b=0;b<input.blocks[s];++b)
      {
        buffer[0] = quantities[0]->block_weights[b];
        for(int q=0;q<quantities.size();++q)
          quantities[q]->pack(b,buffer);
        comm->reduce(buffer);
        if(master)
          for(int q=0;q<quantities.size();++q)
            quantities[q]->unpack(b,buffer);
      }
      // find the minimum number of blocks across ranks
      std::vector<int> nbm_in(1,nblocks_min);
      std::vector<int> nbm_out(1,nblocks_min);
      MPI_Reduce(&nbm_in[0],&nbm_out[0],1,MPI_INT,MPI_MIN,master_rank,comm->getMPI());
      if(master)
        nblocks_min = nbm_out[0];
      app_log()<<"    comm    time: "<<time.elapsed()<< std::endl;
    }
#endif

    //analyze stats and write quantity data
    if(master)  
    {
      app_log()<<"  writing quantities ("<<nblocks_min<<" out of "<<input.blocks[s]<<" blocks)"<< std::endl;
      time.restart();
      for(int q=0;q<quantities.size();++q)
      {
        QuantityAnalyzer& qa = *quantities[q];
        app_log()<<"    writing "<<qa.name<< std::endl;
        qa.analyze();
        qa.write();
      }
      app_log()<<"    write  time: "<<time.elapsed()<< std::endl;
    }
    comm->barrier();
  }

  void finalize_quantities()
  {
    //delete quantities
    delete_iter(quantities.begin(),quantities.end());
    quantities.clear();
  }
};



int main(int argc, char** argv)
{
  Communicate* comm = OHMMS::Controller;

  comm->initialize(argc,argv);

  TraceAnalyzer ta(comm);

  ta.analyze(argv[1]);

  comm->finalize();

  return 0;
}
