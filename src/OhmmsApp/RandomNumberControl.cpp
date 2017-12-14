//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include <Configuration.h>
#include <Message/OpenMP.h>
#include <OhmmsData/AttributeSet.h>
#include <OhmmsApp/RandomNumberControl.h>
#include <Utilities/RandomGeneratorIO.h>
#include <Utilities/Timer.h>
#include <HDFVersion.h>
#include <io/hdf_archive.h>
#include <mpi/collectives.h>
#if defined(HAVE_LIBBOOST)
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>
#endif
#include <Utilities/SimpleParser.h>

namespace qmcplusplus
{
///initialize the static data members
PrimeNumberSet<RandomGenerator_t::uint_type> RandomNumberControl::PrimeNumbers;
std::vector<RandomGenerator_t*>  RandomNumberControl::Children;
RandomGenerator_t::uint_type RandomNumberControl::Offset=11u;

/// constructors and destructors
RandomNumberControl::RandomNumberControl(const char* aname)
  :OhmmsElementBase(aname), NeverBeenInitialized(true), myCur(NULL)//, Offset(5)
{ }

/// generic output
bool RandomNumberControl::get(std::ostream& os) const
{
  if(omp_get_max_threads()>1)
  {
    for(int ip=0; ip<omp_get_max_threads(); ip++)
    {
      Children[ip]->write(os);
      os << std::endl;
    }
  }
  else
  {
    Random.write(os);
  }
  return true;
}

/// generic input
bool RandomNumberControl::put(std::istream& is)
{
  return true;
}

/// reset the generator
void RandomNumberControl::reset()
{
  make_seeds();
}

/// reset the generator
void RandomNumberControl::make_seeds()
{
  int pid = OHMMS::Controller->rank();
  int nprocs = OHMMS::Controller->size();
  uint_type iseed=static_cast<uint_type>(std::time(0))%1024;
  mpi::bcast(*OHMMS::Controller,iseed);
  //OHMMS::Controller->bcast(iseed);//broadcast the seed
  Offset=iseed;
  std::vector<uint_type> mySeeds;
  RandomNumberControl::PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);
  Random.init(pid,nprocs,mySeeds[pid],Offset+pid);
  //change children as well
  make_children();
}

void RandomNumberControl::make_children()
{
  int nthreads=omp_get_max_threads();
  int n=nthreads-Children.size();
  while(n)
  {
    Children.push_back(new RandomGenerator_t);
    n--;
  }
  int rank=OHMMS::Controller->rank();
  int nprocs=OHMMS::Controller->size();
  int baseoffset=Offset+nprocs+nthreads*rank;
  std::vector<uint_type> myprimes;
  PrimeNumbers.get(baseoffset,nthreads,myprimes);
  for(int ip=0; ip<nthreads; ip++)
  {
    int offset=baseoffset+ip;
    Children[ip]->init(rank,nprocs,myprimes[ip],offset);
  }
}

xmlNodePtr
RandomNumberControl::initialize(xmlXPathContextPtr acontext)
{
  xmlXPathObjectPtr rg_request
    = xmlXPathEvalExpression((const xmlChar*)"//random",acontext);
  if(xmlXPathNodeSetIsEmpty(rg_request->nodesetval))
    put(NULL);
  else
    put(rg_request->nodesetval->nodeTab[0]);
  xmlXPathFreeObject(rg_request);
  return myCur;
}

void RandomNumberControl::test()
{
  /* Add random number generator tester
  */
  int nthreads=omp_get_max_threads();
  std::vector<double> avg(nthreads),avg2(nthreads);
  #pragma omp parallel for
  for(int ip=0; ip<nthreads; ++ip)
  {
    const int n=1000000;
    double sum=0.0, sum2=0.0;
    RandomGenerator_t& myrand(*Children[ip]);
    for(int i=0; i<n; ++i)
    {
      double r=myrand.rand();
      sum +=r;
      sum2+= r*r;
    }
    avg[ip]=sum/static_cast<double>(n);
    avg2[ip]=sum2/static_cast<double>(n);
  }
  std::vector<double> avg_tot(nthreads*OHMMS::Controller->size()),avg2_tot(nthreads*OHMMS::Controller->size());
  mpi::gather(*OHMMS::Controller,avg,avg_tot);
  mpi::gather(*OHMMS::Controller,avg2,avg2_tot);
  double avg_g=0.0;
  double avg2_g=0.0;
  for(int i=0,ii=0; i<OHMMS::Controller->size(); ++i)
  {
    for(int ip=0; ip<nthreads; ++ip,++ii)
    {
      app_log() << "RNGTest " << std::setw(4) << i << std::setw(4) << ip
                << std::setw(20) << avg_tot[ii] << std::setw(20) << avg2_tot[ii]-avg_tot[ii]*avg_tot[ii] << std::endl;
      avg_g+=avg_tot[ii];
      avg2_g+=avg2_tot[ii];
    }
  }
  avg_g/=static_cast<double>(nthreads*OHMMS::Controller->size());
  avg2_g/=static_cast<double>(nthreads*OHMMS::Controller->size());
  app_log() << "RNGTest " << std::setw(4) << OHMMS::Controller->size() << std::setw(4) << nthreads
            << std::setw(20) << avg_g << std::setw(20) << avg2_g-avg_g*avg_g<< std::endl;
  app_log().flush();
}

bool RandomNumberControl::put(xmlNodePtr cur)
{
  if(NeverBeenInitialized)
  {
    bool init_mpi = true;
    int offset_in = -1; // default is to generate by Wall-clock
    if(cur != NULL)
    {
      std::string pname("yes");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(pname,"parallel");
      oAttrib.add(offset_in,"seed");
      oAttrib.put(cur);
      if(pname == "0" || pname == "false" || pname == "no")
        init_mpi=false;
    }
    int nprocs = 1;
    int pid = 0;
    if(init_mpi)
    {
      pid = OHMMS::Controller->rank();
      nprocs = OHMMS::Controller->size();
    }
    app_summary() << " Random Number" << std::endl;
    app_summary() << " -------------" << std::endl;
    if(offset_in<0)
    {
      offset_in=static_cast<int>(static_cast<uint_type>(std::time(0))%1024);
      app_summary() << "  Offset for the random number seeds based on time: " << offset_in << std::endl;
      mpi::bcast(*OHMMS::Controller,offset_in);
    }
    else
    {
      offset_in%=1024;
      app_summary() << "  Offset for the random number seeds from input file (mod 1024): " << offset_in << std::endl;
    }
    app_summary() << std::endl;
    Offset=offset_in;
    std::vector<uint_type> mySeeds;
    //allocate twice of what is required
    PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);
    Random.init(pid,nprocs,mySeeds[pid],Offset+pid);
    app_log() <<  "  Range of prime numbers to use as seeds over processors and threads = " << mySeeds[0] <<"-" << mySeeds[nprocs*omp_get_max_threads()] <<  std::endl;
    app_log() << std::endl;

    make_children();
    NeverBeenInitialized = false;
    app_log() << std::endl;
  }
  else
    reset();
  return true;
}

void RandomNumberControl::read_old(const std::string& fname, Communicate* comm)
{
  int nthreads=omp_get_max_threads();
  std::vector<uint_type> vt_tot, vt;
  std::vector<int> shape(2,0),shape_now(2,0);
  shape_now[0]=comm->size()*nthreads;
  shape_now[1]=Random.state_size();

  if(comm->rank()==0)
  {
#if defined(HAVE_LIBBOOST)
    using boost::property_tree::ptree;
    ptree pt;
    std::string xname=fname+".random.xml";
    read_xml(xname, pt);
    if(!pt.empty())
    {
      std::string engname=pt.get<std::string>("random.engine");
      if(engname==Random.EngineName)
      {
        std::istringstream dims(pt.get<std::string>("random.dims"));
        dims >> shape[0] >> shape[1];
        if(shape[0]==shape_now[0] && shape[1]==shape_now[1])
        {
          vt_tot.resize(shape[0]*shape[1]);
          std::istringstream v(pt.get<std::string>("random.states"));
          for(int i=0; i<vt_tot.size(); ++i) v>>vt_tot[i];
        }
        else
          shape[0]=shape[1]=0;
      }
    }
#else
    TinyVector<hsize_t,2> shape_t(0);
    shape_t[1]=Random.state_size();
    hyperslab_proxy<std::vector<uint_type>,2> slab(vt_tot,shape_t);
    std::string h5name=fname+".random.h5";
    hdf_archive hout(comm);
    hout.open(h5name,H5F_ACC_RDONLY);
    hout.push(hdf::main_state);
    hout.push("random");
    std::string engname;
    hout.read(slab,Random.EngineName);
    shape[0]=static_cast<int>(slab.size(0));
    shape[1]=static_cast<int>(slab.size(1));
#endif
  }

  mpi::bcast(*comm,shape);

  if(shape[0]!=shape_now[0] || shape[1] != shape_now[1])
  {
    app_log() << "Mismatched random number generators."
      << "\n  Number of streams     : old=" << shape[0] << " new= " << comm->size()*nthreads
      << "\n  State size per stream : old=" << shape[1] << " new= " << Random.state_size()
      << "\n  Using the random streams generated at the initialization." << std::endl;
    return;
  }

  app_log() << "  Restart from the random number streams from the previous configuration." << std::endl;
  vt.resize(nthreads*Random.state_size());

  if(comm->size()>1)
    mpi::scatter(*comm,vt_tot,vt);
  else
    copy(vt_tot.begin(),vt_tot.end(),vt.begin());

  {
    if(nthreads>1)
    {
      std::vector<uint_type>::iterator vt_it(vt.begin());
      for(int ip=0; ip<nthreads; ip++, vt_it += shape[1])
      {
        std::vector<uint_type> c(vt_it,vt_it+shape[1]);
        Children[ip]->load(c);
      }
    }
    else
      Random.load(vt);
  }
}

void RandomNumberControl::write_old(const std::string& fname, Communicate* comm)
{
  int nthreads=omp_get_max_threads();
  std::vector<uint_type> vt, vt_tot;
  vt.reserve(nthreads*1024);
  if(nthreads>1)
    for(int ip=0; ip<nthreads; ++ip)
    {
      std::vector<uint_type> c;
      Children[ip]->save(c);
      vt.insert(vt.end(),c.begin(),c.end());
    }
  else
    Random.save(vt);
  if(comm->size()>1)
  {
    vt_tot.resize(vt.size()*comm->size());
    mpi::gather(*comm,vt,vt_tot);
  }
  else
    vt_tot=vt;

  if(comm->rank()==0)
  {
#if defined(HAVE_LIBBOOST)
    using boost::property_tree::ptree;
    ptree pt;
    std::ostringstream dims,vt_o;
    dims<<comm->size()*nthreads << " " << Random.state_size();
    std::vector<uint_type>::iterator v=vt_tot.begin();
    for(int i=0; i<comm->size()*nthreads; ++i)
    {
      copy(v,v+Random.state_size(),std::ostream_iterator<uint_type>(vt_o," "));
      vt_o<< std::endl;
      v+=Random.state_size();
    }
    pt.put("random.engine", Random.EngineName);
    pt.put("random.dims",dims.str());
    pt.put("random.states",vt_o.str());
    std::string xname=fname+".random.xml";
    write_xml(xname, pt);
#else
    std::string h5name=fname+".random.h5";
    hdf_archive hout(comm);
    hout.create(h5name);
    hout.push(hdf::main_state);
    hout.push("random");
    TinyVector<hsize_t,2> shape(comm->size()*nthreads,Random.state_size());
    hyperslab_proxy<std::vector<uint_type>,2> slab(vt_tot,shape);
    hout.write(slab,Random.EngineName);
    hout.close();
#endif
  }
}

/*New functions past this point*/
//switch between read functions
void RandomNumberControl::read(const std::string& fname, Communicate* comm)
{
  std::string h5name=fname+".random.h5";
  hdf_archive hin(comm, true); //attempt to read in parallel
  hin.open(h5name,H5F_ACC_RDONLY);
  if(hin.is_parallel())
    read_parallel(hin, comm);
  else
    read_rank_0(hin, comm);
}

//switch between write functions
void RandomNumberControl::write(const std::string& fname, Communicate* comm)
{
  std::string h5name=fname+".random.h5";
  hdf_archive hout(comm, true); //attempt to write in parallel
  hout.create(h5name);
  if(hout.is_parallel())
    write_parallel(hout, comm);
  else
    write_rank_0(hout, comm);
}

//Parallel read
void RandomNumberControl::read_parallel(hdf_archive& hin, Communicate* comm)
{
  int nthreads = omp_get_max_threads();
  std::vector<uint_type> vt, mt;
  TinyVector<int,3> shape_now(comm->size(), nthreads, Random.state_size()); //cur configuration
  TinyVector<int,3> shape_hdf5(3,0); //configuration when file was written

  //grab shape and Random.state_size() used to create hdf5 file
  hin.push(hdf::main_state);
  hin.read(shape_hdf5, "nprocs_nthreads_statesize");

  //if hdf5 file's shape and the current shape don't match, abort read
  if(shape_hdf5[0] != shape_now[0] || shape_hdf5[1] != shape_now[1] || shape_hdf5[2] != shape_now[2])
  {
    app_log() << "Mismatched random number generators."
      << "\n  Number of procs in streams : old=" << shape_hdf5[0] << " new= " << shape_now[0]
      << "\n  Number of threads in streams : old=" << shape_hdf5[1] << " new= " << shape_now[1]
      << "\n  State size per stream : old=" << shape_hdf5[2] << " new= " << shape_now[2]
      << "\n  Using the random streams generated at the initialization.\n";
    return;
  }
  app_log() << "  Restart from the random number streams from the previous configuration.\n";

  TinyVector<int,2> shape(comm->size()*nthreads, Random.state_size()); //global dims of children dataset
  vt.resize(nthreads*Random.state_size()); //buffer for children[ip]
  mt.resize(Random.state_size()); //buffer for single thread Random object of random nums

  TinyVector<int,2> counts(nthreads, Random.state_size()); //local dimensions of dataset
  TinyVector<int,2> offsets(comm->rank() * nthreads, 0); //offsets for each process to read in

  hin.push("random"); //group that holds children[ip] random nums
  hyperslab_proxy<std::vector<uint_type>,2> slab(vt, shape, counts, offsets);
  hin.read(slab,Random.EngineName);

  hin.pop();
  hin.push("random_master"); //group that holds Random_th random nums
  shape[0] = comm->size(); //reset shape, counts and offset for non-multiple threads
  counts[0] = 1;
  offsets[0] = comm->rank();
  hyperslab_proxy<std::vector<uint_type>,2> slab2(mt, shape, counts, offsets);
  hin.read(slab2,Random.EngineName);
  hin.close();

  std::vector<uint_type>::iterator vt_it(vt.begin());
  for(int ip=0; ip<nthreads; ip++, vt_it += shape[1])
  {
    std::vector<uint_type> c(vt_it,vt_it+shape[1]);
    Children[ip]->load(c); //load random nums back to program from buffer
  }
  Random.load(mt); //load random nums back to prog from buffer
}

//Parallel write
void RandomNumberControl::write_parallel(hdf_archive& hout, Communicate* comm)
{
  int nthreads=omp_get_max_threads();
  std::vector<uint_type> vt, mt;
  TinyVector<int,3> shape_hdf5(comm->size(), nthreads, Random.state_size()); //configuration at write time
  vt.reserve(nthreads*Random.state_size()); //buffer for random numbers from children[ip] of each thread
  mt.reserve(Random.state_size()); //buffer for random numbers from single Random object

  for(int ip=0; ip<nthreads; ++ip)
  {
    std::vector<uint_type> c;
    Children[ip]->save(c);
    vt.insert(vt.end(),c.begin(),c.end()); //get nums from each thread into buffer
  }
  Random.save(mt); //get nums for single random object (no threads)

  TinyVector<int,2> shape(comm->size()*nthreads,Random.state_size()); //global dimensions
  TinyVector<int,2> counts(nthreads, Random.state_size()); //local dimensions
  TinyVector<int,2> offsets(comm->rank() * nthreads, 0); //offset for the file write

  hout.push(hdf::main_state);
  hout.write(shape_hdf5, "nprocs_nthreads_statesize"); //save the shape of the data at write

  hout.push("random"); //group for children[ip]
  hyperslab_proxy<std::vector<uint_type>,2> slab(vt, shape, counts, offsets);
  hout.write(slab,Random.EngineName); //write to hdf5file
  hout.pop();

  shape[0] = comm->size(); //adjust shape, counts, offset for just one thread
  counts[0] = 1;
  offsets[0] = comm->rank();
  hout.push("random_master"); //group for random object without threads
  hyperslab_proxy<std::vector<uint_type>,2> slab2(mt, shape, counts, offsets);
  hout.write(slab2,Random.EngineName); //write data to hdf5 file
  hout.close();
}

//Scatter read
void RandomNumberControl::read_rank_0(hdf_archive& hin, Communicate* comm)
{
  int nthreads = omp_get_max_threads();
  std::vector<uint_type> vt, vt_tot, mt, mt_tot;
  TinyVector<int,3> shape_now(comm->size(), nthreads, Random.state_size()); //current configuration
  TinyVector<int,2> shape(comm->size()*nthreads, Random.state_size()); //dimensions of children dataset
  TinyVector<int,3> shape_hdf5(3,0); //configuration when hdf5 file was written

  //grab configuration of threads/procs and Random.state_size() in hdf5 file
  if(comm->rank() == 0)
  {
    hin.push(hdf::main_state);
    hin.read(shape_hdf5, "nprocs_nthreads_statesize");
  }

  mpi::bcast(*comm, shape_hdf5);

  //if hdf5 file's configuration and current configuration don't match, abort read
  if(shape_hdf5[0] != shape_now[0] || shape_hdf5[1] != shape_now[1] || shape_hdf5[2] != shape_now[2])
  {
    app_log() << "Mismatched random number generators."
      << "\n  Number of procs in streams : old=" << shape_hdf5[0] << " new= " << shape_now[0]
      << "\n  Number of threads in streams : old=" << shape_hdf5[1] << " new= " << shape_now[1]
      << "\n  State size per stream : old=" << shape_hdf5[2] << " new= " << shape_now[2]
      << "\n  Using the random streams generated at the initialization.\n";
    return;
  }
  app_log() << "  Restart from the random number streams from the previous configuration.\n";

  vt.resize(nthreads*Random.state_size()); //buffer for random nums in children of each thread
  mt.resize(Random.state_size()); //buffer for random numbers from single Random object

  if(comm->rank() == 0)
  {
    hin.push("random"); //group for children[ip] (Random.object for each thread)
    vt_tot.resize(nthreads*Random.state_size()*comm->size());
    hyperslab_proxy<std::vector<uint_type>,2> slab(vt_tot, shape);
    hin.read(slab,Random.EngineName);
    hin.pop();

    shape[0] = comm->size(); //reset shape to one thread per process
    mt_tot.resize(Random.state_size()*comm->size());
    hin.push("random_master"); //group for single Random object
    hyperslab_proxy<std::vector<uint_type>,2> slab2(mt_tot, shape);
    hin.read(slab2,Random.EngineName);
    hin.close();
  }

  if(comm->size()>1)
  {
    mpi::scatter(*comm,vt_tot,vt); //divide big buffer into on for each proc
    mpi::scatter(*comm,mt_tot,mt);
  }
  else
  {
    copy(vt_tot.begin(),vt_tot.end(),vt.begin());
    copy(mt_tot.begin(),mt_tot.end(),mt.begin());
  }

  std::vector<uint_type>::iterator vt_it(vt.begin());
  for(int i=0; i<nthreads; i++, vt_it += shape[1])
  {
    std::vector<uint_type> c(vt_it,vt_it+shape[1]);
    Children[i]->load(c); //read seeds for each thread from buffer back into object
  }
  Random.load(mt); //read seeds back into object
}

//scatter write
void RandomNumberControl::write_rank_0(hdf_archive& hout, Communicate* comm)
{
  int nthreads = omp_get_max_threads();
  std::vector<uint_type> vt, vt_tot, mt, mt_tot;
  TinyVector<int,2> shape(comm->size()*nthreads, Random.state_size()); //dimensions of children dataset
  TinyVector<int,3> shape_hdf5(comm->size(), nthreads, Random.state_size()); //configuration at write time
  vt.reserve(nthreads*Random.state_size()); //buffer for children[ip] (Random object of seeds for each thread)
  mt.reserve(Random.state_size());  //buffer for single Random object of seeds, one per proc regardless of thread num

  for(int i=0; i<nthreads; ++i)
  {
    std::vector<uint_type> c;
    Children[i]->save(c);
    vt.insert(vt.end(),c.begin(),c.end()); //copy children[nthreads] seeds to buffer
  }
  Random.save(mt); //copy random_th seeds to buffer

  if(comm->size()>1)
  {
    vt_tot.resize(vt.size()*comm->size());
    mt_tot.resize(mt.size()*comm->size());
    mpi::gather(*comm,vt,vt_tot); //gather into one big buffer for master write
    mpi::gather(*comm,mt,mt_tot);
  }
  else
  {
    vt_tot=vt;
    mt_tot=mt;
  }

  if(comm->rank()==0)
  {
    hout.push(hdf::main_state);
    hout.write(shape_hdf5, "nprocs_nthreads_statesize"); //configuration at write time to file

    hout.push("random"); //group for children[ip]
    hyperslab_proxy<std::vector<uint_type>, 2> slab(vt_tot, shape);
    hout.write(slab, Random.EngineName);
    hout.pop();

    shape[0] = comm->size(); //reset dims for single thread use
    hout.push("random_master"); //group for random_th object
    hyperslab_proxy<std::vector<uint_type>,2> slab2(mt_tot, shape);
    hout.write(slab2,Random.EngineName);
    hout.close();
  }
}
}
