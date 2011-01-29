//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <Configuration.h>
#include <Message/OpenMP.h>
#include <OhmmsData/AttributeSet.h>
#include <OhmmsApp/RandomNumberControl.h>
#include <Utilities/RandomGeneratorIO.h>
#include <Utilities/Timer.h>
#include <HDFVersion.h>
#include <io/hdf_archive.h>
#include <mpi/collectives.h>

namespace APPNAMESPACE
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
      for(int ip=0; ip<omp_get_max_threads(); ip++) {
        Children[ip]->write(os); os << endl;
      }
    }
    else
    {
      Random.write(os);
    }
    return true;
  }

  /// generic input
  bool RandomNumberControl::put(std::istream& is) {
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
    vector<uint_type> mySeeds;
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
      Children.push_back(new RandomGenerator_t); n--;
    }

    int rank=OHMMS::Controller->rank();
    int nprocs=OHMMS::Controller->size();
    int baseoffset=Offset+nprocs+nthreads*rank;
    vector<uint_type> myprimes;
    PrimeNumbers.get(baseoffset,nthreads,myprimes);

    for(int ip=0; ip<nthreads; ip++)
    {
      int offset=baseoffset+ip;
      Children[ip]->init(rank,nprocs,myprimes[ip],offset);
    }

    if(nprocs<32)
    {
      ostringstream o;
      o << "  Random seeds Node = " << rank << ":";
      for(int ip=0; ip<nthreads; ip++)
        o << setw(12) << myprimes[ip];
      cout << o.str() << endl;
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

  bool RandomNumberControl::put(xmlNodePtr cur)
  {
    if(NeverBeenInitialized) {
      bool init_mpi = true;
      uint_type offset_in = 0; // default is to generate by Wall-clock
      if(cur != NULL) 
      {
        std::string pname("yes");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(pname,"parallel");
        oAttrib.add(offset_in,"seed");
        oAttrib.put(cur);
        if(pname == "0" || pname == "false" || pname == "no") init_mpi=false;
      }

      int nprocs = 1;
      int pid = 0;
      if(init_mpi) 
      {
        pid = OHMMS::Controller->rank();
        nprocs = OHMMS::Controller->size();
      }

      if(offset_in==0)
      {
        offset_in=static_cast<uint_type>(std::time(0))%1024;
        mpi::bcast(*OHMMS::Controller,offset_in);
      }
      else
        offset_in%=1024;

      Offset=offset_in;
      vector<uint_type> mySeeds;
      //allocate twice of what is required
      PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);

      Random.init(pid,nprocs,mySeeds[pid],Offset+pid);

      app_log() << "  Random number offset = " << Offset 
        <<  "  seeds = " << mySeeds[0] <<"-" << mySeeds[nprocs*omp_get_max_threads()] <<  endl;

      if(nprocs<4)
      {
        int imax=8*(mySeeds.size()/8);
        int jmax=std::min(std::size_t(8),mySeeds.size());
        for(int i=0; i<imax;)
        {
          for(int j=0; j<jmax; j++, i++) app_log() <<  std::setw(12) << mySeeds[i];
          app_log() << endl;
        }
        for(int i=imax; i<mySeeds.size(); i++) app_log() <<  std::setw(12) << mySeeds[i];
        app_log() << endl;
      }

      make_children();
      NeverBeenInitialized = false; 
    }
    else
      reset();
    return true;
  }

  void RandomNumberControl::read(const string& fname, Communicate* comm)
  {
    int nthreads=omp_get_max_threads();

    vector<uint_type> vt_tot(comm->size()*nthreads*Random.state_size()), vt;

    TinyVector<int,2> shape(0);
    {//read it
      string h5name(fname);
      if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");
      hdf_archive hout(comm);
      hout.open(h5name,H5F_ACC_RDONLY);
      hout.push(hdf::main_state);
      hout.push("random");

      TinyVector<hsize_t,2> shape_t(0);
      hyperslab_proxy<vector<uint_type>,2> slab(vt_tot,shape_t);
      hout.read(slab,Random.EngineName);
      shape[0]=slab.size(0);
      shape[1]=slab.size(1);
    }

    //bcast of hsize_t is not working
    mpi::bcast(*comm,shape);

    if(shape[0]!=comm->size()*nthreads || shape[1] != Random.state_size())
    {
      app_warning() << "Mismatched random number generators."
        << "\n  Number of streams     : old=" << shape[0] << " new= " << comm->size()*nthreads 
        << "\n  State size per stream : old=" << shape[1] << " new= " << Random.state_size()
        << "\n  Using the random streams generated at the initialization." << endl;
      return;
    }

    app_log() << "  Restart from the random number streams from the previous configuration." << endl;

    vt.resize(nthreads*Random.state_size());
    if(comm->size()>1)
      mpi::scatter(*comm,vt_tot,vt);
    else
      std::copy(vt_tot.begin(),vt_tot.begin()+vt.size(),vt.begin());

    {
      if(nthreads>1)
      {
        vector<uint_type>::iterator vt_it(vt.begin());
        for(int ip=0; ip<nthreads; ip++, vt_it += shape[1]) 
        {
          vector<uint_type> c(vt_it,vt_it+shape[1]);
          Children[ip]->load(c);
        }
      }
      else
        Random.load(vt);
    }

    /* Add random number generator tester 
    const int n=1000000;
    double avg=0.0,avg2=0.0;
#pragma omp parallel reduction(+:avg,avg2)
    {
      double sum=0.0, sum2=0.0;
      int ip=omp_get_thread_num();
      RandomGenerator_t& myrand(*Children[ip]);
      for(int i=0; i<n; ++i)
      {
        double r=myrand.rand();
        sum +=r; sum2+= r*r;
      }
      avg+=sum; avg2+=sum2;
    }
    comm->allreduce(avg);
    comm->allreduce(avg2);
    avg/=static_cast<double>(nthreads*n);
    avg2/=static_cast<double>(nthreads*n);
    app_log() << " Average = " << avg << "  Variance = " << avg2-avg*avg << endl;
    */
  }

  void RandomNumberControl::write(const string& fname, Communicate* comm)
  {
    int nthreads=omp_get_max_threads();
    vector<uint_type> vt, vt_tot;
    vt.reserve(nthreads*1024);
    if(nthreads>1)
      for(int ip=0; ip<nthreads; ++ip)
      {
        vector<uint_type> c;
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

    string h5name(fname);
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");
    hdf_archive hout(comm);
    hout.open(h5name,H5F_ACC_RDWR);
    hout.push(hdf::main_state);
    hout.push("random");
    TinyVector<hsize_t,2> shape(comm->size()*nthreads,Random.state_size());
    hyperslab_proxy<vector<uint_type>,2> slab(vt_tot,shape);
    hout.write(slab,Random.EngineName);
    hout.close();
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
