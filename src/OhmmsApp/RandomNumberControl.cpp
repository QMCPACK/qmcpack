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

    uint_type iseed=static_cast<uint_type>(std::time(0)%4096);
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
    if(nprocs<4)
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
    const int nthreads=omp_get_max_threads();
    vector<uint_type> vt,vt_tot;

    {//check the size
      std::stringstream otemp;
      if(nthreads>1)
        for(int ip=0; ip<nthreads; ip++) Children[ip]->write(otemp);
      else
        Random.write(otemp);
      std::copy(istream_iterator<uint_type>(otemp), istream_iterator<uint_type>(),back_inserter(vt));
    }

    vt_tot.resize(vt.size()*comm->size());
    TinyVector<hsize_t,2> shape(0);

    {//read it
      string h5name(fname);
      if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");
      hdf_archive hout(comm);
      hout.open(h5name,H5F_ACC_RDWR);
      hout.push(hdf::main_state);
      hout.push("random");

      hyperslab_proxy<vector<uint_type>,2> slab(vt_tot,shape);
      hout.read(slab,Random.EngineName);
      shape[0]=slab.size(0);
      shape[1]=slab.size(1);
      //HDFAttribIO<PooledData<uint_type> > o(vt_tot,shape);
      //o.read(hout.top(),Random.EngineName);
      //shape[0]=o.size(0);
      //shape[1]=o.size(1);
    }
    
    mpi::bcast(*comm,shape);
    if(shape[0]!=comm->size()*nthreads)
    {
      app_error() << "The number of parallel threads has changed from " << shape[0]
        << " to " << comm->size()*nthreads << endl;
      return;
    }
    mpi::scatter(*comm,vt_tot,vt);

    {
      std::stringstream otemp;
      std::copy(vt.begin(),vt.end(),ostream_iterator<uint_type>(otemp," "));
      if(nthreads>1)
        for(int ip=0; ip<nthreads; ip++) Children[ip]->read(otemp);
      else
        Random.read(otemp);
    }
  }

  void RandomNumberControl::write(const string& fname, Communicate* comm)
  {
    const int nthreads=omp_get_max_threads();
    std::stringstream otemp;
    if(nthreads>1)
      for(int ip=0; ip<nthreads; ip++) Children[ip]->write(otemp);
    else
      Random.write(otemp);

    vector<uint_type> vt,vt_tot;
    std::copy(istream_iterator<uint_type>(otemp), istream_iterator<uint_type>(),back_inserter(vt));
    vt_tot.resize(vt.size()*comm->size());
    mpi::gather(*comm,vt,vt_tot);

    //append .config.h5 if missing
    string h5name(fname);
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");
    hdf_archive hout(comm);
    hout.open(h5name,H5F_ACC_RDWR);
    hout.push(hdf::main_state);
    hout.push("random");
    TinyVector<hsize_t,2> shape(comm->size()*nthreads,vt.size()/nthreads);
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
