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
#include "Configuration.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/RandomGeneratorIO.h"
#include "Utilities/Timer.h"
#include "HDFVersion.h"
#include "OhmmsData/HDFAttribIO.h"
#include <fstream>

namespace APPNAMESPACE
{

template<>
struct HDFAttribIO<std::vector<uint32_t> >: public HDFAttribIOBase {

  typedef std::vector<uint32_t> arraytype_t;
  std::vector<hsize_t> dimsL,dimsG,offset;
  arraytype_t&  ref;
  bool replace;

  HDFAttribIO<arraytype_t>(arraytype_t& a, bool over=false): ref(a),replace(over)
  { 
    dimsL.resize(1,a.size());
    dimsG.resize(1,a.size());
    offset.resize(1,0);
  }

  template<unsigned D>
    HDFAttribIO<arraytype_t>(arraytype_t& a, TinyVector<int,D>& gcount, 
        TinyVector<int,D>& count, TinyVector<int,D>& start, bool overwrite=false): 
    ref(a), replace(overwrite)
    { 
      dimsG.resize(D);
      dimsL.resize(D);
      offset.resize(D);
      for(int i=0; i<D; i++) dimsG[i]=static_cast<hsize_t>(gcount[i]);
      for(int i=0; i<D; i++) dimsL[i]=static_cast<hsize_t>(count[i]);
      for(int i=0; i<D; i++) offset[i]=static_cast<hsize_t>(start[i]);
    }

  inline void write(hid_t grp, const char* name) {

    int shape=dimsG.size();
    hid_t dset_id;

    if(replace)
    {
      dset_id=H5Dopen(grp,name);
    }
    else
    {
      hid_t sid1  = H5Screate_simple(shape,&dimsG[0],NULL);
      dset_id=H5Dcreate(grp,name,H5T_NATIVE_UINT,sid1,H5P_DEFAULT);
      H5Sclose(sid1);
    }

    hsize_t stride[]={1,1,1,1};
    hid_t memspace=H5Screate_simple(shape,&dimsL[0],NULL);
    hid_t filespace=H5Dget_space(dset_id);
    herr_t ret=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,&offset[0],stride,&dimsL[0],NULL);
    ret = H5Dwrite(dset_id,H5T_NATIVE_UINT,memspace,filespace,xfer_plist,&ref[0]);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Dclose(dset_id);
  }

  inline void read(hid_t grp, const char* name) {
    hid_t dataset = H5Dopen(grp,name);
    hid_t dataspace = H5Dget_space(dataset);

    vector<hsize_t> gcount(dimsG.size());
    int rank_in = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, &gcount[0], NULL);

    //check the dimensions, accept only if the dimensions match
    if(rank_in == dimsG.size() && gcount[0] == dimsG[0]) 
    {
      for(int i=1;i<dimsG.size(); i++) dimsG[i]=gcount[i];
      hsize_t mreq=dimsL[0];
      for(int i=1;i<dimsG.size(); i++) {mreq *= (dimsL[i]=gcount[i]);}

      ref.resize(mreq);//resize it
      hid_t memspace = H5Screate_simple(rank_in, &dimsL[0], NULL);
      herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, &offset[0],NULL,&dimsL[0],NULL);
      status = H5Dread(dataset, H5T_NATIVE_UINT, memspace, dataspace, xfer_plist, &(ref[0]));
      H5Sclose(memspace);
    }
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
};

  /// constructors and destructors
  RandomNumberControl::RandomNumberControl(const char* aname)
    :OhmmsElementBase(aname), NeverBeenInitialized(true), myCur(NULL), Offset(5)
  { }
    
  /// generic output
  bool RandomNumberControl::get(std::ostream& os) const {
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
    int iseed=static_cast<int>(std::time(0))%4096;
    OHMMS::Controller->bcast(iseed);//broadcast the seed

    Offset=iseed;
    vector<uint_type> mySeeds;
    PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);
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

    ostringstream o;
    o << "  Random seeds Node = " << rank << ":";
    for(int ip=0; ip<nthreads; ip++)
    {
      int offset=baseoffset+ip;
      Children[ip]->init(rank,nprocs,myprimes[ip],offset);
      o << setw(12) << myprimes[ip];
    }
    cout << o.str() << endl;
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
      int offset_in = -1; // default is to generate by Wall-clock
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

      if(offset_in<0)
      {
        offset_in=static_cast<int>(std::time(0))%4096;
        OHMMS::Controller->bcast(offset_in);//broadcast the seed
      }
      else
        offset_in%=4096;

      Offset=offset_in;
      vector<uint_type> mySeeds;
      //allocate twice of what is required
      PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);

      Random.init(pid,nprocs,mySeeds[pid],Offset+pid);

      app_log() << "  Random number offset = " << Offset 
        <<  "  seeds = " << mySeeds[0] <<"-" << mySeeds[nprocs*omp_get_max_threads()] <<  endl;

      int imax=8*(mySeeds.size()/8);
      int jmax=std::min(std::size_t(8),mySeeds.size());
      for(int i=0; i<imax;)
      {
        for(int j=0; j<jmax; j++, i++) app_log() <<  std::setw(12) << mySeeds[i];
        app_log() << endl;
      }
      for(int i=imax; i<mySeeds.size(); i++) app_log() <<  std::setw(12) << mySeeds[i];
      app_log() << endl;

      make_children();
      NeverBeenInitialized = false; 
    }
    else
      reset();
    return true;
  }

  void RandomNumberControl::read(const string& fname, Communicate* comm)
  {
    string h5name(fname);
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");

    hid_t h_plist=H5P_DEFAULT;
    hid_t xfer_plist=H5P_DEFAULT;

#if defined(H5_HAVE_PARALLEL)
    if(comm->size()>1)
    {
      MPI_Info info=MPI_INFO_NULL;
      h_plist = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(h_plist,comm->getMPI(),info);
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
    }
#endif

    hid_t fid =  H5Fopen(h5name.c_str(),H5F_ACC_RDONLY,h_plist);
    herr_t status = H5Eset_auto(NULL, NULL);
    HDFVersion res_version(0,4);
    HDFVersion in_version(0,1);
    in_version.read(fid,hdf::version);
    if(in_version >= res_version)
    {
      hid_t h1 =  H5Gopen(fid,hdf::main_state);
      hid_t h2 = H5Gopen(h1,"random");

      int nthreads=omp_get_max_threads();
      TinyVector<int,2> gDims(comm->size()*nthreads,1);
      TinyVector<int,2> Dims(nthreads,1);
      TinyVector<int,2> offsets(comm->rank()*nthreads,0);

      vector<uint_type> vt;
      HDFAttribIO<std::vector<uint_type> > hin(vt,gDims,Dims,offsets);
      hin.setTransferProperty(xfer_plist);
      hin.read(h2,Random.EngineName.c_str());
      if(vt.size())
      {
        char fname[128];
        sprintf(fname,"random.p%i",comm->rank());
        ofstream fout(fname);
        std::copy(vt.begin(),vt.end(),ostream_iterator<uint_type>(fout,"\n"));

        std::stringstream otemp;
        std::copy(vt.begin(),vt.end(),ostream_iterator<uint_type>(otemp," "));
        if(nthreads>1)
        {
          for(int ip=0; ip<nthreads; ip++) Children[ip]->read(otemp);
        }
        else
          Random.read(otemp);
      }
      H5Gclose(h2);
      H5Gclose(h1);
    }
    else
    {
      app_warning() << "  Old configuration files. Cannot read random states.\n"
        << "  Using new random number seeds generated." << endl;
    }

    H5Fclose(fid);
  }

  void RandomNumberControl::write(const string& fname, Communicate* comm)
  {
    int pid=comm->rank();
    int nthreads=omp_get_max_threads();
    std::stringstream otemp;
    if(nthreads>1)
      for(int ip=0; ip<nthreads; ip++) Children[ip]->write(otemp);
    else
      Random.write(otemp);
    vector<uint_type> vt;
    std::copy(istream_iterator<uint_type>(otemp), istream_iterator<uint_type>(),back_inserter(vt));

    string h5name(fname);
    //append .config.h5 if missing
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");
    hid_t h_plist=H5P_DEFAULT;
    hid_t xfer_plist=H5P_DEFAULT;
#if defined(H5_HAVE_PARALLEL)
    if(comm->size()>1)
    {
      MPI_Info info=MPI_INFO_NULL;
      h_plist = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(h_plist,comm->getMPI(),info);
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
    }
#endif
    hid_t fid =  H5Fopen(h5name.c_str(),H5F_ACC_RDWR,h_plist);
    hid_t h1 =  H5Gopen(fid,hdf::main_state);
    hid_t h2 = H5Gcreate(h1,"random",0);

    TinyVector<int,2> gDims(comm->size()*nthreads,vt.size()/nthreads);
    TinyVector<int,2> Dims(nthreads,vt.size()/nthreads);
    TinyVector<int,2> offsets(pid*nthreads,0);

    HDFAttribIO<std::vector<uint_type> > hout(vt,gDims,Dims,offsets);
    hout.setTransferProperty(xfer_plist);
    hout.write(h2,Random.EngineName.c_str());

    //cleanup H5P
    if(xfer_plist != H5P_DEFAULT) H5Pclose(xfer_plist);
    if(h_plist != H5P_DEFAULT) H5Pclose(h_plist);
    H5Gclose(h2);
    H5Gclose(h1);
    H5Fclose(fid);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
