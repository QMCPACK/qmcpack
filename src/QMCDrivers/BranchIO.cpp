//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/BranchIO.h"
#include "HDFVersion.h"
#include "Message/CommOperators.h"
#include "io/hdf_archive.h"
//#include <boost/archive/text_oarchive.hpp>
//
namespace qmcplusplus 
{
  template<typename T>
    struct h5data_proxy<accumulator_set<T> >: public h5_space_type<T,1>
    {
      enum {CAPACITY=accumulator_set<T>::CAPACITY};
      using h5_space_type<T,1>::dims;
      using h5_space_type<T,1>::get_address;
      typedef accumulator_set<T> data_type;
      data_type& ref_;

      inline h5data_proxy(data_type& a): ref_(a) { dims[0]=CAPACITY;}
      inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
      {
        return h5d_read(grp,aname,get_address(ref_.properties),xfer_plist);
      }

      inline bool write(hid_t grp, const std::string& aname,hid_t xfer_plist=H5P_DEFAULT)
      {
        return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.properties),xfer_plist);
      }

      /** return the start address */
      inline T* begin() { return ref_.properties;}
      /** return the end address */
      inline T* end() { return ref_.properties+CAPACITY;}
    };

  bool BranchIO::write(const string& fname) 
  {
    //append .config.h5 if missing
    string h5name(fname);
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");

    hdf_archive dump(myComm,false);
    hid_t fid = dump.open(h5name);
    //cannot find the file, return false
    if(fid<0) return false;
    dump.push(hdf::main_state);
    bool firsttime=!dump.is_group(hdf::qmc_status);
    dump.push(hdf::qmc_status);
    if(firsttime)
    {
      string v_header("tau:taueff:etrial:eref:branchmax:branchcutoff:branchfilter:sigma:acc_energy:acc_samples");
      string i_header("warmupsteps:energyupdateinterval:counter:targetwalkers:maxwalkers:minwalkers:branchinterval");
      dump.write(v_header,"vparam_def");
      dump.write(i_header,"iparam_def");
    }
    dump.write(ref.vParam,"vparam");
    dump.write(ref.iParam,"iparam");
    dump.write(ref.BranchMode,"branchmode");

    dump.push("histogram");
    dump.write(ref.EnergyHist,"energy");
    dump.write(ref.VarianceHist,"variance");
    dump.write(ref.R2Accepted,"r2accepted");
    dump.write(ref.R2Proposed,"r2proposed");

    //PopHist is not being used in 2010-10-19
    //if(ref.BranchMode[SimpleFixedNodeBranch::B_DMC])
    //{
    //  dump.push("population");
    //  dump.write(ref.PopHist.myData,"histogram");
    //}

    return true;
  }

  bool BranchIO::read(const string& fname) 
  {
    //append .config.h5 if missing
    string h5name(fname);
    if(fname.find("config.h5")>= fname.size()) h5name.append(".config.h5");

    //do not use collective
    hdf_archive prevconfig(myComm,false);
    if(!prevconfig.open(h5name,H5F_ACC_RDONLY)) return false;

    HDFVersion res_version(0,4); //start using major=0 and minor=4
    HDFVersion res_20080624(0,5);//major revision on 2008-06-24 0.5
    HDFVersion in_version(0,1);

    int n=ref.vParam.size()+ref.iParam.size();
    /** temporary storage to broadcast restart data */
    vector<RealType> pdata(n+3+16,-1);
    prevconfig.read(in_version.version,hdf::version);
    myComm->bcast(in_version.version);

    prevconfig.push(hdf::main_state,false);
    if(in_version>= res_20080624)
    {//need a better version control
      prevconfig.push(hdf::qmc_status,false);

      prevconfig.read(ref.vParam,"vparam");
      prevconfig.read(ref.iParam,"iparam");
      prevconfig.read(ref.BranchMode,"branchmode");

      std::copy(ref.vParam.begin(),ref.vParam.end(),pdata.begin());
      std::copy(ref.iParam.begin(),ref.iParam.end(),pdata.begin()+ref.vParam.size());

      int offset=n;
      pdata[offset++]=ref.BranchMode.to_ulong();
      pdata[offset++]=in_version[0];
      pdata[offset++]=in_version[1];

      {//get histogram
        prevconfig.push("histogram",false);

        accumulator_set<RealType> temp;

        prevconfig.read(temp,"energy"); 
        std::copy(temp.properties,temp.properties+4,pdata.begin()+offset);
        offset+=4;

        prevconfig.read(temp,"variance"); 
        std::copy(temp.properties,temp.properties+4,pdata.begin()+offset);
        offset+=4;

        prevconfig.read(temp,"r2accepted"); 
        std::copy(temp.properties,temp.properties+4,pdata.begin()+offset);
        offset+=4;

        prevconfig.read(temp,"r2proposed"); 
        std::copy(temp.properties,temp.properties+4,pdata.begin()+offset);
        offset+=4;
        prevconfig.pop();
      }
      prevconfig.pop();
    }
    prevconfig.pop();

    //broadcast to the nodes : need to add a namespace mpi::
    myComm->bcast(pdata);

    if(myComm->rank())
    {
      //\since 2008-06-24 
      if(in_version>=res_20080624)
      {
        int ii=0;
        for(int i=0; i<ref.vParam.size(); ++i,++ii) ref.vParam[i]=pdata[ii];
        for(int i=0; i<ref.iParam.size(); ++i,++ii) ref.iParam[i]=static_cast<int>(pdata[ii]);
        ref.BranchMode=static_cast<unsigned long>(pdata[ii]);
      }
    }

    {//update historgram
      int ii=n+3;
      ref.EnergyHist.reset(pdata[ii],pdata[ii+1],pdata[ii+2]); ii+=4;
      ref.VarianceHist.reset(pdata[ii],pdata[ii+1],pdata[ii+2]); ii+=4;
      ref.R2Accepted.reset(pdata[ii],pdata[ii+1],pdata[ii+2]); ii+=4;
      ref.R2Proposed.reset(pdata[ii],pdata[ii+1],pdata[ii+2]); ii+=4;
    }

    return true;
  }
}

/***************************************************************************
 * $RCSfile: BranchIO.cpp,v $   $Author: jnkim $
 * $Revision: 2756 $   $Date: 2008-06-23 14:09:25 -0500 (Mon, 23 Jun 2008) $
 * $Id: BranchIO.cpp 2756 2008-06-23 19:09:25Z jnkim $ 
 ***************************************************************************/

