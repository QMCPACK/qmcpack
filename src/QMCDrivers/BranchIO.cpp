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
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "Message/CommOperators.h"

//#include <boost/archive/text_oarchive.hpp>

namespace qmcplusplus 
{
  bool BranchIO::write(const string& fname) 
  {
    hid_t fid =  H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    //cannot find the file, return false
    if(fid<0) return false;

    hid_t h1 =  H5Gopen(fid,hdf::main_state),h2;
    herr_t status = H5Eset_auto(NULL, NULL);
    status = H5Gget_objinfo(h1,hdf::qmc_status,0,NULL);
    bool overwrite=(status == 0);
    if(overwrite)
    {
      h2=H5Gopen(h1,hdf::qmc_status);
    }
    else
    {
      h2=H5Gcreate(h1,hdf::qmc_status,0);

      string v_header("tau:taueff:etrial:eref:branchmax:branchcutoff:branchfilter:sigma:acc_energy:acc_samples");
      HDFAttribIO<string> v_out(v_header);
      v_out.write(h2,"vparam_def");

      string i_header("warmupsteps:energyupdateinterval:counter:targetwalkers:maxwalkers:minwalkers:branchinterval");
      HDFAttribIO<string> i_out(i_header);
      i_out.write(h2,"iparam_def");
    }

    //\since 2008-06-24
    HDFAttribIO<VParamType> vh(ref.vParam,overwrite);
    vh.write(h2,"vparam");

    HDFAttribIO<IParamType> ih(ref.iParam,overwrite);
    ih.write(h2,"iparam");

    HDFAttribIO<BranchModeType> qstatus(ref.BranchMode,overwrite);
    qstatus.write(h2,"branchmode");

    //if(LogNorm.size())//check if collection is done correctly
    //{
    //  HDFAttribIO<vector<RealType> > lh(LogNorm,overwrite);
    //  lh.write(h1,hdf::norm_history);
    //}
    H5Gclose(h2);
    H5Gclose(h1);
    H5Fclose(fid);

    return true;
  }

  bool BranchIO::read(const string& fname) 
  {
    hid_t h_file =  H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
    if(h_file<0)
    {
      app_error() << "  Failed to open " << fname << endl;
      return false;
    }

    HDFVersion res_version(0,4); //start using major=0 and minor=4
    HDFVersion res_20080624(0,5);//major revision on 2008-06-24 0.5
    HDFVersion in_version(0,1);
    vector<RealType> pdata;

    herr_t status = H5Eset_auto(NULL, NULL);
    in_version.read(h_file,hdf::version);
    if(in_version>= res_20080624)
    {//need a better version control
      hid_t h1=H5Gopen(h_file,hdf::main_state);
      hid_t h2=H5Gopen(h1,hdf::qmc_status);
      //branchmode
      HDFAttribIO<BranchModeType> br(ref.BranchMode);
      br.read(h2,"branchmode");
      HDFAttribIO<VParamType> vh(ref.vParam);
      vh.read(h2,"vparam");
      //iParam
      HDFAttribIO<IParamType> ih(ref.iParam);
      ih.read(h2,"iparam");

      //pack these to bcast
      pdata.insert(pdata.end(),ref.vParam.begin(),ref.vParam.end());
      pdata.insert(pdata.end(),ref.iParam.begin(),ref.iParam.end());
      pdata.insert(pdata.end(),ref.BranchMode.to_ulong());
      H5Gclose(h2);
      H5Gclose(h1);
    }
    else if(in_version>=res_version)
    {
      //in_version.read(h_file,hdf::version);
      hid_t h1=H5Gopen(h_file,hdf::main_state);
      HDFAttribIO<vector<RealType> > eh(pdata);
      eh.read(h1,hdf::energy_history);
      //if(LogNorm.size()) 
      //{
      //  HDFAttribIO<vector<RealType> > lh(LogNorm);
      //  lh.read(h1,hdf::norm_history);
      //}
      H5Gclose(h1);
    }
    else 
    { 
      app_log() << "  Missing version. Using old format " << endl;
      hid_t h1 = H5Gopen(h_file,"config_collection");
      herr_t status = H5Eset_auto(NULL, NULL);
      status = H5Gget_objinfo (h1, "Summary", 0, NULL);
      if(status == 0) {
        HDFAttribIO<vector<RealType> > eh(pdata);
        eh.read(h1,"Summary");
        //if(LogNorm.size()) 
        //{
        //  HDFAttribIO<vector<RealType> > lh(LogNorm);
        //  lh.read(h1,"LogNorm");
        //  app_log() << " Normalization factor defined from previous calculation " << endl;
        //}
      } 
      H5Gclose(h1);
    }
    H5Fclose(h_file);

    //broadcast to the nodes : need to add a namespace mpi::
    myComm->bcast(pdata);
    //\since 2008-06-24 
    if(in_version>=res_20080624)
    {
      int ii=0;
      for(int i=0; i<ref.vParam.size(); ++i,++ii) ref.vParam[i]=pdata[ii];
      for(int i=0; i<ref.iParam.size(); ++i,++ii) ref.iParam[i]=static_cast<int>(pdata[ii]);
      ref.BranchMode=static_cast<unsigned long>(pdata[ii]);
    }
    else
    {
      ref.vParam[SimpleFixedNodeBranch::B_ETRIAL]=pdata[0];
      ref.vParam[SimpleFixedNodeBranch::B_EREF]=pdata[0];
    }

    return true;
  }
}

/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.cpp,v $   $Author: jnkim $
 * $Revision: 2756 $   $Date: 2008-06-23 14:09:25 -0500 (Mon, 23 Jun 2008) $
 * $Id: SimpleFixedNodeBranch.cpp 2756 2008-06-23 19:09:25Z jnkim $ 
 ***************************************************************************/

