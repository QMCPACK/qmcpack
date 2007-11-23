//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#include "Particle/HDFWalkerInput_0_4.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIOEngine.h"

namespace qmcplusplus
{

  HDFWalkerInput_0_4::HDFWalkerInput_0_4(MCWalkerConfiguration& W, Communicate* c, const HDFVersion& v)
    : targetW(W), myComm(c), cur_version(0,4)
    {
      i_info.version=v;
    }

  HDFWalkerInput_0_4::~HDFWalkerInput_0_4() {
    if(h_plist != H5P_DEFAULT) H5Pclose(h_plist);
  }

  void HDFWalkerInput_0_4::checkOptions(xmlNodePtr cur) 
  {
    i_info.reset();

    string froot,cfile;
    string collected("no");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(i_info.nprocs,"nprocs");
    pAttrib.add(i_info.rank,"node");
    pAttrib.add(cfile,"href"); pAttrib.add(cfile,"file"); 
    pAttrib.add(froot,"fileroot");
    pAttrib.add(collected,"collected");
    pAttrib.put(cur);

    if(froot.empty()) return;

    int ext=froot.find(hdf::config_ext);
    if(ext<froot.size())
    { //remove extenstion
      froot.erase(froot.begin()+ext,froot.end());
    }
    
    i_info.collected = (collected == "yes");
    if(i_info.collected)
      FileStack.push(froot);
    else
    {
      if(i_info.nprocs>1)//need to process multiple files
      {
        int nprocs_now=myComm->size();
        if(i_info.nprocs>nprocs_now)//using less nodes
        {
          int np=i_info.nprocs/nprocs_now;
          int pid=myComm->rank(), ip=0;
          while(pid<i_info.nprocs && ip<np)
          {
            char *h5name=new char[froot.size()+10];
            sprintf(h5name,"%s.p%03d",froot.c_str(),pid++);
            FileStack.push(h5name);
            delete [] h5name;
            pid += nprocs_now;
          }
        }
        else
        {
          int pid=myComm->rank()%i_info.nprocs;
          char *h5name=new char[froot.size()+10];
          sprintf(h5name,"%s.p%03d",froot.c_str(),pid);
          FileStack.push(h5name);
          delete [] h5name;
        }
      }
      else
      {
        FileStack.push(froot);
      }
    }
  }

  bool HDFWalkerInput_0_4::put(xmlNodePtr cur) 
  {
    checkOptions(cur);

    if(FileStack.empty()) 
    {
      app_error() << "  No valid input hdf5 is found." << endl;
      return false;
    }

    //use collective I/O  for parallel runs with a file 
    i_info.parallel = ((myComm->size()>1) && i_info.collected);
    h_plist=H5P_DEFAULT;

#if defined(H5_HAVE_PARALLEL)
    if(i_info.parallel)
    {
      app_log() << "   HDFWalkerInput_0_4::put in parallel mode " << endl;
      MPI_Info info=MPI_INFO_NULL;
      h_plist = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(h_plist,myComm->getMPI(),info);
    }
#endif

    while(FileStack.size())
    {
      FileName=FileStack.top();
      FileStack.pop();

      string h5name(FileName);
      h5name.append(hdf::config_ext);

      hid_t fid= H5Fopen(h5name.c_str(),H5F_ACC_RDONLY,h_plist);
      if(fid<0)
      {
        app_error() << "  HDFWalkerInput_0_4::put Cannot open " << h5name << endl;
        return false;
      }

      //check if hdf and xml versions can work together
      HDFVersion aversion;
      aversion.setPmode(i_info.parallel);
      aversion.read(fid,hdf::version);
      if(aversion < i_info.version)
      {
        app_error() << " Mismatched version. xml = " << i_info.version << " hdf = " << aversion << endl;
        H5Fclose(fid);
        return false;
      }

      hid_t h1 = H5Gopen(fid,hdf::main_state);
      if(h1<0)
      {
        app_error() << "  HDFWalkerInput_0_4::put Cannot open " << hdf::main_state << " in " << h5name << endl;
        H5Gclose(h1);
        H5Fclose(fid);
        return false;
      }

      HDFWalkerIOEngine win(targetW);
      if(i_info.parallel)
        win.readAll(h1,hdf::walkers,myComm);
      else
        win.read(h1,hdf::walkers);

      H5Gclose(h1);
      H5Fclose(fid);
    }


    //char fname[128];
    //sprintf(fname,"%s.p%03d.xyz",FileName.c_str(),myComm->rank());
    //ofstream fout(fname);
    //MCWalkerConfiguration::iterator it = targetW.begin();
    //int iw=0;
    //while(it != targetW.end())
    //{
    //  fout << iw << endl;
    //  MCWalkerConfiguration::Walker_t& thisWalker(**it);
    //  for(int iat=0; iat<targetW.getTotalNum(); ++iat)
    //    fout << thisWalker.R[iat] << endl;
    //  ++it; ++iw;
    //}

    return true;
  }

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2303 $   $Date: 2007-11-19 13:25:50 -0600 (Mon, 19 Nov 2007) $
 * $Id: HDFWalkerInput_0_4.cpp 2303 2007-11-19 19:25:50Z jnkim $ 
 ***************************************************************************/
