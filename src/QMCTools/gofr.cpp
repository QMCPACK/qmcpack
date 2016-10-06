//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "Configuration.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/IteratorUtility.h"
using namespace qmcplusplus;

int print_help()
{
  std::cout << "Usage: gofr --fileroot std::string --np int " << std::endl;
  return 1;
}

struct GofRObserver
{
  int NumSamples;
  int NumNodes;
  std::string DataSetName;
  Vector<double> gofr;
  Vector<double> gofr2;
  Matrix<double> gofr_dat;

  GofRObserver(const std::string& aname):DataSetName(aname),NumSamples(0)
  {
  }

  inline void accumulate(int count)
  {
    int nbin=gofr_dat.cols();
    for(int i=0; i<count; i++)
    {
      double* restrict dptr=gofr_dat[i];
      for(int k=0; k<nbin; k++,dptr++)
      {
        gofr[k] += *dptr;
        gofr2[k] += (*dptr)*(*dptr);
      }
    }
    NumSamples = count;
  }

  inline void resize(int ns, int nbin)
  {
    NumSamples=ns;
    gofr_dat.resize(ns,nbin);
    gofr_dat=0.0;
    gofr.resize(nbin);
    gofr2.resize(nbin);
    gofr=0.0;
    gofr2=0.0;
  }

  void getData(const char* fname, int nproc);

  void print(const char* fname)
  {
    std::ofstream fout(fname);
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout.setf(std::ios::left,std::ios::adjustfield);
    fout.precision(6);
    //double norm=static_cast<double>(NumNodes)/static_cast<double>(NumSamples);
    double norm=1.0/static_cast<double>(NumSamples);
    double sqrtnorm=sqrt(norm);
    for(int i=0; i<gofr.size(); i++)
    {
      double avg=gofr[i]*norm;
      double var= sqrt(gofr2[i]*norm-avg*avg);
      fout << std::setw(3) << i <<  std::setw(20) << avg << std::setw(20) << var*sqrtnorm << std::setw(20) << var << std::endl;
    }
  }
};

int main(int argc, char** argv)
{
  if(argc<2)
    return print_help();
  int iargc=0;
  int nproc=1;
  std::vector<std::string> gofr_name;
  std::string h5fileroot("0");
  while(iargc<argc)
  {
    std::string c(argv[iargc]);
    if(c == "--fileroot")
    {
      h5fileroot=argv[++iargc];
    }
    else
      if(c == "--np" || c == "-np")
      {
        nproc=atoi(argv[++iargc]);
      }
      else
        if(c == "--help" || c == "-h")
        {
          return print_help();
        }
        else
          if(c == "--dataset" || c == "-d")
          {
            gofr_name.push_back(argv[++iargc]);
          }
    ++iargc;
  }
  for(int i=0; i<gofr_name.size(); i++)
  {
    GofRObserver recorder(gofr_name[i]);
    recorder.getData(h5fileroot.c_str(),nproc);
  }
  return 0;
}

void GofRObserver::getData(const char* froot, int nproc)
{
  NumNodes=nproc;
  char fname[128];
  int count_max=1000000;
  for(int ip=0; ip<nproc; ip++)
  {
    if(nproc>1)
      sprintf(fname,"%s.p%03d.config.h5",froot,ip);
    else
      sprintf(fname,"%s.config.h5",froot);
    std::cout << "Getting data from " << fname << std::endl;
    hid_t fileid =  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t obsid = H5Gopen(fileid,"observables");
    int count=0;
    HDFAttribIO<int> i(count);
    i.read(obsid,"count");
    count_max = std::min(count_max,count);
    std::string gname(DataSetName);
    gname.append("/value");
    //check dimension
    hid_t h1 = H5Dopen(obsid, gname.c_str());
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[2];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    H5Dclose(h1);
    if(gofr_dat.size() ==0)
    {
      resize(dims_out[0],dims_out[1]);
    }
    Matrix<double> tmp(dims_out[0],dims_out[1]);
    HDFAttribIO<Matrix<double> > o(tmp);
    o.read(obsid,gname.c_str());
    gofr_dat += tmp;
    H5Gclose(obsid);
    H5Fclose(fileid);
  }
  std::cout << "Number of blocks " << count_max << std::endl;
  std::accumulate(count_max);
  std::string outfname(froot);
  outfname += "."+ DataSetName+".dat";
  print(outfname.c_str());
}
