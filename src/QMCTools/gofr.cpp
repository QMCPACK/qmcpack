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
using namespace std;

  template<typename IT1, typename IT2>
inline void accumulate_squaredelements(IT1 first, IT1 last, IT2 res)
{
  while(first != last) {*res++ += (*first)*(*first);++first;}
}

int print_help()
{
  cout << "Usage: gofr --fileroot string --np int " << endl;
  return 1;
}

struct GofRObserver {
  int NumSamples;
  string DataSetName;
  Vector<double> gofr;
  Vector<double> gofr2;
  Matrix<double> gofr_dat;

  GofRObserver(const string& aname):DataSetName(aname),NumSamples(0)
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
      //accumulate_elements(gofr_dat[i],gofr_dat[i]+nbin,gofr.begin());
      //accumulate_squaredelements(gofr_dat[i],gofr_dat[i]+nbin,gofr2.begin());
    }
    NumSamples += count;
  }

  inline void resize(int nbin)
  {
    if(gofr.size()<nbin)
    {
      gofr.resize(nbin);
      gofr2.resize(nbin);
      gofr=0.0;
      gofr2=0.0;
    }
  }

  void getData(const char* fname, int nproc);

  void print(const char* fname)
  {
    ofstream fout(fname);
    fout.setf(ios::scientific, ios::floatfield);
    fout.setf(ios::left,ios::adjustfield);
    fout.precision(6);

    double norm=1.0/static_cast<double>(NumSamples);
    double sqrtnorm=sqrt(norm);
    for(int i=0; i<gofr.size(); i++)
    {
      double avg=gofr[i]*norm;
      double var= sqrt(gofr2[i]*norm-avg*avg);
      fout << setw(3) << i <<  setw(20) << avg << setw(20) << var*sqrtnorm << setw(20) << var << endl;
    }
  }
};

int main(int argc, char** argv) {

  if(argc<2) return print_help();

  int iargc=0;
  int nproc=1;
  vector<string> gofr_name;
  string h5fileroot("0");
  while(iargc<argc)
  {
    string c(argv[iargc]);
    if(c == "--fileroot") 
    {
      h5fileroot=argv[++iargc];
    }
    else if(c == "--np" || c == "-np")
    {
      nproc=atoi(argv[++iargc]);
    }
    else if(c == "--help" || c == "-h")
    {
      return print_help();
    }
    else if(c == "--dataset" || c == "-d")
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
  char fname[128];
  for(int ip=0; ip<nproc; ip++)
  {
    if(nproc>1)
      sprintf(fname,"%s.p%03d.config.h5",froot,ip);
    else
      sprintf(fname,"%s.config.h5",froot);

    cout << "Getting data from " << fname << endl;

    hid_t fileid =  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t obsid = H5Gopen(fileid,"observables");

    int count=0;
    HDFAttribIO<int> i(count);
    i.read(obsid,"count");

    string gname(DataSetName);
    gname.append("/value");
    //check dimension
    hid_t h1 = H5Dopen(obsid, gname.c_str());
    hid_t dataspace = H5Dget_space(h1);
    hsize_t dims_out[2];
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    H5Dclose(h1);

    gofr_dat.resize(dims_out[0],dims_out[1]);
    HDFAttribIO<Matrix<double> > o(gofr_dat);
    resize(dims_out[1]);
    o.read(obsid,gname.c_str());
    accumulate(count);

    H5Gclose(obsid);
    H5Fclose(fileid);
  }

  string outfname(froot);
  outfname.append(".gofr.dat");
  print(outfname.c_str());
}
