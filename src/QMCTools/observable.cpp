#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <OhmmsData/OhmmsElementBase.h>
#include <Numerics/HDFSTLAttrib.h>
using namespace qmcplusplus;
using namespace std;

int print_help()
{
  cerr << "Usage: observable [--fileroot string ] list-of-obervables " << endl;
  cerr << " Example: observable --fileroot myhdf --first 10 density " << endl;
  return 1;
}

struct ObservableFile
{
  int first_block;
  string h5name;
  vector<string> olist;
  ObservableFile(int argc, char** argv);
  int dump_averages();

  template<typename IT>
  inline  void accumulate_n(IT first, IT last, IT result)
  {
    for(; first<last;)
      *result++ += *first++;
  }

  template<typename IT, typename T>
  inline  void scale_n(IT first, IT last, T fac)
  {
    for(; first<last;)
      *first++ *=fac;
  }
};

int main(int argc, char** argv)
{
  if(argc<2)
    return print_help();
  ObservableFile in(argc,argv);
  return in.dump_averages();
}

ObservableFile::ObservableFile(int argc, char** argv)
  :first_block(0)
{
  int iargc=1;
  while(iargc<argc)
  {
    string c(argv[iargc]);
    if(c.find("fileroot")<c.size())
    {
      h5name=argv[++iargc];
    }
    else
      if(c.find("first")<c.size())
      {
        first_block=atoi(argv[++iargc]);
      }
      else
        if(c.find("help")<c.size())
        {
          print_help();
        }
        else
        {
          olist.push_back(argv[iargc]);
        }
    ++iargc;
  }
}

int ObservableFile::dump_averages()
{
  if(h5name.empty())
  {
    cerr << "  Invalid file. Provide --fileroot root-of-h5-file " << endl;
    return 1;
  }
  herr_t status = H5Eset_auto(NULL, NULL);
  char fname[256];
  if(h5name.find("h5")<h5name.size())
    sprintf(fname,"%s",h5name.c_str());
  else
    sprintf(fname,"%s.h5",h5name.c_str());
  hid_t fileID =  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
  if(fileID<0)
  {
    cerr << " Cannot find " << fname << endl;
    return 1;
  }
  if(olist.empty())
  {
    cerr << "No observables given." << endl;
    print_help();
    return 1;
  }
  for(int i=0; i<olist.size(); ++i)
  {
    hid_t gid=H5Gopen(fileID,olist[i].c_str());
    if(gid<0)
    {
      cerr << "Cannot find " << olist[i] << " group. Abort" << endl;
      return 1;
    }
    hid_t dataset = H5Dopen(gid,"value");
    hid_t dataspace = H5Dget_space(dataset);
    //check the dimension of the active observable
    int rank = H5Sget_simple_extent_ndims(dataspace);
    vector<hsize_t> in_dims(rank);
    int status_n = H5Sget_simple_extent_dims(dataspace, &in_dims[0], NULL);
    vector<int> dims(rank-1);
    int num_blocks=static_cast<int>(in_dims[0]);
    int nb=num_blocks-first_block;
    //skip it
    if(nb<0)
      continue;
    cout << "<<<< Input dimension  " << num_blocks;
    int tot_size=1;
    for(int dim=0; dim<rank-1; ++dim)
    {
      tot_size *= dims[dim]=static_cast<int>(in_dims[dim+1]);
      cout << " "<< dims[dim] ;
    }
    cout << "\n Total size = " << tot_size << endl;
    vector<double> odata_in;
    odata_in.resize(tot_size*num_blocks);
    hid_t ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(odata_in[0]));
    H5Dclose(dataset);
    H5Sclose(dataspace);
    vector<double> odata_tot(tot_size,0.0), odata2_tot(tot_size,0.0);
    vector<double>::iterator first=odata_in.begin()+first_block*tot_size;
    for(int b=first_block; b<num_blocks; ++b,first+=tot_size)
      accumulate_n(first,first+tot_size,odata_tot.begin());
    double fac=1.0/static_cast<double>(nb);
    scale_n(odata_tot.begin(),odata_tot.end(),fac);
    HDFAttribIO<vector<double> > dummy(odata_tot,dims);
    dummy.write(gid,"average");
  }
  return 0;
}
