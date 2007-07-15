#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/HDFAttribIO.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"

int main(int argc, char** argv) {

  int iargc=0;
  int nproc=1;
  string h5fileroot("0");
  while(iargc<argc)
  {
    string c(argv[iargc]);
    if(c == "--fileroot") 
    {
      h5fileroot=argv[++iargc];
    }
    else if(c == "--np")
    {
      iproc=atoi(argv[++iargc]);
    }
    else
      ++iargc;
  }

  double rmin=-10;
  double rmax=10;
  double delta=0.05;
  int npts=static_cast<int>(rmax-rmin)/delta+1;
  Matrix<double> den(npts,npts);
  typedef TinyVector<double,3> PosType;
  hsize_t dimIn[3],dimTot[3];
  char fname[256],coordname[128];
  for(int ip=0; ip<np; ip++)
  {
    if(np>1)
      sprintf(fname,"%s.p%03d",h5fileroot.c_str(),ip);
    else
      sprintf(fname,"%s",h5fileroot.c_str());

    hid_t fileID =  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    int nconf=1;
    hid_t mastercf = H5Gopen(fileID,"config_collection");
    hid_t h1=H5Dopen(mastercf,"NumOfConfigurations");
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(nconf));
    H5Dclose(h1);

    for(int iconf=0; iconf<nconf; iconf++) { 
      sprintf(coordName,"config%04d/coord",iconf); 
      hid_t dataset = H5Dopen(mastercf,coordName);
      hid_t dataspace = H5Dget_space(dataset);
      int rank = H5Sget_simple_extent_ndims(dataspace);
      int status_n = H5Sget_simple_extent_dims(dataspace, dimTot, NULL);
      cout << "Going to read " << dimTot[0] << " " << dimTot[1] << endl;
      vector<PosType> pos(dimTot[0]*dimTot[1]);
      hid_t ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(ref[0][0]));
      H5Dclose(dataset);
      H5Sclose(dataspace);

      for(int iat=0; iat<pos.size(); iat++)
      {
        cout << pos[iat] << endl;
      }
    }
  }
  return 0;
}
