#include "hdf5.h"
#include <iostream>
#include <string>
#include <sstream>
using namespace std;



int main (int argc, char* argv[])
{
  hid_t           file, dset;    /* Handles */
  herr_t          status;
  if (argc != 2)
  {
    cout << "Program must take as an argument a single input eshdf file to parse" << endl;
    return 1;
  }
  string fname(argv[1]);
  //cout << "Input file is: " << fname << endl;
  string dataset = "electrons/number_of_kpoints";
  file = H5Fopen (fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  // Get the number of kpoints from the eshdf file
  dset = H5Dopen (file, dataset.c_str());
  int data;
  status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
  status = H5Dclose (dset);
  // get all of the reduced kpoints from the file
  float kpt[3];
  for (int i = 0; i < data; i++)
  {
    std::ostringstream os;
    os << "electrons/kpoint_" << i << "/reduced_k";
    dset = H5Dopen (file, os.str().c_str());
    status = H5Dread (dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &kpt);
    cout << kpt[0] << "   " << kpt[1] << "   " << kpt[2] << "\n";
    status = H5Dclose(dset);
  }
  status = H5Fclose (file);
  return 0;
}
