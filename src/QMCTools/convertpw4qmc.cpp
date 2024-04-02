/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <fstream>
#include "config.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#include "XmlRep.h"
#include "WriteEshdf.h"

using namespace std;

void convertToVecStrings(int argc, char* argv[], std::vector<std::string>& vec)
{
  for (int i = 1; i < argc; i++)
  {
    vec.push_back(argv[i]);
  }
}


inline bool file_exists(const std::string& name)
{
  ifstream ifile(name.c_str());
  return (bool)ifile;
}

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  int rank = 0, world_size = 1;
  int provided;
  MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

  // Suppress HDF5 warning and error messages.
  qmcplusplus::hdf_error_suppression hide_hdf_errors;

  vector<string> vecParams;
  convertToVecStrings(argc, argv, vecParams);
  string fname;
  string ofname = "eshdf.h5";
  for (int i = 0; i < vecParams.size(); i++)
  {
    if (vecParams[i] == "--output" || vecParams[i] == "-o")
    {
      ofname = vecParams[i + 1];
      i++;
    }
    else
    {
      fname = vecParams[i];
    }
  }
  if (file_exists(fname) == false)
  {
    cerr << "must specify a valid xml file name as an argument" << endl;
    cerr << fname << " did not work" << endl;
    return (1);
  }

  ifstream is(fname.c_str());

  XmlNode inFile(&is, 0, true);
  string name = inFile.getName();
  string directory("./");
  const int last_slash_index = fname.find_last_of('/');
  if (last_slash_index > -1)
  {
    directory = fname.substr(0, last_slash_index + 1);
  }

  if (name == "qes:espresso")
  {
    // may be good to test different versions of espresso to see
    // if this changes at all
    cout << "xml file comes from quantum espresso" << endl;
    EshdfFile outFile(ofname);
    outFile.writeQEBoilerPlate(inFile);
    outFile.writeQESupercell(inFile);
    outFile.writeQEAtoms(inFile);
    outFile.writeQEElectrons(inFile, directory);
  }
  else if (name == "fpmd:sample")
  {
    cout << "xml file comes from qbox" << endl;
    EshdfFile outFile(ofname);
    outFile.writeQboxBoilerPlate(inFile);
    outFile.writeQboxSupercell(inFile);
    outFile.writeQboxAtoms(inFile);
    outFile.writeQboxElectrons(inFile);
  }
  else
  {
    cerr << "xml format is unrecognized" << endl;
    return (1);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
