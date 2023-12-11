//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <sstream>
#include <array>
#include "Message/Communicate.h"
#include "hdf/hdf_archive.h"

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
#endif
  using namespace qmcplusplus;
  if (argc != 2)
  {
    std::cout << "Program must take as an argument a single input eshdf file to parse" << std::endl;
    return 1;
  }
  // Suppress HDF5 warning and error messages.
  qmcplusplus::hdf_error_suppression hide_hdf_errors;
  std::string fname(argv[1]);
  //cout << "Input file is: " << fname << std::endl;
  hdf_archive hin;
  hin.open(fname, H5F_ACC_RDONLY);
  const std::string dataset = "electrons/number_of_kpoints";
  // Get the number of kpoints from the eshdf file
  int data;
  hin.read(data, dataset);
  // get all of the reduced kpoints from the file
  std::array<float, 3> kpt;
  for (int i = 0; i < data; i++)
  {
    std::ostringstream os;
    os << "electrons/kpoint_" << i << "/reduced_k";
    hin.read(kpt, os.str());
    std::cout << kpt[0] << "   " << kpt[1] << "   " << kpt[2] << "\n";
  }
  hin.close();
  return 0;
}
