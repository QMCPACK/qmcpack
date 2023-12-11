//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file OutputMatrix.cpp
 * @brief Output the Hamiltonian and overlap matrices from linear method
 */

#include "QMCDrivers/WFOpt/OutputMatrix.h"


namespace qmcplusplus
{
OutputMatrix::OutputMatrix() : index_(0) {}

void OutputMatrix::init_file(const std::string& root_name, const std::string& name, int N)
{
  int namelen = root_name.length();
  // Assume that the root_name has the series suffix (".s000").
  // Remove the series suffix to get the project id
  std::string fname = root_name.substr(0, namelen - 5) + "." + name + ".s000.scalar.dat";
  app_log() << "Output matrix file: " << fname << std::endl;

  output_file_.open(fname);
  output_file_ << "# Index ";
  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < N; j++)
    {
      std::string index_name = name + "_" + std::to_string(i) + "_" + std::to_string(j);
      output_file_ << index_name << " ";
    }
  }
  output_file_ << std::endl;
}

void OutputMatrix::output(Matrix<RealType>& mat)
{
  output_file_ << index_ << " ";
  for (int i = 0; i < mat.rows(); i++)
  {
    for (int j = 0; j < mat.cols(); j++)
    {
      output_file_ << mat(i, j) << " ";
    }
  }
  output_file_ << std::endl;
  index_++;
}


} // namespace qmcplusplus
