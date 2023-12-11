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


/** @file OutputMatrix.h
 * @brief Output the Hamiltonian and overlap matrices from linear method
 */
#ifndef QMCPLUSPLUS_OUTPUT_MATRIX_H
#define QMCPLUSPLUS_OUTPUT_MATRIX_H

#include <ostream>
#include <string>
#include "Configuration.h"
#include "OhmmsPETE/OhmmsMatrix.h"


namespace qmcplusplus
{
class OutputMatrix
{
public:
  using RealType = QMCTraits::RealType;

  ///Constructor.
  OutputMatrix();

  /// Open a text-formatted scalar.dat file and print the header line
  void init_file(const std::string& root_name, const std::string& name, int N);

  /// Print matrix to text-formatted scalar.dat file
  void output(Matrix<RealType>& mat);

private:
  int index_;

  std::ofstream output_file_;
};
} // namespace qmcplusplus
#endif
