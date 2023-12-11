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


namespace qmcplusplus
{
// Store the input data and gold data for testing the overlap and Hamiltonian
// matrices for the linear method optimizer.

struct FillData
{
  int numSamples;
  int numParam;
  double sum_wgt;
  double sum_e_wgt;
  double sum_esq_wgt;

  std::vector<double> reweight;
  std::vector<double> energy_new;

  Matrix<double> derivRecords;
  Matrix<double> HDerivRecords;

  // Gold data for comparison
  Matrix<double> ovlp_gold;
  Matrix<double> ham_gold;
};


} // namespace qmcplusplus
