//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "QMCHamiltonians/ObservableHelper.h"
#include "io/hdf/hdf_archive.h"

#include <filesystem>

/*
  -- 05/07/2021 --
  Deterministic unit test for the ObservableHelper class.
*/

namespace qmcplusplus
{

TEST_CASE("ObservableHelper::ObservableHelper(std::vector<std::string>)", "[hamiltonian]")
{
  ObservableHelper oh(hdf_path{"u/v"});
  CHECK(oh.lower_bound == 0);
}

TEST_CASE("ObservableHelper::set_dimensions", "[hamiltonian]")
{
  ObservableHelper oh{hdf_path{"u"}};

  std::vector<int> dims = {10, 10};
  oh.set_dimensions(dims, 1);

  CHECK(oh.lower_bound == 1);
}

TEST_CASE("ObservableHelper::ObservableHelper()", "[hamiltonian]")
{
  std::filesystem::path filename("tmp_ObservableHelper2.h5");
  hdf_archive hFile;
  hFile.create(filename);

  ObservableHelper oh{hdf_path{"u"}};
  std::vector<int> dims = {10, 10};
  float propertyFloat   = 10.f;
  oh.addProperty(propertyFloat, "propertyFloat", hFile);

  Tensor<float, OHMMS_DIM> propertyTensor;
  oh.addProperty(propertyTensor, "propertyTensor", hFile);

  Matrix<float> propertyMatrix;
  oh.addProperty(propertyMatrix, "propertyMatrix", hFile);

  TinyVector<float, OHMMS_DIM> propertyTinyVector;
  oh.addProperty(propertyTensor, "propertyTinyVector", hFile);

  std::vector<float> propertyVector;
  oh.addProperty(propertyVector, "propertyVector", hFile);

  std::vector<TinyVector<float, OHMMS_DIM>> propertyVectorTinyVector;
  oh.addProperty(propertyVectorTinyVector, "propertyVectorTinyVector", hFile);

  hFile.close();
  REQUIRE(std::filesystem::exists(filename));
  REQUIRE(std::filesystem::remove(filename));
}

} // namespace qmcplusplus
