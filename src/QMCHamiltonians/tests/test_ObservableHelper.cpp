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
TEST_CASE("ObservableHelper::ObservableHelper(const std::string&)", "[hamiltonian]")
{
  ObservableHelper oh("u");
  CHECK(oh.lower_bound == 0);
  CHECK(oh.mydims == std::vector<hsize_t>());
  CHECK(oh.maxdims == std::vector<hsize_t>());
  CHECK(oh.curdims == std::vector<hsize_t>());
  CHECK(oh.offsets == std::vector<hsize_t>());
  CHECK(oh.group_name == "u");
  CHECK(oh.isOpened() == false);
}

TEST_CASE("ObservableHelper::ObservableHelper(ObservableHelper&&)", "[hamiltonian]")
{
  std::filesystem::path filename("tmp_ObservableHelper1.h5");
  hdf_archive hFile;
  hFile.create(filename);

  ObservableHelper ohIn("u");
  ohIn.open(hFile);
  CHECK(ohIn.isOpened() == true);
  {
    ObservableHelper oh(std::move(ohIn));
    CHECK(oh.isOpened() == true);
    CHECK(ohIn.isOpened() == false);
  }
  CHECK(ohIn.isOpened() == false);

  hFile.close();
  REQUIRE(std::filesystem::exists(filename));
  REQUIRE(std::filesystem::remove(filename));
}

TEST_CASE("ObservableHelper::set_dimensions", "[hamiltonian]")
{
  ObservableHelper oh("u");

  std::vector<int> dims = {10, 10};
  oh.set_dimensions(dims, 1);

  CHECK(oh.mydims == std::vector<hsize_t>{1, 10, 10});
  CHECK(oh.curdims == std::vector<hsize_t>{1, 10, 10});
  CHECK(oh.maxdims == std::vector<hsize_t>{H5S_UNLIMITED, 10, 10});
  CHECK(oh.offsets == std::vector<hsize_t>{0, 0, 0});
  CHECK(oh.lower_bound == 1);
}

TEST_CASE("ObservableHelper::ObservableHelper()", "[hamiltonian]")
{
  std::filesystem::path filename("tmp_ObservableHelper2.h5");
  hdf_archive hFile;
  hFile.create(filename);

  ObservableHelper oh("u");
  oh.open(hFile);
  std::vector<int> dims = {10, 10};
  float propertyFloat   = 10.f;
  oh.addProperty(propertyFloat, "propertyFloat");

  Tensor<float, OHMMS_DIM> propertyTensor;
  oh.addProperty(propertyTensor, "propertyTensor");

  Matrix<float> propertyMatrix;
  oh.addProperty(propertyMatrix, "propertyMatrix");

  TinyVector<float, OHMMS_DIM> propertyTinyVector;
  oh.addProperty(propertyTensor, "propertyTinyVector");

  std::vector<float> propertyVector;
  oh.addProperty(propertyVector, "propertyVector");

  std::vector<TinyVector<float, OHMMS_DIM>> propertyVectorTinyVector;
  oh.addProperty(propertyVectorTinyVector, "propertyVectorTinyVector");

  hFile.close();
  REQUIRE(std::filesystem::exists(filename));
  REQUIRE(std::filesystem::remove(filename));
}

} // namespace qmcplusplus
