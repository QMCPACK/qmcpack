//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    William F. Godoy, godoywf@ornl.gov,  Oak Ridge National Laboratory
//
// File created by: William F. Godoy, godoywf@ornl.gov,  Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/**@file ObservableHelper.cpp
 *@brief Definition of ObservableHelper class
 */
#include "ObservableHelper.h"

namespace qmcplusplus
{
ObservableHelper::ObservableHelper(const std::string& title)
    : data_id(-1), space1_id(-1), value1_id(-1), group_name(title), isopened(false)
{}

ObservableHelper::ObservableHelper(ObservableHelper&& in) noexcept
    : lower_bound(in.lower_bound),
      data_id(in.data_id),
      space1_id(in.space1_id),
      value1_id(in.value1_id),
      mydims(in.mydims),
      maxdims(in.maxdims),
      curdims(in.curdims),
      offsets(in.offsets),
      group_name(in.group_name),
      isopened(in.isopened)
{
  in.isopened = false;
}

ObservableHelper& ObservableHelper::operator=(ObservableHelper&& in) noexcept
{
  if (this != &in)
  {
    *this       = std::move(in);
    in.isopened = false;
  }
  return *this;
}

ObservableHelper::~ObservableHelper()
{
  if (isopened)
  {
    close();
  }
}

void ObservableHelper::set_dimensions(std::vector<int>& dims, int first)
{
  //rank is increased
  hsize_t rank = dims.size() + 1;
  mydims.resize(rank, 1);
  copy(dims.begin(), dims.end(), mydims.begin() + 1);
  maxdims = mydims;
  curdims = mydims;
  offsets.resize(rank, 0);
  maxdims[0]  = H5S_UNLIMITED;
  lower_bound = first;
}

void ObservableHelper::open(hid_t grp_id)
{
  data_id      = H5Gcreate2(grp_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t rank = mydims.size();
  if (rank)
  {
    //create empty data to write something
    hsize_t nd = 1;
    for (int i = 1; i < rank; ++i)
    {
      nd *= mydims[i];
    }
    std::vector<value_type> zeros(nd, 0.0);
    hid_t p = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(p, rank, &mydims[0]);
    space1_id      = H5Screate_simple(rank, &mydims[0], &maxdims[0]);
    value1_id      = H5Dcreate(data_id, "value", H5T_NATIVE_DOUBLE, space1_id, p);
    hid_t memspace = H5Screate_simple(rank, &mydims[0], NULL);
    herr_t ret     = H5Dwrite(value1_id, H5T_NATIVE_DOUBLE, memspace, space1_id, H5P_DEFAULT, &zeros[0]);
    H5Sclose(memspace);
    H5Pclose(p);
  }
  isopened = true;
}

void ObservableHelper::addProperty(float& p, const std::string& pname)
{
  double p_DP(p);
  HDFAttribIO<double> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::addProperty(Tensor<float, OHMMS_DIM>& p, const std::string& pname)
{
  Tensor<double, OHMMS_DIM> p_DP;
  p_DP = p;
  HDFAttribIO<Tensor<double, OHMMS_DIM>> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::addProperty(Matrix<float>& p, const std::string& pname)
{
  Matrix<double> p_DP;
  p_DP = p;
  HDFAttribIO<Matrix<double>> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::addProperty(TinyVector<float, OHMMS_DIM>& p, const std::string& pname)
{
  TinyVector<double, OHMMS_DIM> p_DP(p);
  HDFAttribIO<TinyVector<double, OHMMS_DIM>> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::addProperty(std::vector<float>& p, const std::string& pname)
{
  std::vector<double> p_DP;
  p_DP.assign(p.begin(), p.end());
  HDFAttribIO<std::vector<double>> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::addProperty(std::vector<TinyVector<float, OHMMS_DIM>>& p, const std::string& pname)
{
  std::vector<TinyVector<double, OHMMS_DIM>> p_DP;
  p_DP.assign(p.begin(), p.end());
  HDFAttribIO<std::vector<TinyVector<double, OHMMS_DIM>>> a(p_DP);
  a.write(data_id, pname.c_str());
}

void ObservableHelper::write(const value_type* first_v, const value_type* /*first_vv*/)
{
  hsize_t rank = mydims.size();
  if (rank)
  {
    H5Sset_extent_simple(space1_id, rank, &curdims[0], &maxdims[0]);
    // According to the HDF5 manual (https://support.hdfgroup.org/HDF5/doc/RM/H5S/H5Sselect_hyperslab.htm)
    // , the fifth argument (count) means the number of hyper-slabs to select along each dimensions 
    // while the sixth argument (block) is the size of each hyper-slab.
    // To write a single hyper-slab of size counts in a dataset, we call
    // H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, ones.data(), counts);
    // The vector "ones" means we want to write one hyper-slab (block) along each dimensions.
    // The result is equivalent to calling 
    // H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    // , but it implies writing count hyper-slabs along each dimension and each hyper-slab is of size one.
    const std::vector<hsize_t> ones(rank, 1);
    H5Sselect_hyperslab(space1_id, H5S_SELECT_SET, &offsets[0], NULL, ones.data(), &mydims[0]);
    H5Dextend(value1_id, &curdims[0]);
    hid_t memspace = H5Screate_simple(rank, &mydims[0], NULL);
    herr_t ret     = H5Dwrite(value1_id, H5T_NATIVE_DOUBLE, memspace, space1_id, H5P_DEFAULT, first_v + lower_bound);
    H5Sclose(memspace);
    curdims[0]++;
    offsets[0]++;
  }
}

bool ObservableHelper::isOpened() const noexcept { return isopened; }

// PRIVATE functions

void ObservableHelper::close()
{
  if (isopened)
  {
    if (space1_id > -1)
    {
      H5Sclose(space1_id);
      space1_id = -1;
    }
    if (data_id > -1)
    {
      H5Gclose(data_id);
      data_id = -1;
    }
    isopened = false;
  }
}

} // namespace qmcplusplus
