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
ObservableHelper::ObservableHelper(hdf_path title) : group_name(std::move(title)) {}

ObservableHelper::~ObservableHelper() = default;

void ObservableHelper::set_dimensions(const std::vector<int>& dims, int first)
{
  //rank is increased
  hsize_t rank = dims.size() + 1;
  mydims.resize(rank, 1);
  copy(dims.begin(), dims.end(), mydims.begin() + 1);
  offsets.resize(rank, 0);
  lower_bound = first;
}

void ObservableHelper::addProperty(float& p, const std::string& pname, hdf_archive& file)
{
  double p_DP(p);
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::addProperty(Tensor<float, OHMMS_DIM>& p, const std::string& pname, hdf_archive& file)
{
  Tensor<double, OHMMS_DIM> p_DP(p);
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::addProperty(Matrix<float>& p, const std::string& pname, hdf_archive& file)
{
  Matrix<double> p_DP;
  p_DP = p;
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::addProperty(TinyVector<float, OHMMS_DIM>& p, const std::string& pname, hdf_archive& file)
{
  TinyVector<double, OHMMS_DIM> p_DP(p);
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::addProperty(std::vector<float>& p, const std::string& pname, hdf_archive& file)
{
  std::vector<double> p_DP;
  p_DP.assign(p.begin(), p.end());
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::addProperty(std::vector<TinyVector<float, OHMMS_DIM>>& p,
                                   const std::string& pname,
                                   hdf_archive& file)
{
  std::vector<TinyVector<double, OHMMS_DIM>> p_DP;
  p_DP.assign(p.begin(), p.end());
  file.push(group_name, true);
  file.write(p_DP, pname);
  file.pop();
}

void ObservableHelper::write(const value_type* const first_v, hdf_archive& file)
{
  hsize_t rank = mydims.size();
  if (rank)
  {
    file.push(group_name, true);
    h5d_append(file.top(), "value", current, rank, mydims.data(), first_v + lower_bound);
    file.pop();
  }
}

} // namespace qmcplusplus
