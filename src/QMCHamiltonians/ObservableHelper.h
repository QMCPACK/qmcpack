//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    William F. Godoy, godoywf@ornl.gov,  Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file ObservableHelper.h
 *@brief Declaration of ObservableHelper and other helper class for observables
 */
#ifndef QMCPLUSPLUS_OBSERVABLEHELPER_H
#define QMCPLUSPLUS_OBSERVABLEHELPER_H

#include "Configuration.h"
#include "hdf/hdf_archive.h"
#include "hdf/hdf_path.h"

namespace qmcplusplus
{
using value_type = QMCTraits::FullPrecRealType;

/** define ObservableHelper
 *
 * This is a helper class to manage a hdf5 dagroup for each collectable.
 * The data handled by an Estimator should be presented by dense N-dim array in C.
 * /observables/title/value
 */
class ObservableHelper
{
public:
  /**
   * Favored constructor
   * \param[in] title         is the ordered hdf5 group path elements of the observable
   */
  ObservableHelper(hdf_path title);

  /**
   * delete copy constructor as hdf5 handlers must have unique owners
   */
  ObservableHelper(const ObservableHelper&)            = delete;
  ObservableHelper& operator=(const ObservableHelper&) = delete;

  /**
   * Move constructor.
   * @param in input object to be moved to this
   */
  ObservableHelper(ObservableHelper&&) noexcept            = default;
  ObservableHelper& operator=(ObservableHelper&&) noexcept = default;

  /**
   * Destructor closes hdf5 remaining resources
   */
  ~ObservableHelper();

  /**
   * set the shape of this observable
   * @param dims dimensions
   * @param first starting index
   */
  void set_dimensions(const std::vector<int>& dims, int first);

  /** add named property to describe the collectable this helper class handles
   * @param p any intrinsic datatype including vector, basic containers
   * @param pname property
   */
  template<typename T>
  inline void addProperty(T& p, const std::string& pname, hdf_archive& file)
  {
    file.push(group_name, true);
    file.write(p, pname);
    file.pop();
  }

  void addProperty(float& p, const std::string& pname, hdf_archive& file);
  void addProperty(Tensor<float, OHMMS_DIM>& p, const std::string& pname, hdf_archive& file);
  void addProperty(Matrix<float>& p, const std::string& pname, hdf_archive& file);
  void addProperty(TinyVector<float, OHMMS_DIM>& p, const std::string& pname, hdf_archive& file);
  void addProperty(std::vector<float>& p, const std::string& pname, hdf_archive& file);
  void addProperty(std::vector<TinyVector<float, OHMMS_DIM>>& p, const std::string& pname, hdf_archive& file);

  void write(const value_type* const first_v, hdf_archive& file);

  ///starting index
  hsize_t lower_bound = 0;

private:
  /// "file pointer" for h5d_append
  hsize_t current = 0;
  /// Path of this observable
  hdf_path group_name;
  ///my dimensions
  std::vector<hsize_t> mydims;
  ///offsets
  std::vector<hsize_t> offsets;
};
} // namespace qmcplusplus
#endif
