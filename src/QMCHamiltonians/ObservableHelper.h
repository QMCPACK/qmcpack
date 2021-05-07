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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file ObservableHelper.h
 *@brief Declaration of ObservableHelper and other helper class for observables
 */
#ifndef QMCPLUSPLUS_OBSERVABLEHELPER_H
#define QMCPLUSPLUS_OBSERVABLEHELPER_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"

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
  ///starting index
  hsize_t lower_bound = 0;
  ///id of this observable
  hid_t data_id = -1;
  ///dataspace for value
  hid_t space1_id = -1;
  ///id of the value dataset
  hid_t value1_id = -1;
  ///my dimensions
  std::vector<hsize_t> mydims;
  ///maximum dimensions
  std::vector<hsize_t> maxdims;
  ///current dimensions while extending data
  std::vector<hsize_t> curdims;
  ///offsets
  std::vector<hsize_t> offsets;
  ///name of this observable
  std::string group_name;

  /**
   * default constructor
   */
  ObservableHelper(const std::string& title = "dummy");

  /**
   * delete copy constructor as hdf5 handlers must have unique owners
   */
  ObservableHelper(const ObservableHelper&) = delete;
  ObservableHelper& operator=(const ObservableHelper&) = delete;

  /**
   * Move constructor. Must properly transfer ownership of resources. Doing copies as these are "cheap" elements
   * @param in input object to be moved to this
   */
  ObservableHelper(ObservableHelper&&) noexcept;
  ObservableHelper& operator=(ObservableHelper&& in) noexcept;

  /**
   * Destructor closes hdf5 remaining resources
   */
  ~ObservableHelper();

  /**
   * set the shape of this observable
   * @param dims dimensions
   * @param first starting index
   */
  void set_dimensions(std::vector<int>& dims, int first);

  /**
   * open a h5 group of this observable
   * Create a group for an observable and dataspace
   */
  void open(hid_t grp_id);


  /** add named property to describe the collectable this helper class handles
   * @param p any intrinsic datatype including vector, basic containers
   * @param pname property
   */
  template<typename T>
  inline void addProperty(T& p, const std::string& pname)
  {
    HDFAttribIO<T> a(p);
    a.write(data_id, pname.c_str());
  }

  void addProperty(float& p, const std::string& pname);
  void addProperty(Tensor<float, OHMMS_DIM>& p, const std::string& pname);
  void addProperty(Matrix<float>& p, const std::string& pname);
  void addProperty(TinyVector<float, OHMMS_DIM>& p, const std::string& pname);
  void addProperty(std::vector<float>& p, const std::string& pname);
  void addProperty(std::vector<TinyVector<float, OHMMS_DIM>>& p, const std::string& pname);

  void write(const value_type* first_v, const value_type* first_vv);

  /**
   * Check if hdf5 handlers are valid after a call to ObservableHelper::open
   * true: opened, false: not opened
   */
  bool isOpened() const noexcept;

private:
  ///true: hdf5 handlers are valid via open, false: otherwise
  bool isopened = false;

  ///closes remaining hdf5 handlers in destructor if isopened = true
  void close();
};
} // namespace qmcplusplus
#endif
