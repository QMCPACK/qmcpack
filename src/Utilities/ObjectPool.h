//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_OBJECTPOOL_HPP
#define QMCPLUSPLUS_OBJECTPOOL_HPP

#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <memory>
#include "StlPrettyPrint.hpp"

namespace qmcplusplus
{

template<class T>
class ObjectPool
{
public:
  using Pool = std::map<std::string, const std::unique_ptr<T>>;

  /** check if an object with the requested name can be found
   * @param name object name to look up
   */
  bool contains(const std::string& name) const { return myPool.find(name) != myPool.end(); }

  bool empty() const { return myPool.empty(); }

  /// add [key, value] pair to map
  void add(const std::string& key, std::unique_ptr<T>&& value)
  {
    if (key.empty())
      throw std::runtime_error("Report bug! Empty object name not allowed in the pool");

    if (contains(key))
      throw std::runtime_error("Report bug! The caller should check key existence by calling contains() and issue "
                               "message to users before starting expensive object construction!");

    myPool.emplace(key, std::move(value));
  }

protected:
  ObjectPool(const std::string& object_type_name) : object_type_name_(object_type_name) {}

  /** look up object by name
   * @param name object name to look up
   * if name is empty and the pool contains one entry, return the only entry
   * if name is not empty and not found in the pool, throw error
   */
  const std::unique_ptr<T>& getObject(const std::string& name) const
  {
    if (myPool.empty())
      throw std::runtime_error(object_type_name_ + " pool is empty! Cannot find " + object_type_name_ + " named \"" +
                               name + "\"!");

    if (name.empty() && myPool.size() == 1)
      return myPool.begin()->second;

    if (auto pit(myPool.find(name)); pit != myPool.end())
      return pit->second;
    else
    {
      std::ostringstream msg;
      msg << object_type_name_ << " pool contains " << myPool << "." << std::endl
          << "Cannot find " << object_type_name_ << " named \"" + name + "\"!";
      throw std::runtime_error(msg.str());
    }
  }

  /// storage of WaveFunctionFactory
  Pool myPool;

private:
  /// object type name
  const std::string object_type_name_;
};
} // namespace qmcplusplus
#endif
