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
#include <stdexcept>
#include <memory>
#include "type_traits/OptionalRef.hpp"

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
  bool size() const { return myPool.size(); }

  /// add [key, value] pair to map
  void add(const std::string& key, std::unique_ptr<T>&& value)
  {
    if (key.empty())
      throw std::runtime_error("Report bug! ObjectPool::add doesn't allow an empty string as the object name.");

    if (contains(key))
      throw std::runtime_error(
          "Report bug! The caller of ObjectPool::add should check key existence by calling contains() and issue "
          "message to users before starting expensive object construction!");

    myPool.emplace(key, std::move(value));
  }

protected:
  /** look up object by name
   * @param name object name to look up. It cannot be an empty string
   */
  OptionalRef<T> getObjectByName(const std::string& name) const
  {
    guardEmptyPool();

    if (name.empty())
      throw std::runtime_error(
          "Report bug! ObjectPool::getObjectByName doesn't allow an empty string as the object name.");

    if (auto pit(myPool.find(name)); pit != myPool.end())
      return *pit->second;
    else
      return std::nullopt;
  }

  /** get the one and only object in the pool
   *
   */
  OptionalRef<T> getTheOneAndOnly() const
  {
    guardEmptyPool();

    if (myPool.size() == 1)
      return *myPool.begin()->second;
    else
      return std::nullopt;
  }

  /// storage of WaveFunctionFactory
  Pool myPool;

private:
  void guardEmptyPool() const
  {
    if (myPool.empty())
      throw std::runtime_error(
          "Report bug! Empty pool! Check pool size before any look-up and issue proper user error!");
  }
};
} // namespace qmcplusplus
#endif
