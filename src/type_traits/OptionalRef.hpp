//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_OPTIONALREF_HPP
#define QMCPLUSPLUS_OPTIONALREF_HPP

#include <optional>

namespace qmcplusplus
{
/** An optional, non-owning reference.
 *
 * This class intentionally stores a pointer instead of using
 * std::optional<std::reference_wrapper<T>>.  In C++17, reference_wrapper
 * examines whether T has legacy callable typedefs.  Doing so for a
 * forward-declared T causes GCC 16's -Wsfinae-incomplete warning when T is
 * later defined.
 */
template<typename T>
class OptionalRef
{
public:
  class Reference
  {
  public:
    constexpr Reference() noexcept = default;
    constexpr Reference(T& object) noexcept : object_(&object) {}

    constexpr Reference& operator=(T& object) noexcept
    {
      object_ = &object;
      return *this;
    }

    constexpr T& get() const noexcept { return *object_; }
    constexpr operator T&() const noexcept { return get(); }

  private:
    T* object_ = nullptr;
  };

  constexpr OptionalRef() noexcept = default;
  constexpr OptionalRef(std::nullopt_t) noexcept {}
  constexpr OptionalRef(T& object) noexcept : reference_(object), has_value_(true) {}

  constexpr OptionalRef& operator=(std::nullopt_t) noexcept
  {
    reset();
    return *this;
  }

  constexpr bool has_value() const noexcept { return has_value_; }
  constexpr explicit operator bool() const noexcept { return has_value(); }

  constexpr Reference& value()
  {
    if (!has_value_)
      throw std::bad_optional_access();
    return reference_;
  }

  constexpr const Reference& value() const
  {
    if (!has_value_)
      throw std::bad_optional_access();
    return reference_;
  }

  constexpr Reference& operator*() { return value(); }
  constexpr const Reference& operator*() const { return value(); }
  constexpr Reference* operator->() { return &value(); }
  constexpr const Reference* operator->() const { return &value(); }

  constexpr void emplace(T& object) noexcept
  {
    reference_ = object;
    has_value_ = true;
  }

  constexpr void reset() noexcept { has_value_ = false; }

private:
  Reference reference_;
  bool has_value_ = false;
};

template<typename T>
constexpr auto makeOptionalRef(T& obj)
{
  return OptionalRef<T>(obj);
}
} // namespace qmcplusplus
#endif
