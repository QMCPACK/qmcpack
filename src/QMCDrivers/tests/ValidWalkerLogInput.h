//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALIDWALKERLOGINPUT_H
#define QMCPLUSPLUS_VALIDWALKERLOGINPUT_H

#include <array>
#include <string_view>
namespace qmcplusplus::testing
{
class WalkerLogInputSections
{
public:
  enum class valid
  {
    DEFAULT = 0
  };

  static std::string_view getXml(valid val) { return xml[static_cast<std::size_t>(val)]; }
  auto begin() { return xml.begin(); }
  auto end() { return xml.end(); }

private:
  static constexpr std::array<std::string_view, 1> xml{
      R"XML(
<walkerlogs/>
)XML"};
};

} // namespace qmcplusplus::testing

#endif
