//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include <variant>
#include <string_view>
#include "catch.hpp"
#include "InputNode.hpp"
#include "type_traits/template_types.hpp"
#include "type_traits/variant_help.hpp"

namespace qmcplusplus
{

struct InputA : public InputNode
{
  InputA() = default;
  InputA(const std::string_view in_name, double in_dA, double in_dB, double in_dC, const std::string_view in_prop)
      : name(in_name), dA(in_dA), dB(in_dB), dC(in_dC), property(in_prop)
  {}
  std::string name{"inputA"};
  double dA{1.1};
  double dB{1.2};
  double dC{1.3};
  std::string property{"some prop"};
};

struct InputB : public InputNode
{
  InputB() = default;
  InputB(const std::string& in_name, short in_sA) : name(in_name), sA(in_sA) {}
  std::string name{"inputB"};
  short sA{1};
};

using ChildInput = std::variant<RefW<InputA>, RefW<InputB>>;

TEST_CASE("InputNode_BasicUse", "[io]")
{
  UPtrVector<InputNode> input_ownership;
  std::vector<ChildInput> inputs;
  InputNode::append<InputA>(input_ownership, inputs, "my_input_a", 2.0, 3.0, 4.0, "my_prop");
  InputNode::append<InputB>(input_ownership, inputs, "my_input_b", 2);
  CHECK(has<InputA>(inputs[0]));
  CHECK(!has<InputB>(inputs[0]));
  CHECK(has<InputB>(inputs[1]));
}

} // namespace qmcplusplus
