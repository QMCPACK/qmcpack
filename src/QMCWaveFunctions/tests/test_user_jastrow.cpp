//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/Jastrow/UserFunctor.h"

namespace qmcplusplus
{
TEST_CASE("UserJastrowFunctor", "[wavefunction]")
{
  using RealType = OptimizableFunctorBase::real_type;

  Communicate* c = OHMMS::Controller;

  UserFunctor<RealType> uf("test_functor");

  // This may or may not need to be present
  uf.setCusp(1.0);

  // Adjust this based on the example XML block
  const char* xmltext = R"(<tmp>
       <var name="B" id="j_B">2.0</var>
        </tmp>)";

  Libxml2Document doc;
  bool okay = doc.parseFromString(xmltext);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  uf.put(root);

  // Adjust these
  CHECK(uf.A == Approx(1.0));
  CHECK(uf.B == Approx(2.0));

  RealType r   = 1.0;
  RealType val = uf.evaluate(r);

  RealType dudr;
  RealType d2udr2;
  RealType val_2 = uf.evaluate(r, dudr, d2udr2);

  CHECK(val == Approx(val_2));

  // Finite difference to verify the spatial derivatives

  RealType h           = 0.01;
  RealType r_plus_h    = r + h;
  RealType r_minus_h   = r - h;
  RealType val_plus_h  = uf.evaluate(r_plus_h);
  RealType val_minus_h = uf.evaluate(r_minus_h);

  RealType approx_dudr = (val_plus_h - val_minus_h) / (2 * h);
  CHECK(dudr == Approx(approx_dudr).epsilon(h));

  RealType approx_d2udr2 = (val_plus_h + val_minus_h - 2 * val) / (h * h);
  CHECK(d2udr2 == Approx(approx_d2udr2).epsilon(h));


  RealType dudr_3;
  RealType d2udr2_3;
  RealType d3udr3_3;
  RealType val_3 = uf.evaluate(r, dudr_3, d2udr2_3, d3udr3_3);

  CHECK(val == Approx(val_3));
  CHECK(dudr == Approx(dudr_3));
  CHECK(d2udr2 == Approx(d2udr2_3));


  // Adjust this based on the number of variational parameters
  const int nparam = 1;

  // Outer vector is over parameters
  // Inner (TinyVector) is the parameter derivative of the value, first, and second derivatives.
  std::vector<TinyVector<RealType, 3>> param_derivs(nparam);

  uf.evaluateDerivatives(r, param_derivs);

  optimize::VariableSet var_param;
  uf.checkInVariablesExclusive(var_param);
  var_param.resetIndex();
  REQUIRE(var_param.size_of_active() == nparam);


  for (int i = 0; i < nparam; i++)
  {
    std::string var_name = var_param.name(i);
    RealType old_param   = std::real(var_param[var_name]);
    var_param[var_name]  = old_param + h;

    uf.resetParametersExclusive(var_param);

    RealType dudr_h;
    RealType d2udr2_h;
    RealType val_h = uf.evaluate(r, dudr_h, d2udr2_h);

    RealType val_dp = (val_h - val) / h;
    CHECK(val_dp == Approx(param_derivs[i][0]).epsilon(h));

    RealType dudr_dp = (dudr_h - dudr) / h;
    CHECK(dudr_dp == Approx(param_derivs[i][1]).epsilon(h));

    RealType d2udr2_dp = (d2udr2_h - d2udr2) / h;
    CHECK(d2udr2_dp == Approx(param_derivs[i][2]).epsilon(h));

    var_param[var_name] = old_param;
    uf.resetParametersExclusive(var_param);
  }

  // Could do finite differences to verify the parameter derivatives
}
} // namespace qmcplusplus
