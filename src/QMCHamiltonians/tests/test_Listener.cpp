//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: DensityMatrices1b.h
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Listener.hpp"
#include <string>
#include <vector>
#include "Configuration.h"

namespace qmcplusplus
{

using QMCT = QMCTraits;
using Real = QMCT::RealType;

namespace testing
{

/** Mock class that collects ListnerVectors as QMCHamiltonian does
 *   and reports ListenerVectors Hamiltonian operators do when they report per particle values.
 */
class MockQMCHamiltonianAndReporter
{
private:
  std::vector<ListenerVector<Real>> listener_vectors_;
  const std::string name_{"Talker"};

public:
  /** why move or not move */
  void registerVector(ListenerVector<Real>&& listener_vector) { listener_vectors_.push_back(std::move(listener_vector)); }
  void reportVector()
  {
    Vector<Real> vec_part(4);
    std::iota(vec_part.begin(), vec_part.end(), 0);
    for (auto& listener : listener_vectors_)
      listener.report(0, name_, vec_part);
  }
};

class MockPerParticleEstimator
{
public:
  /** Return listener frunction that has captured an object data member.
   *  returning a lambda allows access to the listening object controlled but allows a great deal of flexibility
   *  in dealing with the vector report. In this case the receiver just copies the vector it is called with to local storage
   *  which the lambda has captured.
   */
  auto getParticularListener(Vector<Real>& local_vector)
  {
    return [&local_vector](const int walker_index, const std::string& name, const Vector<Real>& values) -> void {
      local_vector = values;
    };
  }
  ListenerVector<Real> makeListener() { return {"kinetic", getParticularListener(receiver_vector_)}; }

  // For purposes of testing this is public.
  Vector<Real> receiver_vector_;
};

TEST_CASE("ListenerVector", "[hamiltonian]")
{
  MockQMCHamiltonianAndReporter mock_ham_report;
  MockPerParticleEstimator mock_estimator;

  mock_ham_report.registerVector(mock_estimator.makeListener());

  mock_ham_report.reportVector();
  CHECK(mock_estimator.receiver_vector_[0] == 0);
  CHECK(mock_estimator.receiver_vector_[3] == 3);
}

} // namespace testing
} // namespace qmcplusplus
