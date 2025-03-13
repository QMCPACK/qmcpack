//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
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
  void registerVector(ListenerVector<Real>&& listener_vector)
  {
    listener_vectors_.push_back(std::move(listener_vector));
  }
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

class MockPerParticleEstimatorCrowd
{
public:
  /** Return listener frunction that has captured an object data member.
   *  returning a lambda allows access to the listening object controlled but allows a great deal of flexibility
   *  in dealing with the vector report. In this case the receiver copies the reported data into a CrowdEnergyValues data boject.
   */
  auto getParticularListener(CrowdEnergyValues<Real>& crowd_energy_values)
  {
    return [&crowd_energy_values](const int walker_index, const std::string& name, const Vector<Real>& inputV) {
      auto& local_values = crowd_energy_values;
      if (walker_index >= local_values[name].size())
        local_values[name].resize(walker_index + 1);
      local_values[name][walker_index] = inputV;
    };
  }
  ListenerVector<Real> makeListener() { return {"combiner", getParticularListener(crowd_energy_values_)}; }
  // For purposes of testing this is public.
  CrowdEnergyValues<Real> crowd_energy_values_;
};

class AnotherMockQMCHamiltonianAndReporter
{
private:
  std::vector<ListenerVector<Real>> listener_vectors_;
  const std::string name_{"PotentialA"};
  const std::string name2_{"PotentialB"};

public:
  /** why move or not move */
  void registerVector(ListenerVector<Real>&& listener_vector)
  {
    listener_vectors_.push_back(std::move(listener_vector));
  }
  void reportVector()
  {
    for (int walk_id : {0, 1, 2, 3})
    {
      Vector<Real> vec_part_(4);
      std::iota(vec_part_.begin(), vec_part_.end(), 0.5);
      for (auto& listener : listener_vectors_)
        listener.report(walk_id, name_, vec_part_);
      std::iota(vec_part_.begin(), vec_part_.end(), 0.75);
      for (auto& listener : listener_vectors_)
        listener.report(walk_id, name2_, vec_part_);
    }
  }
};

TEST_CASE("Listener::CrowdEnergyValues", "[hamiltonian]")
{
  AnotherMockQMCHamiltonianAndReporter another_mock_ham_reporter;
  MockPerParticleEstimatorCrowd mock_estimator;

  another_mock_ham_reporter.registerVector(mock_estimator.makeListener());
  another_mock_ham_reporter.reportVector();

  CHECK(mock_estimator.crowd_energy_values_["PotentialA"][0][0] == 0.5);
  CHECK(mock_estimator.crowd_energy_values_["PotentialA"][0][1] == 1.5);
  CHECK(mock_estimator.crowd_energy_values_["PotentialA"][0][2] == 2.5);
  CHECK(mock_estimator.crowd_energy_values_["PotentialA"][0][3] == 3.5);
  CHECK(mock_estimator.crowd_energy_values_["PotentialB"][1][0] == 0.75);
  CHECK(mock_estimator.crowd_energy_values_["PotentialB"][1][3] == 3.75);

  std::vector<Vector<Real>> reduced_over_reporters;
  combinePerParticleEnergies(mock_estimator.crowd_energy_values_, reduced_over_reporters);
  CHECK(reduced_over_reporters[0][0] == 1.25);
  CHECK(reduced_over_reporters[0][3] == 7.25);
  CHECK(reduced_over_reporters[2][0] == 1.25);
  CHECK(reduced_over_reporters[2][3] == 7.25);
}


} // namespace testing
} // namespace qmcplusplus
