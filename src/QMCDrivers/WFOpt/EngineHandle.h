//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Leon Otis, leon_otis@berkeley.edu, UC Berkeley
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENGINE_HANDLE_HEADER
#define QMCPLUSPLUS_ENGINE_HANDLE_HEADER

#include "Containers/MinimalContainers/RecordArray.hpp"
#include "Concurrency/OpenMP.h"

#include "QMCDrivers/Optimizers/DescentEngine.h"

#ifdef HAVE_LMY_ENGINE
#include "formic/utils/lmyengine/engine.h"
#endif


namespace qmcplusplus
{
class EngineHandle
{
public:
  using Real          = QMCTraits::RealType;
  using Value         = QMCTraits::ValueType;
  using FullPrecReal  = QMCTraits::FullPrecRealType;
  using FullPrecValue = QMCTraits::FullPrecValueType;

  virtual ~EngineHandle() = default;
  /** Function for preparing derivative ratio vectors used by optimizer engines
   *
   *\param[in] num_params           Number of optimizable parameters
   */
  virtual void prepareSampling(int num_params) = 0;
  /** Function for passing derivative ratios to optimizer engines
   *
   * \param[in] energy_list         Vector of local energy values
   * \param[in] dlogpsi_array       Parameter derivatives of log psi
   * \param[in] dhpsioverpsi_array  Parameter derivatives of local energy
   * \param[in] ib                  Sample index
   *
   */
  virtual void takeSample(const std::vector<FullPrecReal>& energy_list,
                          const RecordArray<Value>& dlogpsi_array,
                          const RecordArray<Value>& dhpsioverpsi_array,
                          int ib) = 0;
  /** Function for having optimizer engines execute their sample_finish functions
   */
  virtual void finishSampling() = 0;
};

class NullEngineHandle : public EngineHandle
{
public:
  void prepareSampling(int num_params) override {}
  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int ib) override
  {}
  void finishSampling() override {}
};

class DescentEngineHandle : public EngineHandle
{
private:
  DescentEngine& engine_;
  std::vector<FullPrecValue> der_rat_samp;
  std::vector<FullPrecValue> le_der_samp;

public:
  DescentEngineHandle(DescentEngine& engine) : engine_(engine) {}

  //Retrieve der_rat_samp vector for testing
  const std::vector<FullPrecValue>& getVector() const { return der_rat_samp; }

  void prepareSampling(int num_params) override
  {
    engine_.prepareStorage(omp_get_max_threads(), num_params);

    der_rat_samp.resize(num_params + 1, 0.0);
    le_der_samp.resize(num_params + 1, 0.0);
  }

  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int ib) override
  {
    der_rat_samp[0] = 1.0;
    le_der_samp[0]  = energy_list[ib];

    int num_params = der_rat_samp.size() - 1;
    for (int j = 0; j < num_params; j++)
    {
      der_rat_samp[j + 1] = static_cast<FullPrecValue>(dlogpsi_array.getValue(j, ib));
      le_der_samp[j + 1]  = static_cast<FullPrecValue>(dhpsioverpsi_array.getValue(j, ib)) +
          le_der_samp[0] * static_cast<FullPrecValue>(dlogpsi_array.getValue(j, ib));
    }
    int ip = omp_get_thread_num();
    engine_.takeSample(ip, der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
  }

  void finishSampling() override { engine_.sample_finish(); }
};

class LMYEngineHandle : public EngineHandle
{
#ifdef HAVE_LMY_ENGINE
private:
  cqmc::engine::LMYEngine<Value>& lm_engine_;
  std::vector<FullPrecValue> der_rat_samp;
  std::vector<FullPrecValue> le_der_samp;

public:
  LMYEngineHandle(cqmc::engine::LMYEngine<Value>& lmyEngine) : lm_engine_(lmyEngine){};

  void prepareSampling(int num_params) override
  {
    der_rat_samp.resize(num_params + 1, 0.0);
    le_der_samp.resize(num_params + 1, 0.0);
  }
  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int ib) override
  {
    der_rat_samp[0] = 1.0;
    le_der_samp[0]  = energy_list[ib];

    int num_params = der_rat_samp.size() - 1;
    for (int j = 0; j < num_params; j++)
    {
      der_rat_samp[j + 1] = static_cast<FullPrecValue>(dlogpsi_array.getValue(j, ib));
      le_der_samp[j + 1]  = static_cast<FullPrecValue>(dhpsioverpsi_array.getValue(j, ib)) +
          le_der_samp[0] * static_cast<FullPrecValue>(dlogpsi_array.getValue(j, ib));
    }

    lm_engine_.take_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
  }
  void finishSampling() override { lm_engine_.sample_finish(); }
#endif
};


} // namespace qmcplusplus
#endif
