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
  virtual void prepareSampling(int num_params, int num_samples) = 0;
  /** Function for passing derivative ratios to optimizer engines
   *
   * \param[in] energy_list         Vector of local energy values
   * \param[in] dlogpsi_array       Parameter derivatives of log psi
   * \param[in] dhpsioverpsi_array  Parameter derivatives of local energy
   * \param[in] local_index         Crowd local index
   * \param[in] sample_index        Index of sample on a MPI rank
   *
   */
  virtual void takeSample(const std::vector<FullPrecReal>& energy_list,
                          const RecordArray<Value>& dlogpsi_array,
                          const RecordArray<Value>& dhpsioverpsi_array,
                          int base_sample_index) = 0;
  /** Function for having optimizer engines execute their sample_finish functions
   */
  virtual void finishSampling() = 0;
};

class NullEngineHandle : public EngineHandle
{
public:
  void prepareSampling(int num_params, int num_samples) override {}
  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int base_sample_index) override
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

  void prepareSampling(int num_params, int num_samples) override
  {
    //FIXME it should respect num_samples and avoid relying on threads.
    engine_.prepareStorage(omp_get_max_threads(), num_params);

    der_rat_samp.resize(num_params + 1, 0.0);
    le_der_samp.resize(num_params + 1, 0.0);
  }

  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int base_sample_index) override
  {
    const int current_batch_size = dlogpsi_array.getNumOfEntries();
    for (int local_index = 0; local_index < current_batch_size; local_index++)
    {
      der_rat_samp[0] = 1.0;
      le_der_samp[0]  = energy_list[local_index];

      int num_params = der_rat_samp.size() - 1;
      for (int j = 0; j < num_params; j++)
      {
        der_rat_samp[j + 1] = static_cast<FullPrecValue>(dlogpsi_array[local_index][j]);
        le_der_samp[j + 1]  = static_cast<FullPrecValue>(dhpsioverpsi_array[local_index][j]) +
            le_der_samp[0] * static_cast<FullPrecValue>(dlogpsi_array[local_index][j]);
      }
      //FIXME it should respect base_sample_index and avoid relying on threads.
      int ip = omp_get_thread_num();
      engine_.takeSample(ip, der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
    }
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

  void prepareSampling(int num_params, int num_samples) override
  {
    der_rat_samp.resize(num_params + 1, 0.0);
    le_der_samp.resize(num_params + 1, 0.0);
    if (lm_engine_.getStoringSamples())
      lm_engine_.setUpStorage(num_params, num_samples);
  }
  void takeSample(const std::vector<FullPrecReal>& energy_list,
                  const RecordArray<Value>& dlogpsi_array,
                  const RecordArray<Value>& dhpsioverpsi_array,
                  int base_sample_index) override
  {
    int current_batch_size = dlogpsi_array.getNumOfEntries();
    for (int local_index = 0; local_index < current_batch_size; local_index++)
    {
      const int sample_index = base_sample_index + local_index;
      der_rat_samp[0]        = 1.0;
      le_der_samp[0]         = energy_list[local_index];

      int num_params = der_rat_samp.size() - 1;
      for (int j = 0; j < num_params; j++)
      {
        der_rat_samp[j + 1] = static_cast<FullPrecValue>(dlogpsi_array[local_index][j]);
        le_der_samp[j + 1]  = static_cast<FullPrecValue>(dhpsioverpsi_array[local_index][j]) +
            le_der_samp[0] * static_cast<FullPrecValue>(dlogpsi_array[local_index][j]);
      }


      if (lm_engine_.getStoringSamples())
        lm_engine_.store_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0, sample_index);
      else
        lm_engine_.take_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
    }
  }
  void finishSampling() override
  {
    if (!lm_engine_.getStoringSamples())
      lm_engine_.sample_finish();
  }
#endif
};


} // namespace qmcplusplus
#endif
