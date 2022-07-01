//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ENGINE_HANDLE_HEADER
#define QMCPLUSPLUS_ENGINE_HANDLE_HEADER  

#include "QMCDrivers/Optimizers/DescentEngine.h"

#ifdef HAVE_LMY_ENGINE
#include "formic/utils/lmyengine/engine.h"
#endif


namespace qmcplusplus
{
    using ValueType = QMCTraits::ValueType;
  class EngineHandle
  {
          public:
    virtual void prepareSampling(int num_params) = 0;
    virtual void takeSample(std::vector<ValueType> energy_list, RecordArray<ValueType> dlogpsi_array, RecordArray<ValueType> dhpsioverpsi_array,int ib) = 0;
    virtual void finishSampling() = 0;
  };

  class NullEngineHandle : public EngineHandle
  {
public:
    void prepareSampling(int num_params) override {}
    void takeSample(std::vector<ValueType> energy_list, RecordArray<ValueType> dlogpsi_array, RecordArray<ValueType> dhpsioverpsi_array, int ib) override {}
    void finishSampling() override {}
  };

  //template<typename REAL>
  class DescentEngineHandle : public EngineHandle
  {
      private:
          DescentEngine& engine_;
std::vector<ValueType> der_rat_samp;
          std::vector<ValueType> le_der_samp;

          //cqmc::engine::LMYEngine<ValueType> lm_engine_;
public:
    DescentEngineHandle(DescentEngine& engine) : engine_(engine) {}

    void prepareSampling(int num_params) override 
    {
engine_.prepareStorage(omp_get_max_threads(), num_params);

    der_rat_samp.resize(num_params+1,0.0);
    le_der_samp.resize(num_params+1,0.0);
    }

    void takeSample(std::vector<ValueType> energy_list, RecordArray<ValueType> dlogpsi_array, RecordArray<ValueType> dhpsioverpsi_array,int ib) override 
    {
       der_rat_samp[0] = 1.0;
      le_der_samp[0] = energy_list[ib];

      int num_params = der_rat_samp.size()-1;
      for(int j = 0;j < num_params; j++)
      {
        der_rat_samp[j+1] = std::real(dlogpsi_array.getValue(j, ib));
        le_der_samp[j+1] = std::real(dhpsioverpsi_array.getValue(j, ib)) +le_der_samp[0]*std::real(dlogpsi_array.getValue(j, ib));
      
      }
      int ip = omp_get_thread_num();
      engine_.takeSample(ip, der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
    }

    void finishSampling() override 
{
   engine_.sample_finish(); 
}
  };

//template<typename REAL>
class LMYEngineHandle : public EngineHandle
  {
#ifdef HAVE_LMY_ENGINE
      private:
          cqmc::engine::LMYEngine<ValueType>& lm_engine_;
std::vector<ValueType> der_rat_samp;
          std::vector<ValueType> le_der_samp;

public:
    LMYEngineHandle(cqmc::engine::LMYEngine<ValueType>& lmyEngine):  lm_engine_(lmyEngine){ };

    void prepareSampling(int num_params) override 
{ 

    der_rat_samp.resize(num_params+1,0.0);
    le_der_samp.resize(num_params+1,0.0);

}
    void takeSample(std::vector<ValueType> energy_list, RecordArray<ValueType> dlogpsi_array, RecordArray<ValueType> dhpsioverpsi_array,int ib) override 
{ 
       der_rat_samp[0] = 1.0;
      le_der_samp[0] = energy_list[ib];

      int num_params = der_rat_samp.size()-1;
      for(int j = 0;j < num_params; j++)
      {
        der_rat_samp[j+1] = std::real(dlogpsi_array.getValue(j, ib));
        le_der_samp[j+1] = std::real(dhpsioverpsi_array.getValue(j, ib)) +le_der_samp[0]*std::real(dlogpsi_array.getValue(j, ib));
      
      }

      lm_engine_.take_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, 1.0);
}
    void finishSampling() override 
{ 
    lm_engine_.sample_finish();
}
#endif

  };

/*
  void QMCCostFunction::checkConfigurations(EngineHandle& engine_handle)
  {
    engine_handle.prepareSampling();
    engine_handle.takeSample();
    engine_handle.finishSampling();
  }
  */

} // namespace qmcplusplus
#endif
