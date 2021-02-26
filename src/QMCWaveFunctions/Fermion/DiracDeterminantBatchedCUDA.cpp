
#include "DiracDeterminantBatched.h"

namespace qmcplusplus
{
  
template<>
void DiracDeterminantBatched<MatrixDelayedUpdateCUDA<QMCTraits::ValueType, QMCTraits::QTFull::ValueType>>::mw_invertPsiM(const RefVectorWithLeader<WaveFunctionComponent>& wfc_list,
                                                             const RefVector<const ValueMatrix_t>& logdetT_list)
{
  auto& wfc_leader = wfc_list.getCastedLeader<DiracDeterminantBatched<DET_ENGINE_TYPE>>();
  ScopedTimer inverse_timer(&wfc_leader.InverseTimer);
  const auto nw = wfc_list.size();

  RefVector<LogValueType> log_value_list;

  for (int iw = 0; iw < nw; iw++)
  {
    auto& det = wfc_list.getCastedElement<DiracDeterminantBatched<DET_ENGINE_TYPE>>(iw);
    log_value_list.push_back(det.LogValue);
  }
  wfc_leader.det_inverter_.mw_invert_transpose(*cuda_handles_,logdetT_list, log_value_list);
}

}
