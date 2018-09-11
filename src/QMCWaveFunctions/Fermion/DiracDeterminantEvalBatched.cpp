
#include "QMCWaveFunctions/Fermion/DiracDeterminantEval.h"

namespace qmcplusplus
{

  DiracDeterminantEval<Batching::BATCHED>::DiracDeterminantEval() :
    UpdateJobList_d("DiracDeterminantBase::UpdateJobList_d"),
    srcList_d("DiracDeterminantBase::srcList_d"),
    destList_d("DiracDeterminantBase::destList_d"),
    AList_d("DiracDeterminantBase::AList_d"),
    AinvList_d("DiracDeterminantBase::AinvList_d"),
    newRowList_d("DiracDeterminantBase::newRowList_d"),
    AinvDeltaList_d("DiracDeterminantBase::AinvDeltaList_d"),
    AinvColkList_d("DiracDeterminantBase::AinvColkList_d"),
    gradLaplList_d("DiracDeterminantBase::gradLaplList_d"),
    newGradLaplList_d("DiracDeterminantBase::newGradLaplList_d"),
    AWorkList_d("DiracDeterminantBase::AWorkList_d"),
    AinvWorkList_d("DiracDeterminantBase::AinvWorkList_d"),
    PivotArray_d("DiracDeterminantBase::PivotArray_d"),
    infoArray_d("DiracDeterminantBase::infoArray_d"),
    GLList_d("DiracDeterminantBase::GLList_d"),
    ratio_d("DiracDeterminantBase::ratio_d"),
    gradLapl_d("DiracDeterminantBase::gradLapl_d"),
    iatList_d("DiracDeterminantBase::iatList_d"),
    NLrowBuffer_d("DiracDeterminantBase::NLrowBuffer_d"),
    SplineRowList_d("DiracDeterminantBase::SplineRowList_d"),
    RatioRowList_d("DiracDeterminantBase::RatioRowList_d"),
    NLposBuffer_d("DiracDeterminantBase::NLposBuffer_d"),
    NLAinvList_d("DiracDeterminantBase::NLAinvList_d"),
    NLnumRatioList_d("DiracDeterminantBase::NLnumRatioList_d"),
    NLelecList_d("DiracDeterminantBase::NLelecList_d"),
    NLratioList_d("DiracDeterminantBase::NLratioList_d")
  {
    for(int i = 0; i < 2; ++i)
      NLratios_d[i] = gpu::device_vector<CudaValueType>("DiracDeterminantBase::NLratios_d");
  }
  
  void DiracDeterminantEval<Batching::BATCHED>::reserve (PointerPool<gpu::device_vector<QMCT::CudaValueType> > &pool,
						     SPOSetSingle& Phi)
{
    RowStride = ((NumOrbitals + 31)/32) * 32;
    AOffset           = pool.reserve((size_t)    NumPtcls * RowStride);
    AinvOffset        = pool.reserve((size_t)    NumPtcls * RowStride);
    gradLaplOffset    = pool.reserve((size_t)4 * NumPtcls * RowStride);
    newRowOffset      = pool.reserve((size_t)1            * RowStride);
    AinvDeltaOffset   = pool.reserve((size_t)1            * RowStride);
    AinvColkOffset    = pool.reserve((size_t)1            * RowStride);
    newGradLaplOffset = pool.reserve((size_t)4            * RowStride);
    if (typeid(CudaRealType) == typeid(float))
    {
      AWorkOffset       = pool.reserve((size_t)2 * NumPtcls * RowStride);
      AinvWorkOffset    = pool.reserve((size_t)2 * NumPtcls * RowStride);
    }
    else if (typeid(CudaRealType) == typeid(double))
    {
      AWorkOffset       = pool.reserve((size_t)    NumPtcls * RowStride);
      AinvWorkOffset    = 0;                  // not needed for inversion
    }
    Phi->reserve(pool);
}

}
