#include "hipblas.h"
#include "AFQMC/Numerics/detail/HIP/hipblas_utils.h"


namespace hipblas {

  hipblasStatus_t rocBLASStatusToHIPStatusAFQMC(rocblas_status_ error)
  {
      switch(error)
      {
      case rocblas_status_success:
          return HIPBLAS_STATUS_SUCCESS;
      case rocblas_status_invalid_handle:
          return HIPBLAS_STATUS_NOT_INITIALIZED;
      case rocblas_status_not_implemented:
          return HIPBLAS_STATUS_NOT_SUPPORTED;
      case rocblas_status_invalid_pointer:
          return HIPBLAS_STATUS_INVALID_VALUE;
      case rocblas_status_invalid_size:
          return HIPBLAS_STATUS_INVALID_VALUE;
      case rocblas_status_memory_error:
          return HIPBLAS_STATUS_ALLOC_FAILED;
      case rocblas_status_internal_error:
          return HIPBLAS_STATUS_INTERNAL_ERROR;
      default:
          throw "Unimplemented status";
      }
  }
}
