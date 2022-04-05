#ifndef QMCPLUSPLUS_DRIVERTRAITS_H
#define QMCPLUSPLUS_DRIVERTRAITS_H

namespace qmcplusplus
{
/*! enum for QMC Run Type */
enum class QMCRunType
{
  DUMMY, /*!< dummy */
  VMC,   /**< VMC type: vmc, vmc-ptcl, vmc-multiple, vmc-ptcl-multiple */
  CSVMC,
  DMC,      /**< DMC type: dmc, dmc-ptcl*/
  RMC,      /**< RMC type: rmc, rmc-ptcl */
  VMC_OPT,  /*!< Optimization with vmc blocks */
  LINEAR_OPTIMIZE,
  WF_TEST,
  VMC_BATCH,
  DMC_BATCH,
  LINEAR_OPTIMIZE_BATCH
};

/** enum to set the bit to determine the QMC mode 
 *
 *  Can't be a scoped enum, unsafe out in qmcplusplus scope
 */
enum QMCModeEnum
{
  UPDATE_MODE,    /**< bit for move: walker or pbyp */
  MULTIPLE_MODE,  /**< bit for multple configuration */
  SPACEWARP_MODE, /**< bit for space-warping */
  ALTERNATE_MODE, /**< bit for performing various analysis and weird qmc methods */
  GPU_MODE,       /**< bit to use GPU driver */
  QMC_MODE_MAX = 8
};

} // namespace qmcplusplus
#endif
