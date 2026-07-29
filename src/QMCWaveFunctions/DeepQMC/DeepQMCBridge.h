//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DEEPQMCBRIDGE_H
#define QMCPLUSPLUS_DEEPQMCBRIDGE_H

#include <memory>
#include <string>
#include <vector>

#include "Configuration.h"

namespace qmcplusplus
{

/** Interface between QMCPACK WaveFunctionComponent code and DeepQMC inference.
 *
 *  The bridge is intentionally batch-first.  Implementations are expected to
 *  amortize Python/JAX/device dispatch latency across all walkers in a Crowd.
 *  Returned quantities use QMCPACK's native WaveFunctionComponent convention:
 *  log(psi), grad log(psi), and laplacian log(psi).
 */
class DeepQMCBridge : public QMCTraits
{
public:
  struct BatchResult
  {
    std::vector<RealType> log_values;
    std::vector<RealType> grad_log_values;
    std::vector<RealType> lap_log_values;
  };

  virtual ~DeepQMCBridge() = default;

  virtual BatchResult evaluateLogBatch(const std::vector<RealType>& ion_coords,
                                       const std::vector<RealType>& electron_coords,
                                       int mol_idx,
                                       int batch_size,
                                       int n_elec) const = 0;
};

/** Construct the real Python/JAX-backed DeepQMC inference bridge. */
std::unique_ptr<const DeepQMCBridge> makePythonDeepQMCBridge(std::string model_path, std::string python_module_path);

/** Placeholder useful for tests and clear disabled-path errors. */
std::unique_ptr<const DeepQMCBridge> makeUnavailableDeepQMCBridge(std::string reason);

} // namespace qmcplusplus

#endif
