#ifndef QMCPLUSPLUS_PW_BAND_INFO_H
#define QMCPLUSPLUS_PW_BAND_INFO_H

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace qmcplusplus::pw
{
/** Descriptor for one primitive-cell band used to construct a plane-wave SPO.
 *
 * The projection flags describe how a complex band is expanded when QMCPACK is
 * built for real-valued wave functions.
 */
struct BandInfo
{
  /// Eigenvalue used by energy ordering.
  double energy;
  /// Primitive-cell twist index in the orbital HDF file.
  int twist_index;
  /// Band index within the primitive-cell twist and spin channel.
  int band_index;
  /// Select the imaginary rather than real projection of a complex band.
  bool use_imaginary_part{false};
  /// Expand this complex band into independent real and imaginary SPOs.
  bool make_two_copies{false};
};

/** Compare bands using the bspline sort=1 ordering.
 *
 * Bands are ordered by energy. Energies within 1e-6 are ordered by twist
 * index, then by band index.
 */
inline bool energyOrder(const BandInfo& left, const BandInfo& right)
{
  if (left.energy < right.energy + 1e-6 && left.energy > right.energy - 1e-6)
  {
    if (left.twist_index == right.twist_index)
      return left.band_index < right.band_index;
    return left.twist_index < right.twist_index;
  }
  return left.energy < right.energy;
}

/** Compare bands using the bspline sort=2 ordering.
 *
 * Bands are ordered by band index, then by energy. Equal energies within 1e-6
 * are ordered by twist index.
 */
inline bool indexOrder(const BandInfo& left, const BandInfo& right)
{
  if (left.band_index == right.band_index)
  {
    if (left.energy < right.energy + 1e-6 && left.energy > right.energy - 1e-6)
      return left.twist_index < right.twist_index;
    return left.energy < right.energy;
  }
  return left.band_index < right.band_index;
}

/** Apply the requested bspline-compatible ordering to a band list.
 *
 * sort_mode=0 preserves insertion order, sort_mode=1 orders by energy, and
 * sort_mode=2 orders primarily by band index.
 *
 * @throws std::runtime_error if sort_mode is not 0, 1, or 2.
 */
inline void sortBands(std::vector<BandInfo>& bands, int sort_mode)
{
  if (sort_mode == 1)
    std::sort(bands.begin(), bands.end(), energyOrder);
  else if (sort_mode == 2)
    std::sort(bands.begin(), bands.end(), indexOrder);
  else if (sort_mode != 0)
    throw std::runtime_error("Plane-wave sort must be 0, 1, or 2");
}
} // namespace qmcplusplus::pw

#endif