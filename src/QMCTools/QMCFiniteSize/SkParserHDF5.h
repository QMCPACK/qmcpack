#ifndef SK_PARSER_HDF5_H
#define SK_PARSER_HDF5_H

#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include "hdf/hdf_archive.h"
#include <vector>

namespace qmcplusplus
{
/** Class to handle reading the S(k) directly from stat.h5 file
 *
 * This class parses the kgrid and fluctuation S(k) from skall estimator.
 * In the parse function, it performs a simple equilibration estimate and 
 * block average to create the S(k) to use for the corrections
 */
class SkParserHDF5 : public SkParserBase
{
public:
  void parse(const std::string& fname) override;

private:
  hdf_archive statfile;
};

} // namespace qmcplusplus

#endif
