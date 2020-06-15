#ifndef SK_PARSER_HDF5_H
#define SK_PARSER_HDF5_H

#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include "io/hdf_archive.h"
#include <vector>

namespace qmcplusplus
{
class SkParserHDF5 : public SkParserBase
{
public:
  SkParserHDF5();
  ~SkParserHDF5();
  void parse(const string& fname) override;

private:
  hdf_archive statfile;
};

} // namespace qmcplusplus

#endif
