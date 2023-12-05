#ifndef SK_PARSER_SCALARDAT_H
#define SK_PARSER_SCALARDAT_H

#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include <vector>

namespace qmcplusplus
{
/** Class to handle Sk parsing from scalar.dat
 * This reads a processed scalar.dat from energy.pl
 * format.  Data to be read must have 
 * rhok_e_e val err
 * rhok_e_r val err
 * rhok_e_i val err
 */
class SkParserScalarDat : public SkParserBase
{
public:
  enum data_layout
  {
    KX,
    KY,
    KZ,
    SK,
    SKERR
  };

  void parse(const std::string& fname) override;

private:
  void read_sk_file(const std::string& fname);
};

} // namespace qmcplusplus
#endif
