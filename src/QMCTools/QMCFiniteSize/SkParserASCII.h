#ifndef SK_PARSER_ASCII_H
#define SK_PARSER_ASCII_H

#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include <vector>

namespace qmcplusplus
{
/** Class to handle parsing from ASCII file
 *
 * This assumes you have process your own S(k) file
 * and has file format of
 * # k_x k_y k_z S(k)  err
 *   k0x k0y k0z S(k0) err
 *   k1x k1y k1z S(k1) err
 *   k2x k2y k2z S(k2) err
 *   ...
 */
class SkParserASCII : public SkParserBase
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
  std::vector<std::vector<RealType>> read_sk_file(const std::string& fname);
  std::vector<PosType> get_grid_from_data(std::vector<std::vector<RealType>>& data);
  std::vector<RealType> get_sk_from_data(std::vector<std::vector<RealType>>& data);
  std::vector<RealType> get_skerr_from_data(std::vector<std::vector<RealType>>& data);
};

} // namespace qmcplusplus
#endif
