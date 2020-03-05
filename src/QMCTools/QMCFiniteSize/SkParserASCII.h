#ifndef SK_PARSER_ASCII_H
#define SK_PARSER_ASCII_H

#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Configuration.h"
#include <vector>

// FORMAT OF S(K) ASCII:
// 1 k_x k_y k_z S(k) err
// 2 k0x k0y k0z S(k0) err
// 3 k1x k1y k1z S(k1) err
// ...
//
namespace qmcplusplus
{
class SkParserASCII : public SkParserBase
{
public:
  SkParserASCII(){};
  ~SkParserASCII(){};

  enum data_layout
  {
    KX,
    KY,
    KZ,
    SK,
    SKERR
  };

  bool parse(const string& fname);

private:
  vector<vector<RealType>> read_sk_file(const string& fname);
  vector<PosType> get_grid_from_data(vector<vector<RealType>>& data);
  vector<RealType> get_sk_from_data(vector<vector<RealType>>& data);
  vector<RealType> get_skerr_from_data(vector<vector<RealType>>& data);
};

} // namespace qmcplusplus
#endif
