#ifndef SK_PARSER_SCALARDAT_H
#define SK_PARSER_SCALARDAT_H

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

class SkParserScalarDat: public SkParserBase
{
  public:
    SkParserScalarDat(){};
    ~SkParserScalarDat(){};
    
    enum data_layout { KX, KY, KZ, SK, SKERR};
    
    bool parse(const string & fname);
    
  private:
    void read_sk_file(const string& fname);

};

}
#endif
