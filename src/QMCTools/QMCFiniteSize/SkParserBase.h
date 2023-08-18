#ifndef SK_PARSER_BASE_H
#define SK_PARSER_BASE_H

#include "einspline/bspline.h"
#include "einspline/bspline_structs.h"
#include "Configuration.h"
#include "Numerics/OneDimGridBase.h"
#include <vector>


namespace qmcplusplus
{
/** Base class for Sk parser
 *
 * parse is the only pure virtual function that must be overridden
 * holds various information about S(k) and access functions used by
 * qmcfinitesize
 */
class SkParserBase : public QMCTraits
{
public:
  using Grid_t = LinearGrid<RealType>;

  SkParserBase();
  virtual ~SkParserBase() {}

  virtual void parse(const std::string& fname) = 0;


  void get_grid(Grid_t& xgrid, Grid_t& ygrid, Grid_t& zgrid);
  void get_sk(std::vector<RealType>& sk, std::vector<RealType>& skerr);

  std::vector<PosType> get_grid_raw() { return kgridraw; };
  std::vector<RealType> get_sk_raw() { return skraw; };
  std::vector<RealType> get_skerr_raw() { return skerr_raw; };

  void compute_sk();
  void set_grid(const std::vector<PosType>& gridpoints);
  void set_grid(const std::vector<TinyVector<int, OHMMS_DIM>>& gridpoints);

  void compute_grid();

  inline bool is_normalized() { return isNormalized; }
  inline bool has_grid() { return hasGrid; }

  void setName(std::string in_name) { skname = in_name; }

protected:
  bool isParseSuccess;
  bool isGridComputed;
  bool isSkComputed;
  bool hasGrid;
  bool isNormalized;

  Grid_t xgrid;
  Grid_t ygrid;
  Grid_t zgrid;

  std::vector<RealType> skraw;
  std::vector<RealType> skerr_raw;
  std::vector<PosType> kgridraw;

  std::vector<RealType> sk;
  std::vector<RealType> skerr;
  std::vector<PosType> kgrid;

  std::string skname;
};

} // namespace qmcplusplus
#endif
