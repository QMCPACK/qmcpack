#ifndef SK_PARSER_BASE_H
#define SK_PARSER_BASE_H

#include "einspline/bspline.h" 
#include "einspline/bspline_structs.h"
#include "Configuration.h"
#include "Numerics/OneDimGridBase.h"
#include <vector>


namespace qmcplusplus
{
using namespace std;

class SkParserBase: public QMCTraits
{
  public:
    typedef LinearGrid<RealType> Grid_t;
    
    SkParserBase();
    ~SkParserBase(){};
    
    virtual bool parse(const string& fname)=0;
    
    void get_grid(Grid_t& xgrid, Grid_t& ygrid, Grid_t& zgrid);
    void get_sk(vector<RealType>& sk, vector<RealType>& skerr);
    
    vector<PosType> get_grid_raw(){return kgridraw;};
    vector<RealType> get_sk_raw(){return skraw;};
    vector<RealType> get_skerr_raw(){return skerr_raw;};
    
    void compute_sk();
    void set_grid(const vector<PosType>& gridpoints);
    void set_grid(const vector<TinyVector<int,OHMMS_DIM> >& gridpoints);

    void compute_grid();

    inline bool is_normalized(){return isNormalized;}
    inline bool has_grid(){return hasGrid;}
  protected:

    bool isParseSuccess;
    bool isGridComputed;
    bool isSkComputed;
    bool hasGrid;
    bool isNormalized;
   
    Grid_t xgrid;
    Grid_t ygrid;
    Grid_t zgrid;

    vector<RealType> skraw;
    vector<RealType> skerr_raw;
    vector<PosType> kgridraw;
     
    vector<RealType> sk;
    vector<RealType> skerr;
    vector<PosType> kgrid;

    
};

}
#endif
