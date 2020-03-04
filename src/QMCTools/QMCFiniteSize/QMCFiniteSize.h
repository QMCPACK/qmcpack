#ifndef QMC_FINITE_SIZE_H
#define QMC_FINITE_SIZE_H

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Particle/ParticleSetPool.h"
#include "LongRange/LRCoulombSingleton.h"
#include "einspline/bspline_structs.h"
#include "einspline/nubspline_structs.h"

namespace qmcplusplus
{

//class SkParserBase
class QMCFiniteSize: public QMCAppBase,QMCTraits
{
  public:
    typedef LRCoulombSingleton::LRHandlerType  LRHandlerType;
    typedef LRCoulombSingleton::GridType       GridType;
    typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
    typedef LRHandlerType::mRealType           mRealType;
    typedef SkParserBase::Grid_t               Grid_t;
    typedef QMCTraits::RealType                RealType;
    typedef QMCTraits::FullPrecRealType        FullPrecRealType;
    typedef QMCTraits::PosType                 PosType;
    QMCFiniteSize();
    QMCFiniteSize(SkParserBase* skparser_i);
    ~QMCFiniteSize(){};


    inline void setSkParser(SkParserBase* skparser_i){skparser=skparser_i;};
    bool validateXML();
    bool execute();

    void build_spherical_grid(IndexType mtheta, IndexType mphi);
    void getSkInfo(UBspline_3d_d* spline, vector<RealType>& symmatelem);
    UBspline_3d_d* getSkSpline(RealType limit=1.0);
    RealType sphericalAvgSk(UBspline_3d_d* spline, RealType k);

    RealType integrate_spline(NUBspline_1d_d* spline, RealType a, RealType b, IndexType N);
    NUBspline_1d_d* spline_clamped(vector<RealType>& grid,
                                   vector<RealType>& vals,
                                   RealType lVal,
                                   RealType rVal);
  private:
    RealType h; //this is for finite differencing.
    LRHandlerType* AA;
    GridType* myGrid;
    RadFunctorType* rVs;

    ParticleSet* P;

    SkParserBase* skparser;
    bool processPWH(xmlNodePtr cur);
    bool wfnPut(xmlNodePtr cur);

    void initBreakup();
    ParticleSetPool ptclPool;

    RealType myConst;
    RealType myRcut;

    Grid_t gridx;
    Grid_t gridy;
    Grid_t gridz;

    IndexType mtheta;
    IndexType mphi;
    vector<PosType> sphericalgrid;

};
}

#endif
